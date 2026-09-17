library(here)
library(tidyverse)
library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment
library(ggpubr)
library(viridis)
library(gsl)

# read in data ---------------------------------------------------------------------------

lastyear <- "2025"

#common data to be read in ---------------------------------------------------------------
harvest <- read.csv(here(paste("data/generated/",lastyear,"/harvest-data.csv",sep=""))) |>
  dplyr::rename(stock = population,
                harv_cv = cv)

sp_har <- read.csv(here(paste("data/generated/",lastyear,"/esc-data.csv",sep=""))) |>
  dplyr::rename(spwn = mean,
                spwn_cv = cv) |>
  left_join(harvest, by = c("stock", "year")) |>
  dplyr::rename(CU = stock) |>
  mutate(N = spwn+harv) |>
  arrange(CU)

ages <- read.csv(here(paste("data/generated/",lastyear,"/run-age-comp.csv",sep=""))) |>
  filter(Year >= min(sp_har$year),
         Year <= max(sp_har$year))

A_obs <- ages |>
  select(a4:a7) |>
  as.matrix()

a_min <- 4
a_max <- 7
nyrs <- max(sp_har$year)-min(sp_har$year)+1 #number of calendar years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #total number of fully and partially observed recruitments from observed spawners


# benchmark functions ---
get_Smsy <- function(a, b){
  Smsy <- (1-lambert_W0(exp(1-a)))/b
  if(Smsy <0){Smsy <- 0.001} #dumb hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy){Sgen*a*exp(-b*Sgen)-Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

# Lookup table with CU names in order for plotting
CU_order <- c("NorthernYukonR.andtribs.",
              "Whiteandtribs.", "Stewart",
              "MiddleYukonR.andtribs.","Pelly",
              "Nordenskiold", "Big.Salmon",
              "UpperYukonR.","YukonR.Teslinheadwaters")
CU_prettynames <- c("Northern Yukon R. and tribs.",
                    "White and tribs.", "Stewart",
                    "Middle Yukon R. and tribs.","Pelly",
                    "Nordenskiold", "Big Salmon",
                    "Upper Yukon R.","Yukon R. Teslin Headwaters")

# model fits ---
AR1.fits <- lapply(list.files(here(paste("data/generated/",lastyear,"/model_fits/AR1",sep="")),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)


# process data and fits to make table of benchamrks and plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to
bench.posts <- NULL
par.posts <- NULL
SR.preds <- NULL
AR1.spwn <- NULL
AR1.harv <- NULL
AR1.resids <- NULL
brood.all <- NULL
brood.all.long <- NULL
a.yrs.all <- NULL

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs

  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])

  AR1.spwn <- rbind(AR1.spwn, bind_cols(t(apply(sub_pars$S, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))

  AR1.harv <- rbind(AR1.harv, bind_cols(t(apply(sub_pars$H, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))

  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)], t(spwn.quant), t(rec.quant))) |>
    round(2)
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")

  brood_t <- mutate(brood_t, CU = i)

  brood.all <- rbind(brood.all, brood_t)

  spwn.quant.long <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs)]
  rec.quant.long <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

  brood_t.long <- as.data.frame(cbind(sub_dat$year[1:nyrs],t(spwn.quant.long), rbind(t(rec.quant.long),matrix(NA,4,3)))) |>
    round(2)
  colnames(brood_t.long) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")

  brood_t.long <- mutate(brood_t.long, CU = i)

  brood.all.long <- rbind(brood.all.long, brood_t.long)

  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  SR.pred <- matrix(NA,length(spw), length(sub_pars$lnalpha))
  bench <- matrix(NA,length(sub_pars$lnalpha),8,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("Sgen", "Smsy", "Umsy", "Seq", "Smsr", "S.recent","Smsr.20","Smsr.40")))

  par <- matrix(NA,length(sub_pars$lnalpha),3,
                dimnames = list(seq(1:length(sub_pars$lnalpha)), c("sample","ln_a","beta")))

  # get benchmarks & pars (AR1 model)------------------------------------------------------------------

  for(j in 1:length(sub_pars$lnalpha)){
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR.pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))

    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- 1/b #S_msr
    bench[j,6] <- exp(mean(log(sub_pars$S[j, (nyrs-5):nyrs]))) #S recent - mean spawners in last generation
    bench[j,7] <- (1/b)*0.2
    bench[j,8] <- (1/b)*0.4

    par[j,1] <- j
    par[j,2] <- ln_a
    par[j,3] <- b
  }

  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)

  SR.preds <- rbind(SR.preds, SR.pred)

  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))

  par.posts <- rbind(par.posts, as.data.frame(par) |> mutate(CU = i))

  bench.quant <- apply(bench[,1:8], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()

  mean <- apply(bench[,1:8],2,mean, na.rm=T) #get means of each

  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(CU = i) |>
    relocate('50%', 1)

  #other pars to report
  alpha <- quantile(exp(sub_pars$lnalpha), probs = c(.1, .5, .9))
  beta <- quantile(sub_pars$beta, probs = c(.1, .5, .9))
  sigma <- quantile(sub_pars$sigma_R, probs = c(.1, .5, .9))
  phi <- quantile(sub_pars$phi, probs = c(.1, .5, .9))

  par.quants <- rbind(alpha, beta, sigma, phi)

  #make big table of bench and pars
  par.summary <- as.data.frame(rstan::summary(AR1.fits[[i]])$summary) |>
    select(mean, n_eff, Rhat)

  #summarise not other pars...
  par.summary <- filter(par.summary, row.names(par.summary) %in% c('lnalpha', 'beta',
                                                                   'sigma_R', 'phi')) |>
    mutate(CU = i)
  par.summary[1,1] <- exp(par.summary[1,1]) #exp ln_alpha

  pars <- cbind(par.quants, par.summary)

  sub.bench.par.table <- bind_rows(sub_benchmarks, pars) |>
    mutate(n_eff = round(n_eff, 0),
           Rhat = round(Rhat, 4))

  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table))

  bench.par.table <- bind_rows(bench.par.table, sub.bench.par.table)

  # then residuals---
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs]

  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) |>
    mutate(CU = i)
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")

  AR1.resids <- rbind(AR1.resids, resids)


}  # End data wrangling loop by CU


colnames(SR.preds) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr", "CU")
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(AR1.spwn) <- c("S.25", "S.50", "S.75", "year", "CU")
colnames(AR1.harv) <- c("H.25", "H.50", "H.75", "year", "CU")


SR.preds$CU_f <- factor(SR.preds$CU, levels = CU_order)
AR1.spwn$CU_f <- factor(AR1.spwn$CU, levels = CU_order)
AR1.resids$CU_f <- factor(AR1.resids$CU, levels = CU_order)
AR1.harv$CU_f <- factor(AR1.harv$CU, levels = CU_order)
brood.all$CU_f <- factor(brood.all$CU, levels = CU_order)
brood.all.long$CU_f <- factor(brood.all.long$CU, levels = CU_order)
esc$CU_f <- factor(esc$stock, levels = CU_order)


# write important tables to repo ---------------------------------------------------------
bench.par.table.out <- bench.par.table |> # AR1 SR model pars and benchmarks
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:7, ~round(.,5)) |>
  arrange(bench.par, CU)

write.csv(bench.par.table.out,here(paste("data/generated/",lastyear,"/bench_par_table.csv",sep="")),row.names = F)

write.csv(brood.all.long, here(paste("data/generated/",lastyear,"/brood_table_long.csv",sep="")),
          row.names = FALSE)

