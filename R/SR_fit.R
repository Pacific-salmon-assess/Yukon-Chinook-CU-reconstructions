# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
library(tidyverse)
library(rstan)

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

rm(harvest, ages)

set.seed(2)

# fit AR1 and time varying productivity (TVA) models--------------------------------------
  for(i in unique(sp_har$CU)){

    sp_har1 <- filter(sp_har, CU == i)

    stan.data <- list("nyrs" = nyrs,
                      "a_min" = a_min,
                      "a_max" = a_max,
                      "A" = A,
                      "nRyrs" = nRyrs,
                      "A_obs" = A_obs,
                      "S_obs" = sp_har1$spwn,
                      "H_obs" = sp_har1$harv,
                      "S_cv" = sp_har1$spwn_cv,
                      "H_cv" = sp_har1$harv_cv,
                      "Smax_p" = 0.75*max(sp_har1$spwn), #data for priors in semi_inform models, can tinker based on what assumed Smax is
                      "Smax_p_sig" = 0.75*max(sp_har1$spwn))

    AR1.fit <- stan(file = here("R/stan/SS-SR_AR1.stan"),
                    data = stan.data,
                    cores = 4,
                    seed = 2,
                    iter = 4000)

    saveRDS(AR1.fit, here(paste("data/generated/",lastyear,"/model_fits/AR1/",i,"_AR1.rds",sep="")))

    TV.fit <- stan(file = here("R/stan/SS-SR_TVA.stan"),
                   data = stan.data,
                  cores = 4,
                  seed = 2,
                  iter = 4000,
                  control = list(adapt_delta = 0.999,
                                 max_treedepth = 20))

    saveRDS(AR1.fit, here(paste("data/generated/",lastyear,"/model_fits/TVA/",i,"_TVA.rds",sep="")))
  }



