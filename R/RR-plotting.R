# run-reconstruction-plots.R

# fit of model to daily border passage ----
plotFitI(rpt, folder=here(paste("plots/",lastyear ,sep="")))

# fit of model to annual border passage ----
border <- read.csv(here("data/border-passage.csv")) |>
  filter(year>1984)

par( mar=c(5,5,0,0) )

#t <- !is.na(rpt$I_t)
t <- rep(TRUE,rpt$nT)
yr  <- rpt$years[t]
I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s[1])
E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
sonarN_t <- colSums(rpt$E_dtg[ ,t,1])*1e-3
ymax <- max(I_t,E_t,sonarN_t,na.rm=TRUE)
fw_t <- 1e-3*colSums(rpt$E_dtg[,,2])/exp(rpt$lnqE_tg[ ,2])

png(file=here(paste("plots/",lastyear ,"/CU-RR-fits.PNG",sep="")), width= 7, height = 5,units="in", res =700 )

plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="Year",
      ylab="Total border passage (1000s)", ylim=c(0,1.3*ymax) )
grid()
box()

if( is.finite(rpt$sdrpt[1,5]) )
{
  Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
  segments( x0=yr+0.2, y0=Ese$lCI*1e-3, y1=Ese$uCI*1e-3, col="grey70", lwd=2 )
  segments( x0=yr, y0=border$lwr*1e-3, y1=border$upr*1e-3, col="black", lwd=1.5 )

}

points( x=yr+0.2, y=E_t, pch=16, col="grey40" )
points( x=yr, y=I_t, pch=0, lwd=1.5 )
points( x=yr, y=sonarN_t, pch=16, lwd=1.5, col="red" )
points( x=yr, y=fw_t, pch=2, lwd=1.5, col="green" )

legend( x="bottomleft", bty="n",
        legend=c("CU border passage estimates","Aggregrate JTC border passage estimates","Sonar counts", "Fish wheel counts"),
        pch=c(NA,NA,NA,NA), lwd=c(1,1,NA,NA), col=c("grey70","black",NA,NA), lty=c(1,1,NA,NA), cex=0.75 )
legend( x="bottomleft", bty="n",
        legend=c("","","", ""),
        pch=c(16,0,16,2), lwd=c(1.5), col=c("grey40","black","red","green"), lty=c(0), cex=0.75 )

dev.off()

# Lookup table with CU names in order for plotting ----
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
CU_name_lookup <- data.frame(CU_f = factor(CU_order, levels=CU_order),
                             CU_pretty = factor(CU_prettynames, levels=CU_prettynames))


# trib time series ----
tribs.all <- read.csv(here("data/trib-spwn.csv")) |>
  mutate(
    estimate = case_when(
      system == "whitehorse" ~ estimate*(1-hatch_contrib),
      .default = estimate)) |>
  unite(tributary, c("system_alt", "type"), sep = " ") |>
  select(!hatch_contrib)


dat_text <- tribs.all |>
  group_by(tributary) |>
  slice_head() |>
  select(tributary, CU)


trib_order <- c("Porcupine sonar","Miner aerial","Klondike sonar","Chandindu weir","Tincup aerial" ,"Ross aerial","Pelly aerial",
                "Pelly sonar","Blind creek weir","Tachun foot", "Tachun weir","Little Salmon aerial","Big Salmon aerial",
                "Big Salmon sonar","Tahkini aerial","Tahkini sonar","Whitehorse fishway","Michie foot",
                "Teslin sonar","Nisutlin aerial","Nisutlin sonar","Wolf aerial", "Morley aerial")
tribs.all$tribs_name_ord<- factor(tribs.all$tributary, levels = trib_order)

tribs.all |>
  mutate(tribs_name = gsub("salmon", "Salmon", gsub("_", "-", str_to_sentence(tributary)))) |>
  mutate(CU_f = factor(gsub("BigSalmonR", "Big.Salmon", CU), levels=c(CU_order, "Porcupine"))) |>
  arrange(CU_f) |>
  mutate(tribs_name = factor(tribs_name, levels=unique(tribs_name))) |>
  ggplot(aes(x = year, y = estimate/1000, col=CU_f)) +
  scale_color_viridis_d() +
  geom_line(lwd = 0.8) +
  labs(x="Year", y="Spawners (000s)", col="CU") +
  facet_wrap(~tribs_name, ncol=4, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_sleek() +
  theme(axis.title = element_text(size=12),
        strip.text = element_text(size=8))

ggsave(here(paste("plots/",lastyear ,"/trib-escape.png",sep="")), height = 700*2,
       width=1000*2, units="px", dpi= 240)

# CU escapement plot ----
esc <- as.data.frame(mssr_spwn_2) |>
  mutate(year = as.numeric(year))


esc$CU_f <- factor(esc$stock, levels = CU_order)

esc |>
  left_join(CU_name_lookup, by="CU_f") |>
  ggplot(aes(x = year, y = mean/1000)) +
  geom_ribbon(aes(ymin = lower/1000, ymax = upper/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(lwd = 1.1, col="grey30") +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~CU_pretty, ncol=3, scales = "free_y") +
  theme_sleek() +
  theme(strip.text = element_text(size=10))

ggsave(here(paste("plots/",lastyear ,"/cu-escape.png",sep="")), width=900*2, height=800*2, units="px",
       dpi=240)

# Fishwheel catchability  ----
RR_pars <- rpt[["sdrpt"]]

fw_catch <- as.data.frame(RR_pars) |>
  filter(par=="lnqE_tg") |>
  mutate(mid = exp(val),
         lwr = exp(lCI),
         upr = exp(uCI),
         year = seq(1984,2006)) |>
  filter(year > 1984,
         year < 2005) |>
  select(year,mid,lwr,upr)

ggplot(fw_catch, aes(x = year, y = mid)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0,position=position_dodge(0.9)) +
  theme_sleek() +
  labs(x = "Year", y = "Fish wheel catchability")
ggsave(here(paste("plots/",lastyear ,"/fishwheel-catchability.png",sep="")), height=4.25, width=8)
