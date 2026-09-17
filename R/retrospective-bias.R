
library(here)
library(tidyverse)
library(ggsidekick)

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


esc_2024 <- read.csv(here("data/generated/2024/esc-data.csv"))
esc_2024$rr_run <- "2024"
esc_2025 <- read.csv(here("data/generated/2025/esc-data.csv"))
esc_2025$rr_run <- "2025"



esc_comb <- rbind(esc_2024,esc_2025)

# CU escapement plot ----
esc_comb <- as.data.frame(esc_comb) |>
  mutate(year = as.numeric(year))


esc_comb$CU_f <- factor(esc_comb$stock, levels = CU_order)

esc_comb |>
  left_join(CU_name_lookup, by="CU_f") |>
  ggplot(aes(x = year, y = mean/1000)) +
  geom_line(lwd = 1.1, aes(col=rr_run)) +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~CU_pretty, ncol=3, scales = "free_y") +
  theme_sleek() +
  theme(strip.text = element_text(size=10))+
  scale_color_brewer(palette = "Dark2")+
  labs(color = "RR year")
ggsave("plots/RR-compare.jpeg", width = 9, height=6,units="in", dpi=600)
