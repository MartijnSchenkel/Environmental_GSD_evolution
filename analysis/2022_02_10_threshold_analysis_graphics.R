# ============================================== #
#         Graphics for temperature model         #
# ============================================== #

# Author: Martijn Schenkel
# Purpose: Produce graphs using geno produced by model for temperature and SD gene expression

# ==================== #
# Administrative setup #
# ==================== #
# Clear workspace and load required packages
rm(list = ls())
library(tidyverse)
library(viridis)
library(cowplot)
library(mgcv)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

simdate <- "2020_05_25_threshold"

# ======================= #
# ---- Data processing ---- 
# ======================= #

# Read table, filter into loci-specific datasets.
BG <- read.table(paste(simdate, "_geno_freqs_0_9999.txt", sep = ""), T) %>% as_tibble()

BF <- BG %>% rename(FI = F_IE) %>% 
  mutate(FF = F_SE + F_IU + F_SU) %>%
  select(sex, copy, deme, FI, FF, muActivAM, BetaT, ThresF) %>% 
  filter(sex == "Females", copy == 1)

# Raw data plot. First bin by muActivAM for facetting.
BF$TP <- cut(BF$muActivAM, breaks = c(0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05))

# Plot raw data.
BF %>% ggplot(aes(BetaT, ThresF, color = FI/(FF+FI))) + 
  geom_point() + 
  facet_wrap(~TP) + 
  scale_color_viridis()

# ===================== #
# ---- Data analysis ----
# ===================== #

# FD
GFI0 <- gam(cbind(FI, FF) ~ t2(muActivAM, BetaT, ThresF,
                               bs = "tp", k = 4, m = 2),
            data = BF %>% filter(deme == 0), 
            family = "binomial", full = F, method = "REML")

GFI1 <- gam(cbind(FI, FF) ~ t2(muActivAM, BetaT, ThresF, 
                               bs = "tp", k = 4, m = 2),
            data = BF %>% filter(deme == 10), 
            family = "binomial", full = F, method = "REML")

# Save GAMs
saveRDS(object = GFI0, file = paste(simdate, "_FI_GAM0.rds", sep = ""))
saveRDS(object = GFI1, file = paste(simdate, "_FI_GAM1.rds", sep = ""))

# Load GAMs
GFI0 <- readRDS(paste(simdate, "_FI_GAM0.rds", sep = ""))
GFI1 <- readRDS(paste(simdate, "_FI_GAM1.rds", sep = ""))

# Check GAMs
gam.check(GFI0)
gam.check(GFI1)

# ============================ #
# ---- Supplementary figure ---- 
# ============================ # 

# Generate data set with predictions.
dp <- expand.grid(muActivAM = c(0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03), 
                  BetaT = seq(0, 0.8, length.out = 300), 
                  ThresF = seq(0.8, 1.5, length.out = 300))

dp <- dp %>% 
  transform(FI0 = predict(GFI0, newdata = dp, type = "response"),
            FI1 = predict(GFI1, newdata = dp, type = "response"))

write.table(file = paste(simdate, "_predictions.txt", sep = ""), col.names = T, row.names = F, x = dp)
dp <- read.table(paste(simdate, "_predictions.txt", sep = ""), T)

# Process.
dpp <- dp %>% gather(FI0, FI1, key = "Allele", value =  "Frequency") %>%
  mutate(Temperature = str_sub(Allele, start = 3, end = 3),
         BetaT = 1 + BetaT,
         fT = paste0("italic(T)==", Temperature),
         fM = paste0("italic(µ)[D]==",muActivAM)) # Change ? to mu (ALT+230)



plot <- dpp %>% ggplot(aes(BetaT, ThresF, fill = Frequency)) + 
  geom_tile() + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, limits = c(0,1)) + 
  facet_grid(fM ~ fT, labeller = label_parsed) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = bquote("Temperature-dependent expression rate (1 +"~italic(beta)[T]*")"),
       y = bquote("Feminization threshold ("*italic(theta)[F]*")"),
       fill = bquote("Frequency of"~italic(F)^I)) +
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(0,0,0,0, "lines"),
        panel.spacing = unit(0.75, "lines"), 
        text = element_text(size = 9),
        legend.position = "bottom",
        legend.box.margin = margin(-10,-10,0,-10),
        legend.margin = margin(0,0,0,0)) + 
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = 0.5))

pdf(file = paste(simdate, "_allele_frequencies_v2021_06_15.pdf", sep = ""), height = 8, width = 6)
plot
dev.off()

png(file = paste(simdate, "_allele_frequencies_v2021_06_15.png", sep = ""), height = 8, width = 6, units = "in", res = 1000)
plot
dev.off()

# ==================== #
# ---- Paper figure ----
# ==================== # 

# This part of the script is used to generate Figure 2 from the manuscript.

# Filter relevant info, only need maternal copy in females.
BFF <- BF %>% filter(sex == "Females", copy == 1)

# Add unique identifier per simulation
BFF$n <- rep(seq(1, length(BFF$sex)/11), rep(11, length(BFF$sex)/11))

# Group sims, check if FI frequency exhibits gradient-like behaviour (i.e. goes from <1/2 to >1/2 frequency)
BFF2 <- BFF %>% group_by(n) %>% 
  summarise(gradient = max(FI/(FI+FF)) > 1/2 & min(FI/(FI+FF)) < 1/2) %>% 
  full_join(BFF) %>%
  filter(gradient == T)

# Set seed, this ensures same simulation is sampled all the time. Seed choice is arbitrary obviously.
set.seed(321)
BFS <- BFF2 %>% filter(n == sample(unique(BFF2$n), 1))

# Fish out relevant parameter values.
BFS[1,]$muActivAM
BFS[1,]$BetaT
BFS[1,]$ThresF
dpp$Locality <- c("North", "South")[1 + (dpp$Temperature==1)]
dpp <- dpp %>% mutate(fT = paste0(Locality,"~(italic(T)==", Temperature, ")"))




# Define part A of the figure.
PSA <- dpp %>% 
  filter(muActivAM == 0.001) %>% ggplot(aes(BetaT, ThresF, fill = Frequency)) + 
  geom_tile() + 
  geom_hline(yintercept = BFS[1,]$ThresF, linetype = 2, color = "white") + 
  geom_vline(xintercept = 1 + BFS[1,]$BetaT, linetype = 2, color = "white") + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, limits = c(0,1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  labs(x = bquote("Temperature-dependent expression rate (1 +"~italic(ß)*")"), # beta = ALT + 225 
       y = bquote("Feminization threshold ("*italic(theta)[F]*")"), 
       fill = bquote("Frequency of"~italic(F)[I])) + 
  facet_wrap(~fT, labeller = label_parsed) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(1/2,1/2,1,1/2, "lines"),
        panel.spacing.x = unit(0.75, "lines"), 
        text = element_text(size = 9),
        legend.position = "none",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,1,0, "lines")) + 
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = 0.5))

# Define part B.
PSB <- BFS %>% rename(F = FF) %>%
  gather(F, FI, key = "Allele", value = "Frequency") %>%
  mutate(Allele = factor(Allele, c("F", "FI")),
         Temperature = deme / max(deme)) %>% 
  ggplot(aes(deme, Frequency, fill = Allele)) + 
  geom_col(position = "fill", width = 1) + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, discrete = T, labels = c(bquote(italic(F)), bquote(italic(F)[I]))) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(min(BFS$deme), max(BFS$deme), length.out = 3),
                     labels = c("0", "1/2", "1")) + 
  labs(x = "Temperature") +
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(1/2,1/2,1,1/2, "lines"),
        text = element_text(size = 9),
        legend.position = "none",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,1,0, "lines")) + 
  guides(fill = guide_legend(title.position = "top",
                               title.hjust = 0.5))

# Extract legends.
PSAL <- get_legend(PSA + theme(legend.position = "bottom"))
PSBL <- get_legend(PSB + theme(legend.position = "bottom"))

# Combine into figure.
paper_figure <- plot_grid(PSA, PSB, PSAL, PSBL, ncol = 2, labels = c("A","B"), rel_widths = c(2,1), rel_heights = c(1,0.15))
paper_figure

# Write to output.
pdf(file = paste(simdate, "_paper_figure_v2022_02_10.pdf", sep = ""), width = 7, height = 3.5)
paper_figure
dev.off()

png(file = paste(simdate, "_paper_figure_v2022_02_10.png", sep = ""), width = 7, height = 3.5, units = "in", res = 1000)
paper_figure
dev.off()
