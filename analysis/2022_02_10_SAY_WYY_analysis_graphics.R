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

simdate <- "2020_05_11_SAY_WYY"


# ======================= #
# ---- Data processing ---- 
# ======================= #

BG <- read.table(paste(simdate, "_geno_freqs_0_9999.txt", sep = ""), T) %>% as_tibble()

BF <- BG %>% rename(FI = F_IE) %>% mutate(FF = F_SE + F_IU + F_SU) %>%
  select(sex, copy, deme, FI, FF, muActivAM, SAY, WYY)
BY <- BG %>% rename(YM = Y_E, X = Y_U) %>% 
  select(sex, copy, deme, YM, X, muActivAM, SAY, WYY)
BA <- BG %>% rename(AM = A_E, A = A_U) %>% 
  select(sex, copy, deme, AM, A, muActivAM, SAY, WYY)


# ======================= #
# ---- Data analysis ---- 
# ======================= #

# Fit GAMs.

# FI at T=0 and T=1. 
GFI0 <- gam(cbind(FI, FF) ~ t2(muActivAM, SAY, WYY,  bs = "tp",
                               k = 4, m = 2),
            data = BF %>% filter(sex == "Females", copy == 1, deme == 0), 
            family = "binomial", full = T, method = "REML")

GFI1 <- gam(cbind(FI, FF) ~ t2(muActivAM, SAY, WYY, bs = "tp",
                               k = 4, m = 2),
            data = BF %>% filter(sex == "Females", copy == 1, deme == 10), 
            family = "binomial", full = T, method = "REML")

# YM, dito.
GYM0 <- gam(cbind(YM, X) ~ t2(muActivAM, SAY, WYY, bs = "tp",
                              k = 4, m = 2),
            data = BY %>% filter(sex == "Males", copy == 2, deme == 0),
            family = "binomial", full = T, method = "REML")
GYM1 <- gam(cbind(YM, X) ~ t2(muActivAM, SAY, WYY, bs = "tp",
                              k = 4, m = 2),
            data = BY %>% filter(sex == "Males", copy == 2, deme == 10),
            family = "binomial", full = T, method = "REML")

# AM, dito.
GAM0 <- gam(cbind(AM, A) ~ t2(muActivAM, SAY, WYY, bs = "tp",
                              k = 4, m = 2),
            data = BA %>% filter(sex == "Males", copy == 2, deme == 0),
            family = "binomial", full = T, method = "REML")
GAM1 <- gam(cbind(AM, A) ~ t2(muActivAM, SAY, WYY, bs = "tp",
                              k = 4, m = 2),
            data = BA %>% filter(sex == "Males", copy == 2, deme == 10),
            family = "binomial", full = T, method = "REML")

# Save GAMs 
saveRDS(object = GFI0, file = paste(simdate, "_FI_GAM0.rds", sep = ""))
saveRDS(object = GFI1, file = paste(simdate, "_FI_GAM1.rds", sep = ""))
saveRDS(object = GYM0, file = paste(simdate, "_YM_GAM0.rds", sep = ""))
saveRDS(object = GYM1, file = paste(simdate, "_YM_GAM1.rds", sep = ""))
saveRDS(object = GAM0, file = paste(simdate, "_AM_GAM0.rds", sep = ""))
saveRDS(object = GAM1, file = paste(simdate, "_AM_GAM1.rds", sep = ""))

# Load GAMs 
GFI0 <- readRDS(paste(simdate, "_FI_GAM0.rds", sep = ""))
GFI1 <- readRDS(paste(simdate, "_FI_GAM1.rds", sep = ""))
GYM0 <- readRDS(paste(simdate, "_YM_GAM0.rds", sep = ""))
GYM1 <- readRDS(paste(simdate, "_YM_GAM1.rds", sep = ""))
GAM0 <- readRDS(paste(simdate, "_AM_GAM0.rds", sep = ""))
GAM1 <- readRDS(paste(simdate, "_AM_GAM1.rds", sep = ""))

# Check GAMs
gam.check(GFI0)
gam.check(GFI1)
gam.check(GYM0)
gam.check(GYM1)
gam.check(GAM0)
gam.check(GAM1)

# ========================= #
# ---- Generate graphics ---- 
# ========================= #

# Generate dataframe with predicted frequencies.
dp <- expand.grid(muActivAM = c(0.001, 0.003, 0.005, 0.01, 0.02, 0.03), 
                  SAY = seq(0, 0.1, length.out = 300), 
                  WYY = seq(0, 1, length.out = 300))

# Predict frequencies.
dp <- dp %>% transform(FI0 = predict(GFI0, newdata = dp, type = "response"),
                       FI1 = predict(GFI1, newdata = dp, type = "response"),
                       YM0 = predict(GYM0, newdata = dp, type = "response"),
                       YM1 = predict(GYM1, newdata = dp, type = "response"),
                       AM0 = predict(GAM0, newdata = dp, type = "response"),
                       AM1 = predict(GAM1, newdata = dp, type = "response"))

write.table(file = paste(simdate, "_predictions.txt", sep = ""), col.names = T, row.names = F, x = dp)
dp <- read.table(paste(simdate, "_predictions.txt", sep = ""), T)
# Convert into more usable formatting.
dpp <- dp %>% gather(FI0, FI1, YM0, YM1, AM0, AM1, key = "Allele", value =  "Frequency") %>%
  mutate(Locus = str_sub(Allele, start = 0, end = 2),
         Temperature = str_sub(Allele, start = 3, end = 3)) %>%
  mutate(Locus = factor(Locus, c("FI", "YM", "AM")))

levels(dpp$Locus) <- c("FI", "MY", "MA")


dpp <- dpp %>% as_tibble() %>% mutate(fL = str_c("italic(",str_sub(Locus, 1, 1), ")[", str_sub(Locus, 2,2),"]"),
                                      fT = str_c("italic(T)==", Temperature),
                                      fM = str_c("italic(µ)[D]==", muActivAM)) %>% # substitute Î¼ with mu symbol (ALT + 230) if necessary
  mutate(fL = factor(fL, c("italic(F)[I]", "italic(M)[Y]", "italic(M)[A]")))

dpp2 <- dpp %>% mutate(nT = as.integer(Temperature),
                       fL = factor(fL, c("italic(F)[I]", "italic(M)[Y]", "italic(M)[A]")))
dpp2$fTT <- str_c(c("North", "South")[dpp2$nT + 1], "~(", dpp$fT, ")")


# Define plot
plot <- dpp2 %>% ggplot(aes(SAY, WYY, fill = Frequency)) + 
  geom_tile() + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, limits = c(0,1)) + 
  facet_grid(fM ~ fL * fTT, labeller = label_parsed) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.1, 0.02),
                     labels = seq(0,10,2)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = bquote("Sexually antagonistic effect ("*italic(s)[a]%*%10^2*")"),
       y = bquote("YY survival ("*italic(s)[YY]*")"),
       fill = "Allele frequency") + 
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

plot
# Save plot
pdf(file = paste(simdate, "_allele_frequencies_v2022_02_10.pdf", sep = ""), height = 7, width = 7)
plot
dev.off()

png(file = paste(simdate, "_allele_frequencies_v2022_02_10.png", sep = ""), height = 7, width = 7, units = "in", res = 1000)
plot
dev.off()


# ========================== #
# ---- Musca-like systems ---- 
# ========================== #

# Convert data to frequencies, extract frequencies at T=0 and T=1.
BM <- BF %>% full_join(BY) %>% full_join(BA) %>%
  mutate(PF = FI / (FI + FF),
         PY = YM / (YM + X),
         PA = AM / (AM + A)) %>% 
  filter(deme %in% c(min(BF$deme), max(BF$deme))) %>% 
  group_by(sex, copy, muActivAM, SAY, WYY) %>%
  summarise(PF0 = PF[1], PF1 = PF[2],
            PY0 = PY[1], PY1 = PY[2],
            PA0 = PA[1], PA1 = PA[2])

# Subset datasets by locus and corresponding sex and maternal/paternal allele
# FI uses maternally-inherited allele in females. YM and AM use paternally-inherited allele in males
BMF <- BM %>% filter(sex == "Females", copy == 1) %>% 
  group_by(muActivAM, SAY, WYY) %>% 
  summarise(PF0 = PF0, PF1 = PF1)

BMY <- BM %>% filter(sex == "Males", copy == 2) %>% 
  group_by(muActivAM, SAY, WYY) %>% 
  summarise(PY0 = PY0, PY1 = PY1)

BMA <- BM %>% filter(sex == "Males", copy == 2) %>% 
  group_by(muActivAM, SAY, WYY) %>% 
  summarise(PA0 = PA0, PA1 = PA1)

# Join subsets and define Musca-like system. YM goes from major to minor allele, vice versa for AM 
# and FI. Musca-like system is defined as having all three of these gradients.
BM3 <- BMF %>% full_join(BMY) %>% full_join(BMA) %>% 
  mutate(Musca = PF0 < 1/2 & PF1 > 1/2 & 
           PY0 > 1/2 & PY1 < 1/2 &
           PA0 < 1/2 & PA1 > 1/2)




# Fit GAM to Musca-like system data.
MG <- gam(Musca ~ t2(muActivAM, SAY, WYY,  bs = "tp",
                     m = 2),
           data = BM3, 
           family = "binomial", full = T, method = "REML")
gam.check(MG)


# Fit GAM to Musca-like system data.
MG <- gam(Musca ~ t2(muActivAM, SAY, WYY,  bs = "tp",
                     m = 2),
          data = BM3, 
          family = "binomial", full = T, method = "REML")
gam.check(MG)
saveRDS(MG, file = paste0(simdate, "_MuscaGAM.rds"))
MG <- readRDS(paste0(simdate, "_MuscaGAM.rds"))

# Generate predictions.
dpm <- expand.grid(muActivAM = c(0.001, 0.003, 0.005, 0.01, 0.02, 0.03), 
                   SAY = seq(0, 0.1, length.out = 300), 
                   WYY = seq(0, 1, length.out = 300))
dpm <- dpm %>% transform(M3 = predict(MG, newdata = dpm, type = "response"))
dpm <- dpm %>% mutate(fM = str_c("italic(µ)[D]==", muActivAM)) # substitute Î¼ with mu symbol (ALT + 230) if necessary

# Generate plot.
plot_combined <- dpm %>% 
  gather(M3, key = "Type", value = "Musca") %>% 
  ggplot(aes(SAY, WYY, fill = Musca)) + geom_tile() + 
  geom_tile() + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, limits = c(0,1)) + 
  facet_wrap(~fM, ncol = length(unique(dpm$muActivAM)) / 2, labeller = label_parsed) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.1, 0.02),
                     labels = seq(0,10,2)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = bquote("Sexually antagonistic effect ("*italic(s)[a]%*%10^2*")"),
       y = bquote("YY survival ("*italic(s)[YY]*")"),
       fill = "Frequency housefly-\nlike system") + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(0,0,0,0, "lines"),
        panel.spacing = unit(0.75, "lines"), 
        text = element_text(size = 9),
        legend.position = "bottom",
        legend.box.margin = margin(-10,0,0,0),
        legend.margin = margin(0,0,0,0)) + 
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = 0.5))

plot_combined
# Save plots to output files.
pdf(file = paste(simdate, "_Musca_evolution_v2022_02_10.pdf", sep = ""), height = 3.2, width = 4)
plot_combined
dev.off()

png(file = paste(simdate, "_Musca_evolution_v2022_02_10.png", sep = ""), height = 3.2, width = 4, units = "in", res = 1000)
plot_combined
dev.off()

plot_combined2 <- dpm %>% 
  gather(M3, key = "Type", value = "Musca") %>% 
  ggplot(aes(SAY, WYY, fill = Musca)) + geom_tile() + 
  geom_tile() + 
  scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.95, limits = c(0,1)) + 
  facet_wrap(~fM, ncol = length(unique(dpm$muActivAM)), labeller = label_parsed) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.1, 0.02),
                     labels = seq(0,10,2)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = bquote("Sexually antagonistic effect ("*italic(s[M])%*%10^2*")"),
       y = bquote("YY survival ("*italic(S[YY])*")"),
       fill = "Frequency housefly-\nlike system") + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(0,0,0,0, "lines"),
        panel.spacing = unit(0.75, "lines"), 
        text = element_text(size = 9),
        legend.position = "bottom",
        legend.box.margin = margin(-10,0,0,0),
        legend.margin = margin(0,0,0,0)) + 
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = 0.5))

pdf(file = paste(simdate, "_Musca_evolution_v2021_06_15_wide.pdf", sep = ""), height = 2, width = 6)
plot_combined2
dev.off()

png(file = paste(simdate, "_Musca_evolution_v2021_06_15_wide.png", sep = ""), height = 2, width = 6, units = "in", res = 1000)
plot_combined2
dev.off()