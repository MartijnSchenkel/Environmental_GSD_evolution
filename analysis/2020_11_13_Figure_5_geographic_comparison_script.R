# ============================================== #
#         Graphics for temperature model         #
# ============================================== #

# Author: Martijn Schenkel
# Purpose: Produce graphs using geno produced by model for temperature and SD gene expression

# ============================ #
# ---- Administrative setup ----
# ============================ #
# Clear workspace and load required packages
rm(list = ls())
library(tidyverse)
library(viridis)
library(cowplot)
library(mgcv)
library(maps)
library(mapproj)
library(readxl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Allele status types
isSensitive <- c("I", "S") # Insensitive / Sensitive
isExpressed <- c("U", "E") # Unexpressed / Expressed
isActive <- c("I", "A") # Inactive / Active

# Thresholds for classsification.
TEF <- 0.15 # Threshold below which an F allele is considered non-expressed
TSF <- 0.0 # Threshold below which an F allele is considered insensitive to M
TEM <- (2-1.15)/2 # Threshold below which an M allele is considered non-expressed

# ======================= #
# ---- Genotype data processing ----
# ======================= #

# Read data and parameter file.
dt <- cbind(read.table(dir("2020_06_12_HOTSEX/2020_06_12_HOTSEX/", full.names = T, pattern="geno\\.txt$")[1], header = T))
pt <- cbind(read.table(dir("2020_06_12_HOTSEX/2020_06_12_HOTSEX/", full.names = T, pattern="values\\.txt$")[1], header = T)) 

# Add relevant parameters.
dt <- dt %>% mutate(muActivAM = pt$value[46], BetaT = pt$value[49], ThresF = pt$value[29],
                    SAY = pt$value[17], WYY = pt$value[14],
                    deme = as.numeric(deme), sex = factor(sex))
dt$ID <- rep(seq(1, length(dt$sex)/10), rep(10, length(dt$sex)/10)) # Unique ID for per individual
dt$copy <- rep(rep(seq(1:2)), length(dt$sex)/2) # Maternal (1) and paternal (2) alleles
levels(dt$sex) <- c("Females", "Males", "Intersex")

# Categorize F alleles.
Fs <- dt %>% filter(locus == "F") %>% spread(key = variable, value = value) %>% 
  mutate(Insensitive = isSensitive[(Sensitivity > TSF)+1],
         Expressed = isExpressed[(Expression >= TEF)+1],
         Status = paste("F_", Insensitive, Expressed, sep = "")) %>%
  with(table(sex, copy, deme, muActivAM, BetaT, ThresF, SAY, WYY, factor(Status, levels = c("F_SE", "F_SU", "F_IE", "F_IU")))) %>%
  as.data.frame() %>% spread(key = Var9, value = Freq)

# Join three loci into 1 dataframe.
at <- Fs %>% as_tibble() %>%
  mutate(deme = as.numeric(deme) - 1)

# Observed frequencies map 

d <- read_xlsx("2020_06_15_Kozielska_et_al_data.xlsx", 1)
md <- map_data("world")

# Scale deme to observed latitudes, convert no. F_IE into prop. FI-bearing females
at <- at %>% mutate(lat = max(d$lat) - ((max(d$lat) - min(d$lat)) / max(at$deme)) * deme,
                    temperature = max(d$temperature) - ((max(d$temperature) - min(d$temperature)) / max(at$deme)) * (max(deme) - deme),
                    pFD = F_IE / (F_IE + F_SE + F_IU + F_SU))
d <- d %>% mutate(Type = "Observed")
pd <- at %>% filter(sex == "Females", copy == 1) %>% select(lat, pFD, temperature) %>% 
  mutate(Type = "Predicted") %>% full_join(d)

right <- pd %>% ggplot(aes(x = lat, y = pFD, color = Type)) +
  scale_x_continuous(limits = c(37, 52), breaks = seq(37, 52, 3), expand = c(0,0)) +
  geom_point(size = 2) + 
  coord_flip() +
  scale_color_viridis(discrete = T) + 
  labs(y = "Frequency", x = "Latitude")  + 
  guides(color = guide_legend(title.position = "top")) + 
  theme(legend.position = "none", legend.title.align = 0.5, panel.grid.minor = element_blank())

# Mutate M-factor data to proportion of total M
dp <- d %>% mutate(totM = SC + A1 + A2 + A3 + A4 + A4 + A5) %>%
  mutate(SC = SC / totM,
         A1 = A1 / totM,
         A2 = A2 / totM,
         A3 = A3 / totM,
         A4 = A4/ totM,
         A5 = A5 / totM)
dp$longdev <- c(0, 0, 0, 
                0, -0.3, -0.3, 
                0, 0, 0, 
                0, 0, 0,
                0, 0, 0)
dp$latdev <- c(-1, 0, 0, 
               0, 0, -1, 
               -1, 0, 0, 
               0, 0, 0,
               0, -1, 0)

cols <- viridis_pal(begin = 0.4, end = 0.9, option = "magma", direction = -1)(6)

# Plot F-factor data
F_plot <- ggplot(data = md, aes(x = long, y = lat)) + 
  coord_map(projection = "mercator", xlim = c(6,18), ylim = c(37, 52)) +
  geom_polygon(aes(group = group), fill = rgb(0.95, 0.95, 0.95), colour = "black") + 
  geom_point(data = dp, aes(x = long, y = lat), color = "black", size = 1.5) +
  geom_rect(data = dp, aes(ymin = lat + latdev, ymax = lat + latdev + (1 - pFD), 
                           xmin = long + longdev , xmax = long + longdev + 0.3), color = "black", fill =  cols[1]) + # F
  geom_rect(data = dp, aes(ymin = lat + latdev + (1 - pFD), ymax = lat + latdev + 1, 
                           xmin = long + longdev, xmax = long + longdev + 0.3), color = "black", fill = cols[6]) + # FD
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = rgb(0.8, 0.9, 1))) +
  scale_y_continuous(limits = c(35,55), expand = c(0,0), breaks = seq(37, 52, 3)) + 
  scale_x_continuous(limits = c(0,20), expand = c(0,0), breaks = seq(6, 18, 3)) + 
  labs(y = "Latitude", x = "Longitude")  + 
  theme(legend.position = "none", legend.title.align = 0.5)


# Generate plots to extract legends, not easily done with above plots because of geom_rect for bars.
df <- data.frame(a = c(1, 2), b = c(1,2), Allele = c("F", "FD"))

F_legend_plot <- ggplot(df, aes(a, b, fill = Allele)) + geom_col() + 
  scale_fill_manual(values = cols[c(1,6)], labels = c(bquote(italic(tra^phantom(0))), bquote(italic(tra^D)))) + 
  labs(fill = "Female type") + 
  guides(fill = guide_legend(title.position = "top", nrow = 1)) + 
  theme(legend.position = "bottom", legend.title.align = 0.5)

F_legend <- get_legend(F_legend_plot)
right_legend <- get_legend(right + theme(legend.position = "bottom"))

pdf(file = "2020_11_13_Figure_5.pdf", width = 8, height = 6)
plot_grid(F_plot, right, F_legend, right_legend, ncol = 2, nrow = 2, rel_heights = c(1, 0.1), labels = c("A", "B", "", ""))
dev.off()

png(file = "2020_11_13_Figure_5.png", width = 8, height = 6, res = 1000, units = "in")
plot_grid(F_plot, right, F_legend, right_legend, ncol = 2, nrow = 2, rel_heights = c(1, 0.1), labels = c("A", "B", "", ""))
dev.off()
