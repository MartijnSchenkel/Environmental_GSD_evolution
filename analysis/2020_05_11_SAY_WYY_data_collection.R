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
library(wesanderson)

# library(mgcv)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

startsim <- 0
endsim   <- 9999

# Folder/filename parameters
simdate <- "2020_05_11_SAY_WYY"

# datafolder <- ""
datafolder <- "/data/p275703"
settingsfile <- paste(simdate, "_settings.pdf", sep = "")

# Parameter value ranges
deNovoRate <- c(0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04, 0.05)
SAY <- c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
FitnessYY <- seq(1, 0, -0.1)

length(SAY)
N <- 1

# Scenarios
CYY <- c("", "CYY_")
NXM <- c("", "NXM_")
LME <- c("", "LME_")
effect_types <- c("TFE")#, "TFS", "TME", "TMT")

# Simulation and analysis parameters
r <- c(0, 0.1) # rates for all mutation rates (muExprF, muSensF, muExprM, deNovoRate)
effect <- c(F, T) # whether or not a specific effect (cost of YY, no M on X, and/or location of M) is included in the simulation

# Allele threshold values

# Allele status
isSensitive <- c("I", "S") # Insensitive / Sensitive
isExpressed <- c("U", "E") # Unexpressed / Expressed
isActive <- c("I", "A") # Inactive / Active


# ============= #
# Data analysis #
# ============= #

freq_data <- as_tibble()
DS <- as_tibble()

beta_FD <- matrix(data = NA, ncol = 6)
beta_YM <- matrix(data = NA, ncol = 6)
beta_AM <- matrix(data = NA, ncol = 8)
sexes <- c("Females", "Males")

TEF <- 0.15 # Threshold below which an F allele is considered non-expressed
TSF <- 0.0 # Threshold below which an F allele is considered insensitive to M
TEM <- (2-1.2)/2 # Threshold below which an M allele is considered non-expressed


for(n in startsim:endsim)
{
    print(n); print(Sys.time())
    ddir <- paste(datafolder, simdate, "/", n,"/", sep = "")
    if(length(dir(ddir, full.names = T, pattern="geno\\.txt$")) > 0)
    {
      dt <- cbind(read.table(dir(ddir, full.names = T, pattern="geno\\.txt$")[1], header = T))
      pt <- cbind(read.table(dir(ddir, full.names = T, pattern="values\\.txt$")[1], header = T)) 
      dt <- dt %>% mutate(muActivAM = pt$value[46], SAY = pt$value[17], WYY = pt$value[14])
      
      dt$ID <- rep(seq(1, length(dt$sex)/10), rep(10, length(dt$sex)/10)) # Unique ID for per individual
      dt$copy <- rep(rep(seq(1:2)), length(dt$sex)/2) # Maternal (1) and paternal (2) alleles
      dt <- dt %>% mutate(deme = as.numeric(deme),
                          sex = factor(sex)
                          )
      levels(dt$sex) <- c("Females", "Males", "Intersex")
      
      Fs <- dt %>% filter(locus == "F") %>% spread(key = variable, value = value) %>% 
        mutate(Insensitive = isSensitive[(Sensitivity > TSF)+1],
               Expressed = isExpressed[(Expression >= TEF)+1],
               Status = paste("F_", Insensitive, Expressed, sep = "")) %>%
        with(table(sex, copy, deme, muActivAM, SAY, WYY, factor(Status, levels = c("F_SE", "F_SU", "F_IE", "F_IU")))) %>%
        as.data.frame() %>% spread(key = Var7, value = Freq)
      
      Ys <- dt %>% filter(locus == "YM") %>% spread(key = variable, value = value) %>% 
        mutate(Expressed = isExpressed[(Expression >= TEM)+1],
               Status = paste("Y_", Expressed, sep = "")) %>%
        with(table(sex, copy, deme, muActivAM, SAY, WYY, factor(Status, levels = c("Y_E", "Y_U")))) %>%
        as.data.frame() %>% spread(key = Var7, value = Freq)
      
      As <- dt %>% filter(locus == "AM") %>% spread(key = variable, value = value) %>% 
        mutate(Expressed = isExpressed[(Expression >= TEM)+1],
               Status = paste("A_", Expressed, sep = "")) %>%
        with(table(sex, copy, deme, muActivAM, SAY, WYY, factor(Status, levels = c("A_E", "A_U")))) %>%
        as.data.frame() %>% spread(key = Var7, value = Freq)
      
      at <- full_join(Fs, Ys)
      at <- full_join(at, As)

      
      DS <- rbind(DS, at)
  }
 else {
  print(ddir)
  }
}

write.table(x = DS, file = paste(simdate, "_geno_freqs_", startsim, "_", endsim, ".txt", sep = ""), row.names = F, col.names = T)