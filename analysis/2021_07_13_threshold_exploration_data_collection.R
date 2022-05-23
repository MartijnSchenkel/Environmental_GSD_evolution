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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

startsim <- 0
endsim   <- 9999

# Folder/filename parameters
simdate <- "2021_07_13_thresholds"
# datafolder <- "2019_12_simulations_analyses"
# datafolder <- "/data/p275703"
settingsfile <- paste(simdate, "_settings.pdf", sep = "")

# Allele status
isSensitive <- c("I", "S") # Insensitive / Sensitive
isExpressed <- c("U", "E") # Unexpressed / Expressed


# ============= #
# Data analysis #
# ============= #

freq_data <- as_tibble()
DS <- as_tibble(c(generation = numeric(),
                  deme = numeric(), 
                  sex = factor(), 
                  muActivAM = numeric(), 
                  thetaF = numeric(), 
                  thetaM = numeric(), 
                  ID = numeric(),
                  MatF = factor(),
                  PatF = factor(),
                  MatY = factor(),
                  PatY = factor(),
                  MatA = factor(),
                  PatA = factor()))

sexes <- c("Females", "Males")
for(nt in startsim:endsim)
{
    print(nt); print(Sys.time())
    ddir <- paste("E:/2020_03_temperature_model/data/p275703/", simdate, "/", nt,"/", sep = "")
    if(length(dir(ddir, full.names = T, pattern="geno\\.txt$")) > 0)
    {
      dt <- cbind(read.table(dir(ddir, full.names = T, pattern="geno\\.txt$")[1], header = T))
      pt <- cbind(read.table(dir(ddir, full.names = T, pattern="values\\.txt$")[1], header = T)) 
      
      
      if(length(dt[,1] > 0))
      {
        muActivAM <- pt$value[46]
        thetaF <- pt$value[29]
        thetaM <- pt$value[30]
        
        dt <- dt %>% mutate(muActivAM = muActivAM, thetaF = thetaF, thetaM = thetaM)
        
        TEF <- thetaM/2 # Threshold below which an F allele is considered non-expressed
        TSF <- 0.0 # Threshold below which an F allele is considered insensitive to M
        TEM <- (2-thetaF)/2 # Threshold below which an M allele is considered non-expressed
        
        dt$ID <- rep(seq(1, length(dt$sex)/10), rep(10, length(dt$sex)/10)) # Unique ID for per individual
        dt$copy <- rep(rep(seq(1:2)), length(dt$sex)/2) # Maternal (1) and paternal (2) alleles
        dt <- dt %>% mutate(deme = as.numeric(deme),
                            sex = factor(sex),
                            copy = factor(copy)
                            )
        levels(dt$sex) <- c("Females", "Males", "Intersex")
        levels(dt$copy) <- c("Maternal", "Paternal")
        Fs <- dt %>% filter(locus == "F") %>% spread(key = variable, value = value) %>% 
          mutate(Insensitive = isSensitive[(Sensitivity > TSF)+1],
                 Expressed = isExpressed[(Expression >= TEF)+1],
                 Status = paste(Insensitive, Expressed, sep = ""))
        
        Fs2 <- Fs %>% select(-Expression, -Sensitivity, -Insensitive, -Expressed, -locus) %>% 
          spread(key = copy, value = Status) %>%
          rename(MatF = Maternal, PatF = Paternal) %>% 
          mutate(MatF2 = ifelse(MatF == "SE", "F", 
                             ifelse(MatF == "IE", "FI", "FU")),
                 PatF2 = ifelse(PatF == "SE", "F", 
                                ifelse(PatF == "IE", "FI", "FU")))
  
        Ys <- dt %>% filter(locus == "YM") %>% spread(key = variable, value = value) %>% 
          mutate(Expressed = isExpressed[(Expression >= TEM)+1],
                 Status = Expressed) 
        
        Ys2 <- Ys %>% select(-Expression, -Expressed, -locus) %>% 
          spread(key = copy, value = Status) %>% 
          rename(MatY = Maternal, PatY = Paternal)
        
        
        As <- dt %>% filter(locus == "AM") %>% spread(key = variable, value = value) %>% 
          mutate(Expressed = isExpressed[(Expression >= TEM)+1],
                 Status = Expressed)
        
        As2 <- As %>% select(-Expression, -Expressed, -locus) %>% 
          spread(key = copy, value = Status) %>% 
          rename(MatA = Maternal, PatA = Paternal)
        
        
        at <- full_join(Fs2, Ys2, by = c("generation", "deme", "sex", "muActivAM", "thetaF", "thetaM", "ID"))
        at <- full_join(at, As2, by = c("generation", "deme", "sex", "muActivAM", "thetaF", "thetaM", "ID")) %>% 
              mutate(MatF = factor(MatF), PatF = factor(PatF),
                 MatY = factor(MatY), PatY = factor(PatY),
                 MatA = factor(MatA), PatA = factor(PatA),
                 Replicate = nt) %>% 
          mutate(MatY = ifelse(MatY == "U",0,1),
                 PatY = ifelse(PatY == "U",0,1),
                 MatA = ifelse(MatA == "U",0,1),
                 PatA = ifelse(PatA == "U",0,1),
                 MatF = ifelse(MatF2 == "F", 1, ifelse(MatF2 == "FI", 3, 0)),
                 PatF = ifelse(PatF2 == "F", 1, ifelse(PatF2 == "FI", 3, 0))) %>%
          mutate(SumF = MatF + PatF,
                 SumM = MatY + PatY + MatA + PatA) %>% 
          mutate(genoF = ifelse(SumF == 0, "UU",
                                ifelse(SumF == 1, "FU",
                                       ifelse(SumF == 2, "FF",
                                              ifelse(SumF == 3, "IU",
                                                     ifelse(SumF == 4, "IF", "II"))))))
          
        males <- at %>% filter(sex == "Males") %>% 
          with(table(by = list(genoF, SumM))) %>% 
          as_tibble() %>%
          rename(Fs = by.1, Ms = by.2) %>%
          mutate(Sex = "Males")
  
        females <- at %>% filter(sex == "Females") %>% 
          with(table(by = list(genoF, SumM))) %>% 
          as_tibble() %>%
          rename(Fs = by.1, Ms = by.2) %>%
          mutate(Sex = "Females")
        
        gt <- males %>% full_join(females, by = c("Fs", "Ms", "n", "Sex")) %>%
          mutate(muActivAM = muActivAM, thetaF = thetaF, thetaM = thetaM, replicate = nt)
        
        DS <- rbind(DS, gt)
      } else {
        print(ddir)
        missing <- c(missing, nt)
      } 
    } else {
      print(ddir)
      missing <- c(missing, nt)
    }
}

rm(missing)
for(nt in startsim:endsim)
{
  ddir <- paste("D:/2020_03_temperature_model/data/p275703/", simdate, "/", nt,"/", sep = "")
  if(length(dir(ddir, full.names = T, pattern="geno\\.txt$")) > 0)
  {
    dt <- cbind(read.table(dir(ddir, full.names = T, pattern="geno\\.txt$")[1], header = T))
    if(length(dt$generation) == 0) 
    { 
      missing <- c(missing, nt)
      print(nt)
    }
    if(length(dir(ddir, full.names = T, pattern="values\\.txt$")) == 0) 
    { 
      missing <- c(missing, nt)
      print(nt)
    }
  } else {
    missing <- c(missing, nt)
    print(nt)
  }
}

cbind(missing)
length(missing)

# Save to output
write.table(x = DS, file = paste0(simdate, "_geno_freqs_", startsim, "_", endsim, ".txt"), row.names = F, col.names = T)

# Save failed sims to output for checking
fails <- cbind(missing) %>% matrix()
write.table(x = fails, file = paste0(simdate, "_geno_freqs_", startsim, "_", endsim, "_MISSING.txt"), row.names = F,col.names = T)

# 
# DS3 <- read.table(paste0(simdate, "_geno_freqs_0_7499.txt"), T)
# DS3 %>% as_tibble()
# DST <- DS %>% mutate(Ms = as.numeric(as.character(Ms)))
# 
# DS4 <- DS3 %>% full_join(DST) %>% as_tibble()
# 
# 
# write.table(x = DS4, file = paste0(simdate, "_geno_freqs_", startsim, "_", endsim, ".txt"), row.names = F,col.names = T)
# 
# 
# missing


# ============= #
# DATA ANALYSIS #
# ============= #

DS <- read.table(paste0(simdate, "_geno_freqs_", startsim, "_", endsim, ".txt"), T)

DS <- DS %>% group_by(replicate) %>% mutate(Proportion = n / sum(n))


DS2 <- DS %>% group_by(replicate, Sex) %>% arrange(-Proportion) %>% slice(1)
DS3 <- DS2 %>% group_by(replicate) %>% mutate(Fgeno = str_c(Fs[1], "/", Fs[2]),
                                       Mgeno = str_c(Ms[1], "/", Ms[2]),
                                       PropF = Proportion[1],
                                       PropM = Proportion[2]) 

DS4 <- DS3 %>% mutate(RelActiv = ifelse(thetaF < 1, "MFA", ifelse(thetaM < 1, "MAF", "AMF")))
DS4 %>% with(table(Fgeno, Mgeno, RelActiv)) / 2
mean(DS2$PropF)
mean(DS2$PropM)


t <- DS2 %>% with(table(Fgeno, Mgeno, RelActiv)) / 2
t2 <- data.frame(t) %>% tibble()
sum(t2$Freq)


length(unique(DS$replicate))
seq(0,9999)[!(seq(0, 9999) %in% unique(DS$replicate))]

