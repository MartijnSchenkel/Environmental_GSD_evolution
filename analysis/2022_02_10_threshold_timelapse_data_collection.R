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

startgen <- 10
endgen   <- 10000
genskip <- 10
gens <- seq(startgen, endgen, genskip)

# Folder/filename parameters
simdate <- "2021_06_22_thresholds"
# datafolder <- "2019_12_simulations_analyses"
# datafolder <- "/data/p275703"
# settingsfile <- paste(simdate, "_settings.pdf", sep = "")

# Allele status
isSensitive <- c("I", "S") # Insensitive / Sensitive
isExpressed <- c("U", "E") # Unexpressed / Expressed

# ============= #
# Data analysis #
# ============= #


AS <- tibble(Status = factor(),
                copy = factor(), 
                sex = factor(), 
                locus = factor(), 
                Freq = numeric(), 
                generation = numeric())

GS <- tibble(Fs = character(),
                Ms = character(),
                n = numeric(),
                Sex = character(),
                generation = numeric())
                
                
pt <- read.table("D:/Work/HOTSEX_every_x_gens/HOTSEX_every_x_gens/2021_07_08_07_48_58_parameter_values.txt", T)
muActivAM <- pt$value[46]
thetaF <- pt$value[29]
thetaM <- pt$value[30]

sexes <- c("Females", "Males")
for(gent in gens)
{
  print(gent); print(Sys.time())

  # Load data
  dt <- read.table(str_c("D:/Work/HOTSEX_every_x_gens/HOTSEX_every_x_gens/2021_07_08_07_48_58", gent, "_geno.txt"), T)
  
  # Add parameter data.
  dt <- dt %>% mutate(muActivAM = muActivAM, thetaF = thetaF, thetaM = thetaM)
  
  
  # Define critical expression/sensitivity values for determining which alleles have which traits
  TEF <- thetaM/2 # Threshold below which an F allele is considered non-expressed
  TSF <- 0.0 # Threshold below which an F allele is considered insensitive to M
  TEM <- (2-thetaF)/2 # Threshold below which an M allele is considered non-expressed
  
  # Add unique IDs and copy info to identify different alleles in different individuals
  dt$ID <- rep(seq(1, length(dt$sex)/10), rep(10, length(dt$sex)/10)) # Unique ID for per individual
  dt$copy <- rep(rep(seq(1:2)), length(dt$sex)/2) # Maternal (1) and paternal (2) alleles
  dt <- dt %>% mutate(deme = as.numeric(deme),
                      sex = factor(sex),
                      copy = factor(copy)
                      )
  
  # Relevel factors to be easier to interpret.
  levels(dt$sex) <- c("Females", "Males", "Intersex")
  levels(dt$copy) <- c("Maternal", "Paternal")
  
  # Categorize all F alleles based on sensitivity and expression
  Fs <- dt %>% filter(locus == "F") %>% spread(key = variable, value = value) %>% 
    mutate(Insensitive = isSensitive[(Sensitivity > TSF)+1],
           Expressed = isExpressed[(Expression >= TEF)+1],
           Status = paste(Insensitive, Expressed, sep = ""))
  
  # Relabel F alleles for easier genotype categorization
  Fs2 <- Fs %>% select(-Expression, -Sensitivity, -Insensitive, -Expressed, -locus) %>% 
    spread(key = copy, value = Status) %>%
    rename(MatF = Maternal, PatF = Paternal) %>% 
    mutate(MatF2 = ifelse(MatF == "SE", "F", 
                       ifelse(MatF == "IE", "FI", "FU")),
           PatF2 = ifelse(PatF == "SE", "F", 
                          ifelse(PatF == "IE", "FI", "FU")))
  
  # Categorize all YM alleles based on expression
  Ys <- dt %>% filter(locus == "YM") %>% spread(key = variable, value = value) %>% 
    mutate(Expressed = isExpressed[(Expression >= TEM)+1],
           Status = Expressed) 
  
  # Relabel for easier genotyping
  Ys2 <- Ys %>% select(-Expression, -Expressed, -locus) %>% 
    spread(key = copy, value = Status) %>% 
    rename(MatY = Maternal, PatY = Paternal)
  
  # Categorize all AM alleles based on expression
  As <- dt %>% filter(locus == "AM") %>% spread(key = variable, value = value) %>% 
    mutate(Expressed = isExpressed[(Expression >= TEM)+1],
           Status = Expressed)
  
  # Relabel for easier genotyping
  As2 <- As %>% select(-Expression, -Expressed, -locus) %>% 
    spread(key = copy, value = Status) %>% 
    rename(MatA = Maternal, PatA = Paternal)
  
  # Combine and reprocess into genotype info
  at <- full_join(Fs2, Ys2, by = c("generation", "deme", "sex", "muActivAM", "thetaF", "thetaM", "ID"))
  at <- full_join(at, As2, by = c("generation", "deme", "sex", "muActivAM", "thetaF", "thetaM", "ID"))
  at <- at %>% mutate(MatF = factor(MatF), PatF = factor(PatF),
                      MatY = factor(MatY), PatY = factor(PatY),
                      MatA = factor(MatA), PatA = factor(PatA)) %>% 
    
    # NOTE: F alleles  are assigned numeric values so that later on we can identify the entire genotype based on unique
    # sums of these scores (e.g. F = 1 -> FF = 2).
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
    
  # Collect allele frequencies
  Fst <- data.frame(with(Fs, table(Status, copy, sex, locus)))
  Yst <- data.frame(with(Ys, table(Status, copy, sex, locus)))
  Ast <- data.frame(with(As, table(Status, copy, sex, locus)))
  allele_freqs <- Fst %>% full_join(Yst,  by = c("Status", "copy", "sex", "locus", "Freq")) %>% 
    full_join(Ast,  by = c("Status", "copy", "sex", "locus", "Freq")) %>% 
    mutate(generation = gent)

  # collect genotype frequencies
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
  
  intersex <- at %>% filter(sex == "Intersex") %>% 
    with(table(by = list(genoF, SumM))) %>% 
    as_tibble() %>%
    rename(Fs = by.1, Ms = by.2) %>%
    mutate(Sex = "Intersex")
  
  genotype_freqs <- males %>% full_join(females, by = c("Fs", "Ms", "n", "Sex")) %>%
    full_join(intersex, by = c("Fs", "Ms", "n", "Sex")) %>%
    mutate(generation = gent)
  
  AS <- full_join(AS, allele_freqs, by = c("Status", "copy", "sex", "locus", "Freq", "generation"))
  GS <- full_join(GS, genotype_freqs, by = c("Fs", "Ms", "n", "Sex", "generation"))
# } else {
#   print(gent)
#   missing <- c(missing, gent)
  # }
}

write.table(x = AS, file = paste0(simdate, "_2021_07_08_07_48_58_timelapse_allele_freqs.txt"), col.names = T, row.names = F)
write.table(x = GS, file = paste0(simdate, "_2021_07_08_07_48_58_timelapse_geno_freqs.txt"), col.names = T, row.names = F)

AS <- read.table(paste0(simdate, "_2021_07_08_07_48_58_timelapse_allele_freqs.txt"), T) %>% as_tibble()
GS <- read.table(paste0(simdate, "_2021_07_08_07_48_58_timelapse_geno_freqs.txt"), T) %>% as_tibble()

# --------------------------- #
# Allele frequencies graphics #
# --------------------------- #

AS <- AS %>% group_by(locus, copy, sex, generation) %>% mutate(RelFreq = Freq / sum(Freq))
AS2 <- AS %>% group_by(generation, sex) %>% summarize(tot = sum(Freq)) %>% 
  group_by(generation) %>% mutate(PropF = tot[1] / sum(tot),
                                  PropI = tot[2] / sum(tot),
                                  PropM = tot[3] / sum(tot),
                                  OSRF = tot[1] / (tot[1] + tot[3]),
                                  OSRM = tot[3] / (tot[1] + tot[3]) 
                                  )
AS2 %>% filter(generation <= 4000) %>% ggplot(aes(generation)) + 
  geom_line(aes(y = PropM), color = "dodgerblue") +
  # geom_line(aes(y = IntersexFreq), color = "grey5") +
  geom_line(aes(y = PropF), color = "hotpink")



pf <- AS %>% filter(locus == "F", generation <= 4000, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = str_sub(Status, 1,1), linetype = str_sub(Status, 2,2))) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("hotpink", "grey5"), labels = c("Insensitive", "Sensitive")) +
  scale_linetype(labels = c("Expressed", "Unexpressed")) +
  
  labs(color = "Sensitivity", linetype = "Expression", 
       x = NULL, y = "Frequency", title = bquote(italic(F))) +  
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

py <- AS %>% filter(locus == "YM", generation <= 4000, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("dodgerblue", "grey5"), labels = c("Expressed", "Unexpressed")) + 
  labs(x = NULL,
       y = "Frequency",
       title = bquote(italic(Y^M))) + 
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

pa <- AS %>% filter(locus == "AM", generation <= 4000, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("dodgerblue", "grey5"), labels = c("Expressed", "Unexpressed")) + 
  labs(x = "Generation",
       y = "Frequency",
       title = bquote(italic(A^M))) + 
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  


p_alleles <- plot_grid(pf,py,pa, ncol = 1)

pdf(file = paste0(simdate, "_alleles_timelapse.pdf"), height = 9, width = 6)
p_alleles
dev.off()

png(file = paste0(simdate, "_alleles_timelapse.png"), height = 9, width = 6, res = 1000, units = "in")
p_alleles
dev.off()

pf2 <- AS %>% filter(locus == "F", generation <= 250, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = str_sub(Status, 1,1), linetype = str_sub(Status, 2,2))) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("hotpink", "grey5"), labels = c("Insensitive", "Sensitive")) +
  scale_linetype(labels = c("Expressed", "Unexpressed")) +
  
  labs(color = "Sensitivity", linetype = "Expression", 
       x = NULL, y = NULL, title = bquote(italic(F))) +  
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

py2 <- AS %>% filter(locus == "YM", generation <= 250, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("dodgerblue", "grey5"), labels = c("Expressed", "Unexpressed")) + 
  labs(x = NULL,
       y = NULL,
       title = bquote(italic(Y^M))) + 
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

pa2 <- AS %>% filter(locus == "AM", generation <= 250, sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("dodgerblue", "grey5"), labels = c("Expressed", "Unexpressed")) + 
  labs(x = NULL,
       y = NULL, 
       title = bquote(italic(A^M))) + 
  facet_grid(copy ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  


p_alleles2 <- plot_grid(pf2,py2,pa2, ncol = 1)

pdf(file = paste0(simdate, "_alleles_timelapse_detail.pdf"), height = 9, width = 6)
p_alleles2
dev.off()

png(file = paste0(simdate, "_alleles_timelapse_detail.png"), height = 9, width = 6, res = 1000, units = "in")
p_alleles2
dev.off()



pf_legend <- get_legend(pf)
py_legend <- get_legend(py)
pa_legend <- get_legend(pa)

p_alleles3a <- plot_grid(pf + theme(legend.position = "none"),
                         py + theme(legend.position = "none"),
                         pa + theme(legend.position = "none"), ncol = 1)

p_alleles3b <- plot_grid(pf2 + theme(legend.position = "none"),
                         py2 + theme(legend.position = "none"),
                         pa2 + theme(legend.position = "none"), ncol = 1)
p_alleles3legend <- plot_grid(pf_legend, py_legend, pa_legend, ncol = 1)

p_alleles3 <- plot_grid(p_alleles3a, p_alleles3b, p_alleles3legend, 
                        ncol = 3, rel_widths = c(1,1,.5), 
                        labels = c("A", "B", ""))


pdf(file = paste0(simdate, "_alleles_timelapse_combined.pdf"), height = 8, width = 9)
p_alleles3
dev.off()

png(file = paste0(simdate, "_alleles_timelapse_combined.png"), height = 8, width = 9, res = 1000, units = "in")
p_alleles3
dev.off()



# ----------------------------- #
# Genotype frequencies graphics #
# ----------------------------- #

GS <- GS %>% group_by(Sex, generation) %>% 
  mutate(RelFreq = n / sum(n),
         Ms2 = ifelse(Ms == 0, "None", ifelse(Ms == 1, "Single", "Multiple")),
         Fs = factor(Fs)) 

levels(GS$Fs) <- c("F/F", "F/0", "FI/F", "FI/FI", "FI/0", "0/0")
pg <- GS %>% group_by(Fs, Sex, generation, Ms2) %>% summarize(RelFreq = sum(RelFreq)) %>% 
  filter(generation <= 4000, Sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Fs, linetype = factor(Ms2, levels = c("None", "Single", "Multiple")))) +
  geom_line(size = 1.25) + 
  facet_wrap(~Sex) + 
  labs(linetype = "Expressed M alleles",
       color = "F genotype",
       x = "Generation",
       y = "Frequency") + 
  scale_color_viridis(option = "magma", begin = 0.05, end = 0.9, discrete = T,
                      labels = c(bquote(italic(F/F)),
                                 bquote(italic(F/0)),
                                 bquote(italic(F^I/F)),
                                 bquote(italic(F^I/F^I)),
                                 bquote(italic(F^I/0)),
                                 bquote(italic(0/0)))) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

pg  

# GS %>% filter(generation <= 4000, Sex != "Intersex") %>% 
#   ggplot(aes(generation, RelFreq, color = Fs, linetype = factor(Ms))) +
#   geom_line(size = 1.25) + 
#   facet_wrap(~Sex) + 
#   labs(linetype = "Ms") + 
#   scale_color_viridis(option = "magma", begin = 0.05, end = 0.95, discrete = T)


pdf(file = paste0(simdate, "_geno_timelapse.pdf"), height = 4, width = 8)
pg
dev.off()

png(file = paste0(simdate, "_geno_timelapse.png"), height = 4, width = 8, res = 1000, units = "in")
pg
dev.off()


pg2 <- GS %>% group_by(Fs, Sex, generation, Ms2) %>% summarize(RelFreq = sum(RelFreq)) %>% 
  filter(generation <= 250, Sex != "Intersex") %>% 
  ggplot(aes(generation, RelFreq, color = Fs, linetype = factor(Ms2, levels = c("None", "Single", "Multiple")))) +
  geom_line(size = 1.25) + 
  facet_wrap(~Sex) + 
  labs(linetype = "Expressed M alleles",
       color = "F genotype",
       x = "Generation",
       y = "Frequency") + 
  scale_color_viridis(option = "magma", begin = 0.05, end = 0.9, discrete = T,
                      labels = c(bquote(italic(F/F)),
                                 bquote(italic(F/0)),
                                 bquote(italic(F^I/F)),
                                 bquote(italic(F^I/F^I)),
                                 bquote(italic(F^I/0)),
                                 bquote(italic(0/0)))) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines")) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  

pg2

# GS %>% filter(generation <= 4000, Sex != "Intersex") %>% 
#   ggplot(aes(generation, RelFreq, color = Fs, linetype = factor(Ms))) +
#   geom_line(size = 1.25) + 
#   facet_wrap(~Sex) + 
#   labs(linetype = "Ms") + 
#   scale_color_viridis(option = "magma", begin = 0.05, end = 0.95, discrete = T)


pdf(file = paste0(simdate, "_geno_timelapse_detail.pdf"), height = 4, width = 8)
pg2
dev.off()

png(file = paste0(simdate, "_geno_timelapse_detail.png"), height = 4, width = 8, res = 1000, units = "in")
pg2
dev.off()

pg_legend <- get_legend(pg)


pg3 <- plot_grid(pg + theme(legend.position = "none"),
                 pg2 + theme(legend.position = "none"), ncol = 1,
                 labels = c("A", "B"))

pg4 <- plot_grid(pg3, pg_legend, ncol = 2, rel_widths = c(1, 0.25))

pdf(file = paste0(simdate, "_geno_timelapse_combined.pdf"), height = 6, width = 8)
pg4
dev.off()

png(file = paste0(simdate, "_geno_timelapse_combined.png"), height = 6, width = 8, res = 1000, units = "in")
pg4
dev.off()
pg4

# ------------------- #
# Simplified graphics #
# ------------------- #

pf3 <- AS %>% group_by(Status, sex, locus, generation) %>% 
  summarise(RelFreq = sum(RelFreq) / 2) %>%
  mutate(locus = factor(locus, c("F", "YM","AM"))) %>%
  filter(locus == "F", sex != "Intersex", generation <= 4000) %>% 
  ggplot(aes(generation, RelFreq, color = str_sub(Status, 2,2), linetype = str_sub(Status, 1,1))) +
  geom_line(size = 1) + 
  scale_color_viridis(option = "inferno", begin = 0.1, end = 0.9, discrete = T, labels = c("Expressed", "Unexpressed"), direction = -1) + 
  scale_linetype(labels = c("Insensitive", "Sensitive")) +
  
  labs(color = "Expression", linetype = bquote(Sensitivity~(italic(F)~only)), 
       x = "", y = "Frequency", title = bquote(italic(F))) +  
  facet_wrap( ~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "lines"),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))
pf3
py3 <- AS %>% group_by(Status, sex, locus, generation) %>% 
  summarise(RelFreq = sum(RelFreq) / 2) %>%
  mutate(locus = factor(locus, c("F", "YM","AM"))) %>%
  filter(locus == "YM", sex != "Intersex", generation <= 4000) %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_viridis(option = "inferno", begin = 0.1, end = 0.9, discrete = T, c("Expressed", "Unexpressed"), direction = -1) + 
  labs(x = "",
       y = "Frequency",
       title = bquote(italic(M)[Y])) + 
  facet_wrap(~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "lines"),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01))  
py3

pa3 <- AS %>% group_by(Status, sex, locus, generation) %>% 
  summarise(RelFreq = sum(RelFreq) / 2) %>%
  mutate(locus = factor(locus, c("F", "YM","AM"))) %>%
  filter(locus == "AM", sex != "Intersex", generation <= 4000) %>% 
  ggplot(aes(generation, RelFreq, color = Status)) +
  geom_line(size = 1) + 
  scale_color_viridis(option = "inferno", begin = 0.1, end = 0.9, discrete = T, c("Expressed", "Unexpressed"), direction = -1) + 
  labs(x = "Generation",
       y = "Frequency",
       title = bquote(italic(M)[A])) + 
  facet_wrap(~ sex) + 
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing = unit(1, "lines"), 
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "lines"),
        panel.grid.minor = element_blank()) + 
  scale_y_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) 
pa3


pf_legend <- get_legend(pf3 + theme(legend.position = "bottom",
                                    legend.direction = "vertical"))



pf3d <- pf3 + scale_x_continuous(limits = c(0,250), expand = c(0.01, 0.01)) + labs(y = "")
py3d <- py3 + scale_x_continuous(limits = c(0,250), expand = c(0.01, 0.01)) + labs(y = "")
pa3d <- pa3 + scale_x_continuous(limits = c(0,250), expand = c(0.01, 0.01)) + labs(y = "")


p_alleles3a <- plot_grid(pf3, py3, pa3, ncol = 1)

p_alleles3b <- plot_grid(pf3d, py3d, pa3d, ncol = 1)
p_alleles3legend <- plot_grid(NULL,pf_legend, NULL, ncol = 1)

p_alleles3 <- plot_grid(NULL,p_alleles3a, NULL, p_alleles3b, NULL,
                        ncol = 5, labels = c( "","A", "", "B",  ""), rel_widths = c(0.025, 1, 0.05,1, 0.025))

p_alleles4 <- plot_grid(p_alleles3, p_alleles3legend,
              ncol = 1, rel_heights = c(1,0.2))
p_alleles4
pdf(file = paste0(simdate, "_alleles_timelapse_combined_simplified.pdf"), height = 6, width = 7)
p_alleles4
dev.off()

png(file = paste0(simdate, "_alleles_timelapse_combined_simplified.png"), height = 6, width = 7, res = 1000, units = "in")
p_alleles4
dev.off()

