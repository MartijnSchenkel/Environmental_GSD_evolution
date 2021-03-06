/2020_05_11_SAY_WYY_geno_freqs_0_9999.txt

  Allele frequencies (as counts) for simulations in which Y-chromosomal fitness effects were varied.
  Generated by ../analysis/2020_05_11_SAY_WYY_data_collection.R after performing associated 
  simulations and used to generate Figures 3 and 5 using 
  ../analysis/2022_02_10_SAY_WYY_analysis_graphics.R. Variables are as follows:
 
  sex		= Sex for which allele counts are reported
  copy		= Copy for which allele counts are reported (1 = maternal, 2 = paternal)
  deme		= Deme for which allele counts are reported (0 to N-1)
  muActivAM	= mu_D, rate at which AM evolves de novo
  SAY		= s_M = -s_F, SA effect of Y-chromosome
  WYY		= S_YY, survival of YY homozygotes
  F_SE		= Number of sensitive (S) expressed (E) F alleles
  F_SU		= Number of sensitive (S) unexpressed (U) F alleles
  F_IE		= Number of insensitive (I) expressed (U) F alleles - corresponds here to frequency
		  of dominant female-determining F^I allele as discussed in manuscript.
  F_IU		= Number of insensitive (I) unexpressed (U) F alleles
  Y_E		= Number of expressed YM alleles (AKA: Y-chromosome)
  Y_U		= Number of unexpressed YM alleles (AKA: X-chromosome)
  A_E		= Number of expressed AM alleles
  A_U		= Number of unexpressed AM alleles

/2020_05_11_SAY_WYY_predictions.txt

  Predicted values from GAMs generated in ../analysis/2022_02_10_SAY_WYY_analysis_graphics.R. Used to
  generate Figure 3 and 5 (see manuscript for analysis details). 

  muActivAM	= mu_D, rate at which AM evolves de nobo
  SAY		= s_M = -s_F, SA effect of Y-chromosome
  WYY		= S_YY, survival of YY homozygotes
  FI0		= GAM-predicted frequency of FI in first deme
  FI1		= GAM-predicted frequency of FI in last deme
  YM0		= GAM-predicted frequency of YM in first deme
  YM1		= GAM-predicted frequency of YM in last deme
  AM0		= GAM-predicted frequency of AM in first deme
  AM1		= GAM-predicted frequency of AM in last deme

/2020_05_25_threshold_predictions.txt

  Predicted values from GAMs generated in ../analysis/2022_02_10_threshold_analysis_graphics.R. Used to
  generate Figure 2 (see manuscript for analysis details). 
  
  muActivAM	= mu_D, rate at which AM evolves de nobo
  BetaT		= beta_T, extent to which F expression is affected by temperature
  ThresF	= theta_F, feminization threshold
  FI0		= GAM-predicted frequency of FI in first deme
  FI1		= GAM-predicted frequency of FI in last deme

/2020_05_25_threshold_geno_freqs_0_9999.txt

  Genotype frequencies for simulations in which the feminization thresholds and temperature effects 
  were varied. Generated by ../analysis/2020_05_25_threshold_data_collection.R after performing
  associated simulations and used to generate Figure 2 using 
  ../analysis/2022_02_10_threshold_analysis_graphics.R. Variables are identical to those described 
  for 2020_05_11_SAY_WYY_geno_freqs_0_9999.txt but without SAY and WYY, which are replaced by the 
  following parameters:

  BetaT		= beta_T, as used in the F overexpression rate 1 + beta.
  ThresF	= theta_F, the lower boundary of F expression required for feminization.


/2020_06_12_13_23_16_geno.txt
  
  Genotypes generated in a single simulation run, used to generate Figure 6 along with 
  2020_06_15_Kozielska_et_al_data.xls. The parameter settings for said simulation are included in
  2020_06_12_13_23_16_parameter_values.txt. This datafile is generated using the source code found 
  in ../source_code/figure_5/. Note that the individual ID and the maternal/paternal allele are NOT
  provided in this dataset but that these can be generated anew. Each individual has 10 rows of data
  corresponding to 2 entries for F sensitivity, F expression, YM expression, AM expression, and sex
  chromosome type (5 traits x 2 = 10 entries). The first entry for each trait is for the maternal
  allele, and the second entry is the paternal allele. This applies to ALL raw genotype datasets 
  generated with the model.

  generation	= Generation at which output was generated
  deme		= Deme in which individual is found
  sex		= Sex of individual whose genes are reported
  locus		= Locus whose trait is reported
  variable	= Which trait is reported
  value		= Trait value


/2020_06_12_13_23_16_parameter_values.txt

  Parameter output file associated with 2020_06_12_13_23_16_geno.txt. This file is generated using 
  the source code found in ../source_code/figure_5/. Descriptions of parameter values are all 
  provided in the source code.


/2020_06_15_Kozielska_et_al_data.xlsx

  Data mined from:
  Kozielska, M., Feldmeyer, B., Pen, I., Weissing, F.J. & Beukeboom, L.W. (2008). Are autosomal sex-
  determining factors of the housefly (Musca domestica) spreading north? Genet. Res. Cambridge 
  90:157???165.

  Specifically, this is the data on FD frequencies and M-factors in Table 3, combined with the 
  geographic info in table 2 (lat, long, temperature - location labels not included here).

  lat		= Latitude (on Northern hemisphere)
  long		= Longitude
  females	= Number of females assayed per location
  nFD		= Number of females found to carry FD
  pFD		= Proportion of females found to carry FD versus
  males		= Number of males assayed per location
  SC		= Frequency of M-factor on sex chromosomes
  A1-A5		= Frequency of M-factor on autosomes 1 through 5.
  temperature	= Average yearly temperature

/2021_06_22_thresholds_2021_07_08_07_48_58_timelapse_allele_freqs.txt
/2021_06_22_thresholds_2021_07_08_07_48_58_timelapse_geno_freqs.txt

  Timelapse data used to generate Supplementary Figure 3 using 
  ../analysis/2021_07_13_threshold_timelapse_data_collection.R. See associated analysis file for 
  further details.

