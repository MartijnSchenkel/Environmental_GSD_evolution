/HOTSEX_TFE_new_2020_05_11_SAY_WYY.sh

  Used to run simulations in which the fitness effects of the Y-chromosome are included. Assumes 
  that the source code with the exception of the required "parameters.h" file are located in a 
  folder:

  $HOME/2019_02_HOTSEX/source_min_parameters

  This source code is combined n = 10000 times with a unique "parameters.h" file which can be 
  generated using ../generate_parameter_files/main_SAY_WYY.cpp + associated files, which also 
  generates the required folder structure (referred to as $SIMDIR) which looks like this:

  $HOME/2019_02_HOTSEX/data/$DATE/$n

  In which $DATE is an identifier for the simulation batch (here: 2020_05_11_SAY_WYY) and $n is the 
  unique identifier for a separate simulation in said batch (here: 0 through 9999). The source code
  combined with the unique "parameters.h" file is compiled using GCC (v. 8.2.0-2.31.1, see line 13 
  in case you want to change this) into an executable HOTSEX located in $SIMDIR. A SLURM job script 
  named TFE$DATE$n.sh is then generated and executed which runs the HOTSEX executable to perform the
  simulation. Each simulation yields a file with the genotype information of all individuals in the
  population in the final generation (ending in "geno.txt") and a parameter value file (ending in
  "parameter_values.txt") which are used for data analysis.

/HOTSEX_TFE_new_2020_05_25_threshold

  Similar to the HOTSEX_TFE_new_2020_05_11_SAY_WYY.sh script, but with $DATE = 2020_05_25_threshold
  instead. The unique parameter files and folder structure must be generated with the file 
  ../generate_parameter_files/main_threshold_beta.cpp + associated files.

