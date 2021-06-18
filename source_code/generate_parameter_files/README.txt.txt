This folder contains the files necessary to generate the parameters.h files for the simulations. It
contains two main files which correspond to the following analyses:

/main_SAY_WYY.cpp 
  Used for simulations in which the Y-chromosome has SA effects and/or costs of homozygosity. For 
  each value of n in 0:9999, a folder /n is generated, and within it a parameter.h file containing
  a unique combination of parameter values for:
  
  name in code		name in paper   Range	
  * muActivAM		mu_D		[0, 0.05]
  * sSAM (=-sSAF)	s_M = -s_F	[0, 0.1]
  * survivalYY		S_YY		[0,1]

  Each of these parameter values is sampled from a uniform distribution with the ranges specified
  above. The resulting parameters.h files are used to run the simulations for Figures 4 and 6.


/main_threshold_beta.cpp 
  Used for simulations in which the feminization threshold and temperature effect are varied. For
  each value of n in 0:9999, a folder /n is generated, and within it a parameter.h file containing
  a unique combination of parameter values for:

  name in code		name in paper   Range	
  * muActivAM		mu_D		[0, 0.05]
  * thresF		theta_F	[0.8, 1.5]
  * betaT		1 + beta_T	[1.0, 1.8]

  Each of these parameter values is sampled from a uniform distribution with the ranges specified
  above. The resulting parameters.h files are used to run the simulations for Figure 3.


These files were compiled and executed separately along with randomnumbers.cpp and randomnumbers.h 
in Microsoft Visual Studio 2017 (v15.9.21) using the C++17 standard. This generates a folder 
/2020_05_11_SAY_WYY/ (for main_SAY_WYY.cpp) or /2020_05_25_threshold/ (for main_threshold_beta.cpp),
each of which contains folders 0-9999, which in turn contain a unique parameters.h file for a single
simulation. The /2020_05_11_SAY_WYY/ and /2020_05_25_threshold/ folders were transferred to the HPC
Peregrine computer cluster and stored in a folder $HOME/2019_02_HOTSEX/data/ after which the job 
scripts found in ../job_files/ could be used to run the required simulations. Note that attempts to
run this procedure on the HPC Peregrine cluster were unsuccesful for reasons unclear, and hence were
this step was executed locally instead.
