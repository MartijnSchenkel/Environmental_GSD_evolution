/figure_5B
  
  Contains source code of the model as used to generate the data in Figure 5B which is based on a 
  single simulation. For details see documentation on /model below and comments in source code. 
  

/generate_parameter_files

  Contains source code for a simple C++ program that generates unique parameter files where certain
  variables are altered. Two versions are included corresponding to the key analyses presented in 
  the manuscript. To run, compile with either main*.cpp file and the randomnumbers.cpp and 
  randomnumbers.h files. The two main*.cpp files correspond to the following analyses/variables 
  affected:

  * main_threshold_beta.cpp - used for simulations on [Figure 3] 
    Variables affected (code/paper): 
    * thresF/theta_F
    * FE_diff/(1+beta_T)
    * muActivAM/[mu_A]
  * main_SAY_WYY.cpp - used for simulations on [Figure 4, Figure 6]
    Variables affected (code/paper):
    * sSAM/s_M (NB: sSAF = -sSAM / s_F = -s_M)
    * survivalYY/S_YY
    * muActivAM/[mu_D]


/job_files

  Contains the job files used for running simulations. NB: Job files all assume a certain folder 
  structure, where unique parameter.h files are located in separate folders and all source code 
  except for the parameters.h file are included in another folder. Unique executables are compiled 
  by combining these source code files in the folder containing the parameter.h folder. Upon 
  completing the simulation, a data file and a separate file listing all parameter values in a R-
  friendly format have been generated, and the entire folder is moved to a data storage location 
  because of memory constraints on cluster partitions. 


/model

  Contains source code in C++ for the model described in the manuscript. Can be compiled and 
  executed locally using e.g. Microsoft Visual Studio to run a single instance of the model. 
  Simulations on the HPC Peregrine were ran using this source code and a unique parameters.h file 
  generated using the code provided in the /generate_par_files folder. All files include comments to
  ensure proper documentation.