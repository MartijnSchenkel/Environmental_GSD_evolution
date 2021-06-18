#include "Population.h"
#include "randomnumbers.h"
#include "utils.h"
#include <iostream>
#include <algorithm>
#include "parameters.h"
#include <random>
#include <iterator>
#include <vector>

int main()
{	
	// SET UP SIMULATION //
	long seed = randomize();
	std::string date; printDate(date);
	std::string time; printTimeFileName(time);

	// Create parameter file and save parameter values and seed to file
	std::ofstream parameter_file(date + time + "_parameter_values.txt");
	parameter_file << "variable" << "\t" << "value" << std::endl;
	printParameters(parameter_file, seed);

	// Initialize population
	Population Pop;

	// Perform burnin and warmup period
	std::cout << "Perform burnin" << std::endl;
	for (int t = 0; t < burnin; ++t)
	{
		Pop.popReproduce();
		if(t % genskip == 0)
		{
			std::cout << "Burnin: " << t << " / " << burnin << std::endl;
		}
	}
	std::cout << "Burnin complete" << std::endl << "Initiate warming up" << std::endl;
	Pop.increaseTemperatures(warmup_time);
	std::cout << "Warmup complete" << std::endl << "Initiate simulation" << std::endl;	

	// RUN SIMULATION //
	for (int t = 1; t <= runTime; ++t)
	{
		
		Pop.popReproduce();

		if (t % genskip == 0)
		{
			std::cout << "Simulation: " << t << " / " << runTime << std::endl;
        }
    }




	// Save all genotypes in simulation.
	std::ofstream geno(date + time + "_geno.txt");
	geno << "generation" << "\t" << "deme" << "\t" << "sex" << "\t" << "locus" << "\t" << "variable" << "\t" << "value" << "\t" << "evolvable" << "\t" << "wildtype" << std::endl;
	Pop.popPrintGenos(geno, runTime);


	return(0);
}
