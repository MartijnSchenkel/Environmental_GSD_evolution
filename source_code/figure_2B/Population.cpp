#include "Population.h"
#include "parameters.h"
#include <iostream>

Population::Population()
{
	Demes.reserve(nDemes);
    for (int deme = 0; deme < nDemes; ++deme) Demes.emplace_back(Deme{0.0});
    // ** so intially all demes have temperature zero
}

void Population::increaseTemperatures(const size_t warmup_time)
{
	std::array<double, nDemes> deltaT; // rate of change in temperature per generation for each deme
	std::array<double, nDemes> currentT{ 0.0 }; // current temperature of each deme
	
	// Calculate rates of change in temperature
    for (size_t d = 0; d < nDemes; ++d){
        deltaT[d] = ((TempMax - TempMin)*d / (nDemes - 1))/warmup_time;
    }


	// Increase temperature.
	for (size_t t = 0; t < warmup_time; ++t)
	{
		if (t % genskip == 0)
		{
			std::cout << "Warmup: " << t << " / " << warmup_time << std::endl;
		}
		for (size_t d = 0; d < nDemes; ++d)
		{			
			Demes[d].setTemperature(currentT[d] + deltaT[d]); // Increase deme temperature
			currentT[d] = currentT[d] + deltaT[d]; // Update current temperature.
		}
		popReproduce();
	}
}

void Population::increaseTemperatures(const size_t warmup_time, std::string date, std::string time, int& tot_t)
{
	std::array<double, nDemes> deltaT; // rate of change in temperature per generation for each deme
	std::array<double, nDemes> currentT{ 0.0 }; // current temperature of each deme

	// Calculate rates of change in temperature
	for (size_t d = 0; d < nDemes; ++d) {
		deltaT[d] = ((TempMax - TempMin) * d / (nDemes - 1)) / warmup_time;
	}


	// Increase temperature.
	for (size_t t = 0; t < warmup_time; ++t)
	{
		if (t % genskip == 0)
		{
			std::cout << "Warmup: " << t << " / " << warmup_time << std::endl;
			tot_t = tot_t + genskip;

			std::ofstream geno_t(date + time + std::to_string(tot_t) + "_geno.txt");
			geno_t << "generation" << "\t" << "deme" << "\t" << "sex" << "\t" << "locus" << "\t" << "variable" << "\t" << "value" << std::endl;
			popPrintGenos(geno_t, runTime);
		}
		for (size_t d = 0; d < nDemes; ++d)
		{
			Demes[d].setTemperature(currentT[d] + deltaT[d]); // Increase deme temperature
			currentT[d] = currentT[d] + deltaT[d]; // Update current temperature.
		}
		popReproduce();
	}
}


void Population::popReproduce()
{
	// Within-deme reproduction, followed by dispersal.
    for (size_t i = 0; i < nDemes; ++i)
	{ 	
		// Reproduction within each deme.		
		Demes[i].demeReproduction(); 

		// Dispersal 
        for (size_t j = 0, sz =Demes[i].Babies.size(); j < sz; ++j)
		{	
            size_t toDeme = Demes[i].Babies[j].disperse(i);
            Demes[toDeme].BabiesDispersed.push_back(std::move(Demes[i].Babies[j]));
		}	

        Demes[i].Babies.clear();
    }



	// Survival to adults and sorting into males and females.
    for (auto& d : Demes)
	{
		// Kill all adults from previous generation.
        d.Females.clear();
        d.Males.clear();
        d.Intersex.clear();

        // might as well shuffle them once rather than picking random indices over and over again
        std::shuffle(d.BabiesDispersed.begin(),d.BabiesDispersed.end(),rng);
        size_t sz{d.BabiesDispersed.size()};

        for (size_t num_moved{0}, fplusm{0}; num_moved<sz && fplusm<PopSize; ++num_moved)
		{

            if (d.BabiesDispersed[num_moved].sex == Sex::Female){
                d.Females.push_back(std::move(d.BabiesDispersed[num_moved]));
                ++fplusm;
            }
            else if (d.BabiesDispersed[num_moved].sex == Sex::Male){
                d.Males.push_back(std::move(d.BabiesDispersed[num_moved]));
                ++fplusm;
            }
            else d.Intersex.push_back(std::move(d.BabiesDispersed[num_moved]));
        }

        d.BabiesDispersed.clear();
	}
}

void Population::popPrintGenos(std::ofstream &data, int Time)
{
    for (size_t d = 0; d < nDemes; ++d)
	{
		Demes[d].demePrintGenos(data, d, Time);
	}
}
