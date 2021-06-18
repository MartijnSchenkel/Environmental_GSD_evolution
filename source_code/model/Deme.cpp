#include "Deme.h"
#include <algorithm>
#include <iostream>

// Constructor:
Deme::Deme(const double setDemeTemp) : demeTemperature(setDemeTemp)
{
    // Reserve space
	Females.reserve(PopSize);
	Males.reserve(PopSize);
	Babies.reserve(PopSize);
	BabiesDispersed.reserve(PopSize);

	// Fill up until population is saturated (max popsize is reached)
	for (int i = 0; i < PopSize;)
	{
        Individual Baby; // default constructor
		Baby.setSex(demeTemperature);
			
        if (Baby.sex == Sex::Female)
		{
            Females.push_back(std::move(Baby));
			++i;
        }
        else if (Baby.sex == Sex::Male)
		{
            Males.push_back(std::move(Baby));
			++i;
		}
		else
		{
            Intersex.push_back(std::move(Baby));
			++i;
		}
    }
}


void Deme::demeReproduction()
{
	Babies.clear();
	
	//std::cout << Males.size() << " " << Females.size() << std::endl;
    if ((Males.size() > 0) & (Females.size() > 0)) {
		
		// Calculate fitness scores
		std::vector<double> WM;
		std::vector<double> WF;
		
        for(auto& male : Males) WM.push_back(male.getFitness());
        for(auto& female : Females) WF.push_back(female.getFitness());

		std::discrete_distribution<> MaleW(WM.begin(), WM.end());
		std::discrete_distribution<> FemaleW(WF.begin(), WF.end());

        while (Babies.size() < PopSize)
		{						
            Individual Dad = Males[MaleW(rng)];
            Individual Mom = Females[FemaleW(rng)];

			Individual Baby(Mom, Dad);
			
            // Check if baby has two Y-chromosomes
            // ** why not have function Individual::has2Y()? **

            bool survived{true};
            if (Baby.homozygousYY())
                if (ru()>survivalYY) survived = false;
            if (survived){
                 Baby.mutateSD(demeTemperature);
                 Baby.setSex(demeTemperature);
                 Babies.push_back(std::move(Baby));
            }
        } // while #babies < Popsize
    } // enough moms & dads
}


void Deme::demePrintGenos(std::ofstream &data, size_t Deme, int Time) const {
    for (auto& male : Males) male.printStatistics(data, Deme, Time);
    for (auto& female : Females) female.printStatistics(data, Deme, Time);
	for (auto& intersexual : Intersex) intersexual.printStatistics(data, Deme, Time);
}

