#pragma once
#include "Individual.h"
#include <vector>
#include <fstream>
#include "utils.h"

class Deme
{
public:
    Deme(const double setDemeTemp);
	void demeReproduction(); // generate babies from females and males.
    void demePrintGenos(std::ofstream &data, size_t Deme, int Time) const;
    void setTemperature(const double temp) { demeTemperature = temp; }

	double demeTemperature;
	   
	std::vector<Individual> Females;
	std::vector<Individual> Males;
	std::vector<Individual> Babies;
	std::vector<Individual> BabiesDispersed;
    std::vector<Individual> Intersex;
	std::vector<double> MaleFitness;
	std::vector<double> FemaleFitness;
};
