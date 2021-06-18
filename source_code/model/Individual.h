#pragma once
#include "F_factor.h"
#include "M_factor.h"
#include <array>
#include <fstream>
#include <algorithm>
#include <optional>

// there should be explanatory comments

enum class Sex { Female, Male, Intersex };

class Individual
{
public:
    Individual();
    Individual(const Individual& Mom, const Individual& Dad);
    void setSex(const double localTemp);
    void mutateSD_old(const double Temperature);
	void mutateSD(const double Temperature);
    size_t disperse(const size_t demeNumber) const;
    void printStatistics(std::ofstream &data, const size_t Deme, const int Time) const;
    Sex sex = Sex::Female;
	bool homozygousYY();
	double getFitness();


private:
    std::array<bool, 2> YChromosomes; // ** what's the use of this? true if YM has active M; better return from function? 
	// MS: Track if M-factor is on Y-chromosome (ancestral) or X-chromosome (transposed). Ancestral Y-linked M is linked to SA allele. Also sometimes want to exclude transposition to X.
   
	std::array<F_factor, 2> F;
    std::array<std::optional<M_factor>, 2> YM; 
    std::array<std::optional<M_factor>, 2> AM;
    double Fitness = 1.0;
};
