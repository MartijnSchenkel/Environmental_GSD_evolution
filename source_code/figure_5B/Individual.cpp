#include "Individual.h"
#include "parameters.h"
#include "utils.h"
#include <iostream>

// Constructors (** without faster initializer lists **)
Individual::Individual()
{
	const bool YM1 = r2(pMY[0]);
	const bool YM2 = r2(pMY[1]);
	const bool AM1 = r2(pMA[0]);
	const bool AM2 = r2(pMA[1]);

	
	if (YM1) { YM[0] = M_factor(baseExprM); }
	if (YM2) { YM[1] = M_factor(baseExprM); }
	if (AM1) { AM[0] = M_factor(baseExprM); }
	if (AM2) { AM[1] = M_factor(baseExprM); }
	
	YChromosomes = { YM1, YM2 };
} // F default

Individual::Individual(const Individual& Mom, const Individual& Dad)
{ 
	const bool MomSC = r2();
	const bool DadSC = r2();

	F = { F_factor{ Mom.F[r2()] }, F_factor{ Dad.F[r2()] } };
    YM = { Mom.YM[MomSC], Dad.YM[DadSC] }; // ** using default copy constructor of M_factor **
    AM = { Mom.AM[r2()], Dad.AM[r2()] };   // ** same here **
    YChromosomes = { Mom.YChromosomes[MomSC], Dad.YChromosomes[DadSC] }; // ** seems a little weird at first **
}

// Get (not set???) sex of individual
void Individual::setSex(const double localTemp)
{
    // Calculate total M-factor expression
	double YM_expression = YM[0]->getExpression(localTemp) + YM[1]->getExpression(localTemp);
	double AM_expression = AM[0]->getExpression(localTemp) + AM[1]->getExpression(localTemp);

    if (location_M_expression) // ** what does this mean??? currently false **
	{
        double a1 = AM[0]->getExpression(localTemp) * ME_YA_factor; // ** partially repeating calculation... **
        double a2 = AM[1]->getExpression(localTemp) * ME_YA_factor; // ** currently * 1.0 **
		a1 = clip(a1, 0, 1);
		a2 = clip(a2, 0, 1);
		AM_expression = a1 + a2;
	}

	double produceM = YM_expression + AM_expression;
	// Calculate F-factor expression levels
	// I.e. determine product produced by each allele, and then degrade product based on its interaction with M
	// This also accounts for differences in sensitivity between the two alleles

	double deviate = rnorm(0, FSD); // calculate random deviate for F expression due to local differences in temperature.
	double Fprime0 = std::max((F[0].getExpression(localTemp) + deviate) * ( 1 - F[0].getSensitivity(localTemp) * produceM), 0.0); // net F activity of first allele
	double Fprime1 = std::max((F[1].getExpression(localTemp) + deviate) * ( 1 - F[1].getSensitivity(localTemp) * produceM), 0.0); // dito second allele
    
	double Fprime = Fprime0 + Fprime1; // total net F activity
	
	 // Set thresholds depending on whether thresholds are temperature-dependent.
	double demeThresF;
	double demeThresM;
	if (temperature_F_threshold)
	{
		 demeThresF = thresF / ((localTemp - TempMin) * threshold_Slope + 1);
		 demeThresM = thresM / ((localTemp - TempMin) * threshold_Slope + 1);
	}
	else {
		demeThresF = thresF;
		demeThresM = thresM;
	}

	// Compare net F activity to thresholds to determine the sex, and calculate fitness (mating success) based on number of Y-chromosomes.
	if (Fprime > demeThresF)
	{
        sex = Sex::Female;
        size_t nSC = YChromosomes[0] + YChromosomes[1];
		Fitness = wF[nSC];		
	}
	else {
		if (Fprime < demeThresM)
		{
            sex = Sex::Male;
            size_t nSC = YChromosomes[0] + YChromosomes[1];
			Fitness = wM[nSC];
		}
		else
		{
            sex = Sex::Intersex;
		}
	}
}

void Individual::mutateSD(const double Temperature)
{
	// Mutate F expression
	if (ru() < muExprF) { F[0].mutateExpr(); }
	if (ru() < muExprF) { F[1].mutateExpr(); }

	// Mutate F sensitivity
	if (ru() < muSensF) { F[0].mutateSens(); }
	if (ru() < muSensF) { F[1].mutateSens(); }

	// Mutate YM expression
	if (YM[0]) if (ru() < muExprM) { YM[0]->mutateExpr(); }
	if (YM[1]) if (ru() < muExprM) { YM[1]->mutateExpr(); }
	
	// Mutate de novo AM
	if (deNovoAM)
	{
		if (!AM[0]) if (ru() < muActivAM) AM[0] = rnorm(baseExprM, baseExprMSD);
		if (!AM[1]) if (ru() < muActivAM) AM[1] = rnorm(baseExprM, baseExprMSD);
	} 

	// Mutate AM expression
	if (AM[0]) if (ru() < muExprM) { AM[0]->mutateExpr(); }
	if (AM[1]) if (ru() < muExprM) { AM[1]->mutateExpr(); }

	// Transposition of M-factors (simplified - only YM to autosome transposition)
	if (transposingYM)
	{
		double TR = transpositionRate;

		// Calculate M transposition rate if it is affected by temperature
		if (temperature_M_transposition)
		{
			// Calculate expression multiplier as a result of the local temperature.
			const double TR_mult = (Temperature - TempMin) * MT_slope + 1.0;

			// Multiply base expression by expression multiplier.
			TR = TR * TR_mult;
		}

		// Logical array to keep track of which YM transposes and which does not
		std::array<bool, 2> willTranspose = { false };

		// Check if YM is expressed and on Y-chromosome (latter should be true by default)
		if (YM[0] && YChromosomes[0]) if (ru() < TR) willTranspose[0] = true;
		if (YM[1] && YChromosomes[1]) if (ru() < TR) willTranspose[1] = true;

		// Determine order in which YMs will transpose (only relevant if both alleles do so!)
		std::vector<size_t> order = { 0, 1 };
		std::shuffle(order.begin(), order.end(), rng);

		// Determine to which linkage group and which allele M will transpose
		std::array<bool, 2> targetAllele = { false };
		std::array<std::optional<M_factor>, 2> newAM = AM;

		
		for (size_t i = 0; i < 2; ++i)
		{
			if (willTranspose[order[i]])
			{
				targetAllele[order[i]] = r2(); // 0 = maternal allele, 1 = paternal allele
				newAM[targetAllele[order[i]]] = YM[order[i]];
			}
		}
	}	
	return;
};

size_t Individual::disperse(const size_t demeNumber) const
{
	if (demeNumber == 0)
	{
        if (r2(dispersalRate / 2)) return 1;
	}
    else if (demeNumber == (nDemes - 1))
	{
        if (r2(dispersalRate / 2)) return nDemes-2;
	}
    else
	{
        if (r2(dispersalRate)) { if (r2()) return demeNumber-1; else return demeNumber+1; }
	}
    return demeNumber;
}

bool Individual::homozygousYY() {
	if (YChromosomes[0] & YChromosomes[1])
	{
		return(true);
	}
	else {
		return(false);
	};
}

double Individual::getFitness()
{
	return(Fitness);
}

void Individual::printStatistics(std::ofstream &data, const size_t Deme, const int Time) const
{
    const int sx = static_cast<std::underlying_type<Sex>::type>(sex);
    data << Time << "\t" << Deme << "\t" << sx << "\tF\tSensitivity\t" << F[0].getSensitivity(0.0) << std::endl;
    data << Time << "\t" << Deme << "\t" << sx << "\tF\tSensitivity\t" << F[1].getSensitivity(0.0) << std::endl;

    data << Time << "\t" << Deme << "\t" << sx << "\tF\tExpression\t" << F[0].getExpression(0.0) << std::endl;
    data << Time << "\t" << Deme << "\t" << sx << "\tF\tExpression\t" << F[1].getExpression(0.0) << std::endl;
	
    data << Time << "\t" << Deme << "\t" << sx << "\tYM\tExpression\t" << YM[0]->getExpression(0.0) << std::endl;
    data << Time << "\t" << Deme << "\t" << sx << "\tYM\tExpression\t" << YM[1]->getExpression(0.0) << std::endl;
	
    data << Time << "\t" << Deme << "\t" << sx << "\tAM\tExpression\t" << AM[0]->getExpression(0.0) << std::endl;
    data << Time << "\t" << Deme << "\t" << sx << "\tAM\tExpression\t" << AM[1]->getExpression(0.0) << std::endl;

    data << Time << "\t" << Deme << "\t" << sx << "\tSexChromosome\tType\t" << YChromosomes[0] << std::endl;
    data << Time << "\t" << Deme << "\t" << sx << "\tSexChromosome\tType\t" << YChromosomes[1] << std::endl;
}
