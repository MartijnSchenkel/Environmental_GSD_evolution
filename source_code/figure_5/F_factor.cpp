#include "F_factor.h"
#include "utils.h"

// print values
double F_factor::getExpression(const double Temperature) const
{	
	// Calculate F expression if it is affected by temperature
    if (temperature_F_expression) // currently true
	{
		// Calculate expression multiplier as a result of the local temperature.
		const double ExprMult = 1 + betaT * Temperature;
		
		// Multiply base expression by expression multiplier.
		double ExprRealized = Expression * ExprMult;
		ExprRealized = clip_positive(ExprRealized);
		return ExprRealized;
	}

    return Expression;
}

// Print sensitivity
double F_factor::getSensitivity(const double Temperature) const
{
	// Calculate F sensitivity if it is affected by temperature
	if (temperature_F_sensitivity)
	{
		// Calculate expression multiplier as a result of the local temperature.
		const double SensMult = (Temperature - TempMin) * FS_slope + 1;

		// Multiply base expression by expression multiplier.
		double SensRealized = Sensitivity * SensMult;
		return clip(SensRealized, 0, 1);
	}
	
	return(Sensitivity);
}

// mutate values
void F_factor::mutateExpr()
{
	// Check if allele can still mutate
    if (canMutateExpr & (sigmaExprF > 0.0))
	{
		// Check if it will develop a null mutation (loss-of-function); if yes set to 0 and avoid further mutation of trait.
		if (ru() < muExprFNull)
		{
			Expression = 0.0;
			canMutateExpr = false;
		}

		else {
			// Sample change
			double delta = rnorm(0, sigmaExprF);

			// Set new expression level
			if (Expression + delta <= 0.0) { canMutateExpr = false; }
			Expression = clip(Expression + delta, 0, 1);
		}
	}
}

void F_factor::mutateSens()
{
	// Check if allele can still mutate
    if (canMutateSens & (sigmaSensF > 0.0))
	{
		// Check if it will develop a null mutation (loss-of-function); if yes set to 0 and avoid further mutation of trait.
		if (ru() < muSensFNull)
		{
			Sensitivity = 0.0;
			canMutateSens = false;
		}

		else {
			// Sample change
			double delta = rnorm(0, sigmaSensF);

			// Set new expression level
			if (Sensitivity + delta <= 0.0) { canMutateSens = false; } // 
            Sensitivity = clip(Sensitivity + delta, 0, 1);
		}
	}
}

