#include "M_factor.h"
#include "utils.h"
#include <iostream>

M_factor::M_factor(const double Expr)
{
	Expression = Expr;
	isWildType = false;
	return;
}

double M_factor::getExpression(const double Temperature) const
{
	// Calculate M expression if it is affected by temperature
    if (temperature_M_expression) // currently false
	{
		// Calculate expression multiplier as a result of the local temperature.
        const double ExprMult = (Temperature - TempMin) * ME_slope + 1;

		// Multiply base expression by expression multiplier.
		double ExprRealized = Expression * ExprMult;
		//return clip(ExprRealized, 0, 1);
		return clip_positive(ExprRealized);
	}
	else
	{
		return(Expression);
	}
}

// mutate values
void M_factor::mutateExpr() 
{	

	// Check if allele can still mutate
    if (canMutateExpr & (sigmaExprM > 0))
	{
		// Check if it will develop a null mutation (loss-of-function); if yes set to 0 and avoid further mutation of trait.
		if (ru() < muExprMNull)
		{
			Expression = 0.0;
			canMutateExpr = false;
		}

		else {
			// Sample change
			double delta = rnorm(0, sigmaExprM);
			
			// Set new expression level
			if (Expression + delta <= 0.0) { canMutateExpr = false; }
            Expression = clip(Expression + delta, 0, 1);
		}
	}
}

double M_factor::getCanMutateExpr() const
{
	return(canMutateExpr);
}

bool M_factor::getWildTypeStatus() const
{
	return(isWildType);
}