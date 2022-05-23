#pragma once
#include "parameters.h"
#include "randomnumbers.h"

class F_factor
{
public:
	// constructor
	F_factor() = default;
    F_factor(const double Expr, const double Sens) : Expression{Expr}, Sensitivity{Sens} {}
	// print values
    double getExpression(const double Temperature) const;
    double getSensitivity(const double Temperature) const;
	
	// mutate values
	void mutateExpr();
	void mutateSens();

    // set expression level and sensitivity to M product
    void setExpression(double Expr) { Expression = Expr; }
    void setSensitivity(double Sens) { Sensitivity = Sens; }

private:
	double Expression = baseExprF;
	double Sensitivity = baseSensF;
	bool canMutateExpr = true;
	bool canMutateSens = true;
};
