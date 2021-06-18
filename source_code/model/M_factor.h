#pragma once
#include "parameters.h"
#include "randomnumbers.h"

class M_factor
{
public:
	M_factor() = default;
	M_factor(const double Expr);
    double getExpression(const double Temperature) const; // print values
	void mutateExpr(); // mutate values	
    void setExpression(const double Expr) { Expression = Expr; } // set locus expression level
	double getCanMutateExpr() const;
	bool getWildTypeStatus() const;
    
private:
    double Expression = baseExprM; // initial expression level
    bool canMutateExpr = true;  // false after null-mutation
	bool isWildType = true;
};
