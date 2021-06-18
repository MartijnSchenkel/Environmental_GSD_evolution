#pragma once
#include "Deme.h"
#include <fstream>

class Population
{
public:
	Population();
    void increaseTemperatures(const size_t warmup_time);
	void popReproduce();
	void popPrintGenos(std::ofstream &data, int Time);
private:
	std::vector<Deme> Demes;
};
