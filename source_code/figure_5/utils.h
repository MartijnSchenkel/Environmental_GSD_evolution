#include <string>
#include <vector>
#include "parameters.h"

// Function to get the date, used to write output with the starting date of the simulation to simplify data management.
void printDate(std::string &sDate);
void printTime(std::string &sTime);
void printTimeFileName(std::string &sTime);

// Calculate mean and SD for a vector of doubles.
double vAvg(std::vector<double>& v);
double vSD(std::vector<double>& v, double avg);

// Clip function
double clip(double x, const double bottom, const double top);
double clip_positive(double x);

// Print parameter values to file.
void printParameters(std::ofstream &pf, const long seed);
