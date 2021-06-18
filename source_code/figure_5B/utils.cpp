#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <vector>
#include <cmath>
#include "utils.h"
#include <numeric>
#include <cassert>
#define _CRT_SECURE_NO_WARNINGS

void printDate(std::string &sDate)
{
	time_t now = time(0);
	tm *local_time = localtime(&now);

	int iY = local_time->tm_year + 1900;
	int iM = local_time->tm_mon + 1;
	int iD = local_time->tm_mday;

	std::string sYear = std::to_string(iY);
	std::string sMonth = std::to_string(iM);
	std::string sDay = std::to_string(iD);


	std::ostringstream YearOut;
	YearOut << std::setw(4) << std::setfill('0') << sYear;
	sYear = YearOut.str();

	std::ostringstream MonthOut;
	MonthOut << std::setw(2) << std::setfill('0') << sMonth;
	sMonth = MonthOut.str();

	std::ostringstream DayOut;
	DayOut << std::setw(2) << std::setfill('0') << sDay;
	sDay = DayOut.str();

	sDate = (sYear + "_" + sMonth + "_" + sDay + "_");
	return;
}

void printTime(std::string &sTime)
{
	time_t now = time(0);
	tm *local_time = localtime(&now);

	int iH = local_time->tm_hour;
	int iM = local_time->tm_min;
	int iS = local_time->tm_sec;


	std::string sHour = std::to_string(iH);
	std::string sMinute = std::to_string(iM);
	std::string sSecond = std::to_string(iS);


	std::ostringstream HourOut;
	HourOut << std::setw(2) << std::setfill('0') << sHour;
	sHour = HourOut.str();

	std::ostringstream MinuteOut;
	MinuteOut << std::setw(2) << std::setfill('0') << sMinute;
	sMinute = MinuteOut.str();

	std::ostringstream SecondOut;
	SecondOut << std::setw(2) << std::setfill('0') << sSecond;
	sSecond = SecondOut.str();

	sTime = (sHour + ":" + sMinute + ":" + sSecond);
	return;
}

void printTimeFileName(std::string &sTime)
{
	time_t now = time(0);
	tm *local_time = localtime(&now);

	int iH = local_time->tm_hour;
	int iM = local_time->tm_min;
	int iS = local_time->tm_sec;

	std::string sHour = std::to_string(iH);
	std::string sMinute = std::to_string(iM);
	std::string sSecond = std::to_string(iS);

	std::ostringstream HourOut;
	HourOut << std::setw(2) << std::setfill('0') << sHour;
	sHour = HourOut.str();

	std::ostringstream MinuteOut;
	MinuteOut << std::setw(2) << std::setfill('0') << sMinute;
	sMinute = MinuteOut.str();

	std::ostringstream SecondOut;
	SecondOut << std::setw(2) << std::setfill('0') << sSecond;
	sSecond = SecondOut.str();

	sTime = (sHour + "_" + sMinute + "_" + sSecond);
	return;
};


double vAvg(std::vector<double>& v)
{
    size_t sz{v.size()};
    assert(sz>0);

    return std::accumulate(v.begin(),v.end(),0.0)/sz;
};

double vSD(std::vector<double>& v, double avg)
{
    size_t sz{v.size()};
    assert(sz>0);

    double variance = std::accumulate(v.begin(),v.end(),0.0,[&avg,&sz](auto acc, const auto& x){ return acc+(x-avg)*(x-avg)/sz; });

    return(sqrt(variance));
};

double clip(double x, const double bottom, const double top)
{
    if(x < bottom) return(bottom);
    if (x > top) return(top);
    return(x);
}

double clip_positive(double x)
{
    if (x < 0.0) return(0.0);
    return(x);
}

// Print parameter values to file.
void printParameters(std::ofstream& pf, const long seed)
{
	// To Ido: Is there a way to just write an entire header file as output? I didn't find one, hence this very ugly function...
	
	// I guess I would first store all constants in a struct. Then using some magic you can iterate over the struct members, maybe...
        // https://github.com/apolukhin/magic_get

	pf << "temperature_F_expression" << "\t" << temperature_F_expression << std::endl;
	pf << "temperature_F_threshold" << "\t" << temperature_F_threshold << std::endl;
	pf << "deNovoAM" << "\t" << deNovoAM << std::endl;
	pf << "transposingYM" << "\t" << transposingYM << std::endl;
	pf << "burnin" << "\t" << burnin << std::endl;
	pf << "warmup_time" << "\t" << warmup_time << std::endl;
	pf << "runTime" << "\t" << runTime << std::endl;
	pf << "genskip" << "\t" << genskip << std::endl;
	pf << "PopSize" << "\t" << PopSize << std::endl;
	pf << "nDemes" << "\t" << nDemes << std::endl;
	pf << "dispersalRate" << "\t" << dispersalRate << std::endl;
	pf << "TempMax" << "\t" << TempMax << std::endl;
	pf << "TempMin" << "\t" << TempMin << std::endl;
	pf << "survivalYY" << "\t" << survivalYY << std::endl;
	pf << "dSAF" << "\t" << dSAF << std::endl;
	pf << "dSAM" << "\t" << dSAM << std::endl;
	pf << "sSAM" << "\t" << sSAM << std::endl;
	pf << "sSAF" << "\t" << sSAF << std::endl;
	pf << "wF_XX" << "\t" << wF[0] << std::endl;
	pf << "wF_XY" << "\t" << wF[1] << std::endl;
	pf << "wF_YY" << "\t" << wF[2] << std::endl;
	pf << "wM_XX" << "\t" << wM[0] << std::endl;
	pf << "wM_XY" << "\t" << wM[1] << std::endl;
	pf << "wM_YY" << "\t" << wM[2] << std::endl;
	pf << "pMY" << "\t" << pMY[0] << std::endl;
	pf << "pMY" << "\t" << pMY[1] << std::endl;
	pf << "pMA" << "\t" << pMA[0] << std::endl;
	pf << "pMA" << "\t" << pMA[1] << std::endl;
	pf << "thresF" << "\t" << thresF << std::endl;
	pf << "thresM" << "\t" << thresM << std::endl;
	pf << "FSD" << "\t" << FSD << std::endl;
	pf << "baseExprF" << "\t" << baseExprF << std::endl;
	pf << "baseSensF" << "\t" << baseSensF << std::endl;
	pf << "baseExprM" << "\t" << baseExprM << std::endl;
	pf << "baseExprMSD" << "\t" << baseExprMSD << std::endl;
	pf << "muNull" << "\t" << muNull << std::endl;
	pf << "muExprF" << "\t" << muExprF << std::endl;
	pf << "sigmaExprF" << "\t" << sigmaExprF << std::endl;
	pf << "muExprFNull" << "\t" << muExprFNull << std::endl;
	pf << "muSensF" << "\t" << muSensF << std::endl;
	pf << "sigmaSensF" << "\t" << sigmaSensF << std::endl;
	pf << "muSensFNull" << "\t" << muSensFNull << std::endl;
	pf << "muExprM" << "\t" << muExprM << std::endl;
	pf << "sigmaExprM" << "\t" << sigmaExprM << std::endl;
	pf << "muExprMNull" << "\t" << muExprMNull << std::endl;
	pf << "muActivAM" << "\t" << muActivAM << std::endl;
	pf << "transpositionRate" << "\t" << transpositionRate << std::endl;
	pf << "FE_diff" << "\t" << FE_diff << std::endl;
	pf << "betaT" << "\t" << betaT << std::endl;
	pf << "thresholdScale" << "\t" << thresholdScale << std::endl;
	pf << "threshold_Slope" << "\t" << threshold_Slope << std::endl;
	pf << "temperature_F_sensitivity" << "\t" << temperature_F_sensitivity << std::endl;
	pf << "temperature_M_expression" << "\t" << temperature_M_expression << std::endl;
	pf << "temperature_M_transposition" << "\t" << temperature_M_transposition << std::endl;
	pf << "location_M_expression" << "\t" << location_M_expression << std::endl;
	pf << "FS_diff" << "\t" << FS_diff << std::endl;
	pf << "FS_slope" << "\t" << FS_slope << std::endl;
	pf << "ME_diff" << "\t" << ME_diff << std::endl;
	pf << "ME_slope" << "\t" << ME_slope << std::endl;
	pf << "MT_diff" << "\t" << MT_diff << std::endl;
	pf << "MT_slope" << "\t" << MT_slope << std::endl;
	pf << "ME_YA_factor" << "\t" << ME_YA_factor << std::endl;
	pf << "M_on_X" << "\t" << M_on_X << std::endl;
}
