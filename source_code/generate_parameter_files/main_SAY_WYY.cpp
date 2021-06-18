#include <iostream>
#include <fstream>
#include "randomnumbers.h"
#include <string>
#include <filesystem>
int main() 
{
	// SET UP SIMULATION //
	long seed = randomize();

	// Setting up shotgun-sampled SAY and WYY effects...

	constexpr double maxSA = 0.1;
	constexpr double maxTP = 0.05;
	
	int n_tot = 10000;
	for (int n = 0; n < n_tot; ++n)
	{	
		std::filesystem::create_directory("2020_05_11_SAY_WYY");
		std::filesystem::create_directory("2020_05_11_SAY_WYY/" + std::to_string(n));
		std::ofstream pf("2020_05_11_SAY_WYY/" + std::to_string(n) + "/parameters.h");
		pf << "#pragma once" << std::endl;
		pf << "#include <array>" << std::endl << std::endl;

		pf << "// perhaps cleaner to store constants in classes that can become static class members downstream" << std::endl << std::endl;

		pf << "constexpr long seed_generate = " << seed << ";" << std::endl << std::endl;

		pf << "// Controls for which procedures affect SD in what ways" << std::endl;
		pf << "constexpr bool temperature_F_expression = true; // F expression temperature-sensitive" << std::endl;
		pf << "constexpr bool temperature_F_threshold = !temperature_F_expression; // F expression thresholds temperature-sensitive" << std::endl << std::endl;

		pf << "// Controls for how AM expression evolves (de novo or through YM transposition)" << std::endl;
		pf << "constexpr bool deNovoAM = true;" << std::endl;
		pf << "constexpr bool transposingYM = !deNovoAM;" << std::endl << std::endl;

		pf << "// Runtime parameters" << std::endl;
		pf << "constexpr int burnin = 10000; // Number of generations without temperature differences in model" << std::endl;
		pf << "constexpr int warmup_time = 10000; // Number of generations during which temperature increases" << std::endl;
		pf << "constexpr int runTime = 200000; // Number of generations to run the model after temperature is increased" << std::endl;
		pf << "constexpr int genskip = 100; // Number of generations between writing output" << std::endl << std::endl;

		pf << "// Population parameters" << std::endl;
		pf << "constexpr int PopSize = 10000; // Population size per deme" << std::endl;
		pf << "constexpr int nDemes = 11; // Number of demes" << std::endl;
		pf << "constexpr double dispersalRate = 0.1; // Rate at which dispersal occurs (halved for first and last deme, see disperse() function in Individual.cpp)" << std::endl << std::endl;

		pf << "// Temperature parameters" << std::endl;
		pf << "constexpr double TempMax = 1; // Maximum temperature" << std::endl;
		pf << "constexpr double TempMin = 0; // Minimum temperature" << std::endl << std::endl;

		pf << "// Fitness parameters for Y-chromosomal fitness effects" << std::endl;
		pf << "constexpr double survivalYY = "<< ru() << "; // proportion of YY individuals that survive" << std::endl << std::endl;

		pf << "constexpr double dSAF = 0.5; // dominance of SA effect in females" << std::endl;
		pf << "constexpr double dSAM = 0.5; // dominance of SA effect in males" << std::endl;
		pf << "constexpr double sSAM = "<< ru() * maxSA <<"; // fitness effect of SA in males" << std::endl;
		pf << "constexpr double sSAF = -sSAM; // fitness effect of SA in females" << std::endl << std::endl;

		pf << "// Fitness arrays for SA locus on XY, dependent on number of Y chromosomes (w*[0] = XX, w*[1] = XY, w*[2] = YY)" << std::endl;
		pf << "constexpr std::array<double, 3> wF{ { 1.0, 1.0 + dSAF * sSAF, 1.0 + sSAF } };" << std::endl;
		pf << "constexpr std::array<double, 3> wM{ { 1.0, 1.0 + dSAM * sSAM, 1.0 + sSAM } };" << std::endl << std::endl;

		pf << "// Initial allele frequencies of M-factors" << std::endl;
		pf << "constexpr std::array<double, 2> pMY{ { 0.0, 0.5 } }; // Initial frequency of M on Y on maternal and paternal copies" << std::endl;
		pf << "constexpr std::array<double, 2> pMA{ { 0.0, 0.0 } }; // Initial frequency of M on A on maternal and paternal copies" << std::endl << std::endl;

		pf << "// Sex determination parameters" << std::endl;
		pf << "constexpr double thresF = 1.2; // Threshold for F expression, higher => female. Below or above 1 (ALWAYS well below 2)." << std::endl;
		pf << "constexpr double thresM = 0.3; // Threshold for F expression, lower => male; if thresM < expression < thresF then => intersex." << std::endl;
		pf << "constexpr double FSD = 0.05; // Random error term for F expression, introduces a bit of noise in model." << std::endl << std::endl;

		pf << "// F-factor parameters" << std::endl;
		pf << "constexpr double baseExprF = 0.7; // Initial expression level of F" << std::endl;
		pf << "constexpr double baseSensF = 0.95; // Initial sensitivity of F" << std::endl << std::endl;

		pf << "// M-factor parameters" << std::endl;
		pf << "constexpr double baseExprM = 0.9; // Initial expression level of M" << std::endl;
		pf << "constexpr double baseExprMSD = 0.05; // SD for initial M expression, only used for de novo evolution of AM" << std::endl << std::endl;

		pf << "// Mutation rate parameters" << std::endl;
		pf << "constexpr double muNull = 0.01; // Rate at which null mutations arise." << std::endl << std::endl;

		pf << "// F expression mutation" << std::endl;
		pf << "constexpr double muExprF = 0.1; // Mutation rate for F expression level" << std::endl;
		pf << "constexpr double sigmaExprF = 0.1; // SD for mutation effect for F expression level" << std::endl;
		pf << "constexpr double muExprFNull = muNull; // Proportion of F expression mutations that is a null mutation" << std::endl << std::endl;

		pf << "// F sensitivity mutation" << std::endl;
		pf << "constexpr double muSensF = 0.1; // Mutation rate for F sensitivity" << std::endl;
		pf << "constexpr double sigmaSensF = 0.1; // SD for mutation effect for F sensitivity" << std::endl;
		pf << "constexpr double muSensFNull = muNull;  // Proportion of F sensitivity mutations that is a null mutation" << std::endl << std::endl;

		pf << "// M expression mutation" << std::endl;
		pf << "constexpr double muExprM = 0.1; // Mutation rate for M expression level" << std::endl;
		pf << "constexpr double sigmaExprM = 0.1; // Mutation rate for M expression level" << std::endl;
		pf << "constexpr double muExprMNull = 0.0;  // Proportion of M expression mutations that is a null mutation" << std::endl << std::endl;

		pf << "constexpr double muActivAM = "<< ru() * maxTP << "; // Rate at which AM becomes expressed (de novo M-factor evolution)" << std::endl;
		pf << "constexpr double transpositionRate = muActivAM; // Rate at which M-factors will attempt to transpose to different genomic location" << std::endl << std::endl;

		pf << "// F expression ~ temperature parameters" << std::endl;
		pf << "constexpr double FE_diff = 1.5; // Fold difference in F expression between minimal and maximal temperature." << std::endl;
		pf << "constexpr double betaT = (FE_diff - 1) / (TempMax - TempMin); // Slope at which F expression changes when temperature increases by 1 unit." << std::endl << std::endl;

		pf << "// F threshold ~ temperature parameters" << std::endl;
		pf << "constexpr double thresholdScale = 1.5; // Rate by which threshold is reduced between minimal and maximal temperature" << std::endl;
		pf << "constexpr double threshold_Slope = (thresholdScale - 1) / (TempMax - TempMin); // slope for reduction in threshold per 1 unit temperature" << std::endl << std::endl;

		pf << "// NB: Following parameters are associated with other versions of the model, primarily different mechanisms of temperature effects" << std::endl;
		pf << "// However, these models did not produce any noteworthy results, and are only included here for completeness' sake" << std::endl << std::endl;

		pf << "// Controls for which effects on SD are included in the model (temperature effects and/or M location effect)" << std::endl;
		pf << "constexpr bool temperature_F_sensitivity = false; // F sensitivity temperature-sensitive" << std::endl;
		pf << "constexpr bool temperature_M_expression = false; // M expression temperature-sensitive" << std::endl;
		pf << "constexpr bool temperature_M_transposition = false; // M transposition temperature-sensitive" << std::endl;
		pf << "constexpr bool location_M_expression = false; // M expressed multiplied on autosome" << std::endl << std::endl;

		pf << "// F sensitivity ~ temperature parameters" << std::endl;
		pf << "constexpr double FS_diff = 1.0; // Fold difference in F sensitivity" << std::endl;
		pf << "constexpr double FS_slope = (FS_diff - 1) / (TempMax - TempMin); // Slope for F sensitivity" << std::endl << std::endl;

		pf << "// M expression ~ temperature parameters" << std::endl;
		pf << "constexpr double ME_diff = 1.0; // Fold difference in M expression." << std::endl;
		pf << "constexpr double ME_slope = (ME_diff - 1) / (TempMax - TempMin);" << std::endl << std::endl;

		pf << "// M transposition ~ temperature parameters" << std::endl;
		pf << "constexpr double MT_diff = 1.0; // Fold difference in M transposition rate." << std::endl;
		pf << "constexpr double MT_slope = (MT_diff - 1) / (TempMax - TempMin);" << std::endl << std::endl;

		pf << "// M expression ~ location parameters" << std::endl;
		pf << "constexpr double ME_YA_factor = 1.0; // Factor by which expression of an M-factor on A is amplified relative to its expression if it were on Y" << std::endl << std::endl;

		pf << "// Logical control whether M can transpose to the X-chromosome (used in mutateSD_old() function)" << std::endl;
		pf << "constexpr bool M_on_X = 0; // logical control if M can be on X or not" << std::endl << std::endl;

		pf.close();
	} 
	return(0);
}