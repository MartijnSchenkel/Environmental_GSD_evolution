#pragma once
#include <array>
#include "randomnumbers.h"
// perhaps cleaner to store constants in classes that can become static class members downstream

// Controls for which procedures affect SD in what ways
constexpr bool temperature_F_expression = false; // F expression temperature-sensitive
constexpr bool temperature_F_threshold = !temperature_F_expression; // F expression thresholds temperature-sensitive

// Controls for how AM expression evolves (de novo or through YM transposition)
constexpr bool deNovoAM = true;
constexpr bool transposingYM = !deNovoAM;

// Runtime parameters
constexpr int burnin = 1000; // Number of generations without temperature differences in model
constexpr int warmup_time = 1000; // Number of generations during which temperature increases
constexpr int runTime = 20000; // Number of generations to run the model after temperature is increased
constexpr int genskip = 100; // Number of generations between writing output

// Population parameters
constexpr int PopSize = 1000; // Population size per deme
constexpr int nDemes = 11; // Number of demes
constexpr double dispersalRate = 0.1; // Rate at which dispersal occurs (halved for first and last deme, see disperse() function in Individual.cpp)

// Temperature parameters
constexpr double TempMax = 1; // Maximum temperature
constexpr double TempMin = 0; // Minimum temperature

// Fitness parameters for Y-chromosomal fitness effects
constexpr double survivalYY = 0.9; // proportion of YY individuals that survive

constexpr double dSAF = 0.5; // dominance of SA effect in females
constexpr double dSAM = 0.5; // dominance of SA effect in males
constexpr double sSAM = ru_const(); // fitness effect of SA in males
constexpr double sSAF = -sSAM; // fitness effect of SA in females

// Fitness arrays for SA locus on XY, dependent on number of Y chromosomes (w*[0] = XX, w*[1] = XY, w*[2] = YY)
constexpr std::array<double, 3> wF{ { 1.0, 1.0 + dSAF * sSAF, 1.0 + sSAF } };
constexpr std::array<double, 3> wM{ { 1.0, 1.0 + dSAM * sSAM, 1.0 + sSAM } };

// Initial allele frequencies of M-factors
constexpr std::array<double, 2> pMY{ { 0.0, 0.5 } }; // Initial frequency of M on Y on maternal and paternal copies
constexpr std::array<double, 2> pMA{ { 0.0, 0.0 } }; // Initial frequency of M on A on maternal and paternal copies

// Sex determination parameters
constexpr double thresF = 1.2; // Threshold for F expression, higher => female. Below or above 1 (ALWAYS well below 2).
constexpr double thresM = 0.3; // Threshold for F expression, lower => male; if thresM < expression < thresF then => intersex.
constexpr double FSD = 0.05; // Random error term for F expression, introduces a bit of noise in model.

// F-factor parameters
constexpr double baseExprF = 0.7; // Initial expression level of F
constexpr double baseSensF = 0.95; // Initial sensitivity of F

// M-factor parameters
constexpr double baseExprM = 0.9; // Initial expression level of M
constexpr double baseExprMSD = 0.05; // SD for initial M expression, only used for de novo evolution of AM

// Mutation rate parameters
constexpr double muNull = 0.01; // Rate at which null mutations arise.

// F expression mutation
constexpr double muExprF = 0.1; // Mutation rate for F expression level
constexpr double sigmaExprF = 0.1; // SD for mutation effect for F expression level
constexpr double muExprFNull = muNull; // Proportion of F expression mutations that is a null mutation

// F sensitivity mutation
constexpr double muSensF = 0.1; // Mutation rate for F sensitivity
constexpr double sigmaSensF = 0.1; // SD for mutation effect for F sensitivity
constexpr double muSensFNull = muNull;  // Proportion of F sensitivity mutations that is a null mutation

// M expression mutation
constexpr double muExprM = 0.1; // Mutation rate for M expression level
constexpr double sigmaExprM = 0.1; // Mutation rate for M expression level
constexpr double muExprMNull = 0.0;  // Proportion of M expression mutations that is a null mutation

constexpr double muActivAM = 0.002; // Rate at which AM becomes expressed (de novo M-factor evolution)
constexpr double transpositionRate = 0.01; // Rate at which M-factors will attempt to transpose to different genomic location

// F expression ~ temperature parameters
constexpr double FE_diff = 1.5; // Fold difference in F expression between minimal and maximal temperature.
constexpr double betaT = (FE_diff - 1) / (TempMax - TempMin); // Slope at which F expression changes when temperature increases by 1 unit.

// F threshold ~ temperature parameters
constexpr double thresholdScale = 1.5; // Rate by which threshold is reduced between minimal and maximal temperature
constexpr double threshold_Slope = (thresholdScale - 1) / (TempMax - TempMin); // slope for reduction in threshold per 1 unit temperature
																   
// NB: Following parameters are associated with other versions of the model, primarily different mechanisms of temperature effects
// However, these models did not produce any noteworthy results, and are only included here for completeness' sake

// Controls for which effects on SD are included in the model (temperature effects and/or M location effect)
constexpr bool temperature_F_sensitivity = false; // F sensitivity temperature-sensitive
constexpr bool temperature_M_expression = false; // M expression temperature-sensitive
constexpr bool temperature_M_transposition = false; // M transposition temperature-sensitive
constexpr bool location_M_expression = false; // M expressed multiplied on autosome

// F sensitivity ~ temperature parameters
constexpr double FS_diff = 1.0; // Fold difference in F sensitivity
constexpr double FS_slope = (FS_diff - 1) / (TempMax - TempMin); // Slope for F sensitivity

// M expression ~ temperature parameters
constexpr double ME_diff = 1.0; // Fold difference in M expression.
constexpr double ME_slope = (ME_diff - 1) / (TempMax - TempMin);

// M transposition ~ temperature parameters
constexpr double MT_diff = 1.0; // Fold difference in M transposition rate.
constexpr double MT_slope = (MT_diff - 1) / (TempMax - TempMin);

// M expression ~ location parameters
constexpr double ME_YA_factor = 1.0; // Factor by which expression of an M-factor on A is amplified relative to its expression if it were on Y

// Logical control whether M can transpose to the X-chromosome (used in mutateSD_old() function)
constexpr bool M_on_X = 0; // logical control if M can be on X or not 