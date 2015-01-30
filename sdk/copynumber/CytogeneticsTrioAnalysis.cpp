////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License (version 2) as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program;if not, write to the
//
// Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/**
 * @file CytogeneticsTrioAnalysis.cpp
 *
 * @brief This file contains some Trio Analysis functions for Cytogenetics.
 */

/**
 * @brief  Some Trio Analysis functions for Cytogenetics.
 *
 */

#include "copynumber/CytogeneticsTrioAnalysis.h"
//

#include "util/Util.h"
//

/**
 * Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vKnownParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 * @return - The calculated LOD score.
 */
double CytogeneticsTrioAnalysis::calculateLODPaternity(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownMotherGenotypeCalls, std::vector<char>& vAllegedFatherGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, double dErrorRate)
{
    if (vSampleGenotypeCalls.size() != vKnownMotherGenotypeCalls.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}
    if (vSampleGenotypeCalls.size() != vAllegedFatherGenotypeCalls.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}
    if (vSampleGenotypeCalls.size() != vProbabilityAAlleles.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}
    if (vSampleGenotypeCalls.size() != vProbabilityBAlleles.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}

    std::vector<float> vLikelihoodsH1(vSampleGenotypeCalls.size());
    std::vector<float> vLikelihoodsH2(vSampleGenotypeCalls.size());
    calculateLikelihoodH1(vSampleGenotypeCalls, vKnownMotherGenotypeCalls, vAllegedFatherGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, vLikelihoodsH1, dErrorRate);
    calculateLikelihoodH2(vSampleGenotypeCalls, vKnownMotherGenotypeCalls, vAllegedFatherGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, vLikelihoodsH2, dErrorRate);
    double dLOD = 0;
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        dLOD += log(vLikelihoodsH1[uiIndex] / vLikelihoodsH2[uiIndex]);
    }
    return dLOD;
}

/**
 * Calculate the LOD (Likelihood of odds) value for the paternity test with one parent.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 * @return - The calculated LOD score.
 */
double CytogeneticsTrioAnalysis::calculateLODPaternity(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, double dErrorRate)
{
    if (vSampleGenotypeCalls.size() != vAllegedParentGenotypeCalls.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}
    if (vSampleGenotypeCalls.size() != vProbabilityAAlleles.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}
    if (vSampleGenotypeCalls.size() != vProbabilityBAlleles.size()) {Err::errAbort("CytogeneticsTrioAnalysis::calculateLODPaternity() All vector parameters must be of the same size.");}

    std::vector<float> vLikelihoodsH1(vSampleGenotypeCalls.size());
    std::vector<float> vLikelihoodsH2(vSampleGenotypeCalls.size());
    calculateLikelihoodH1(vSampleGenotypeCalls, vAllegedParentGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, vLikelihoodsH1, dErrorRate);
    calculateLikelihoodH2(vSampleGenotypeCalls, vAllegedParentGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, vLikelihoodsH2, dErrorRate);
    double dLOD = 0;
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        dLOD += log(vLikelihoodsH1[uiIndex] / vLikelihoodsH2[uiIndex]);
    }
    return dLOD;
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vKnownParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param vLikelihoods - std::vector<float> The likelihoods calculated by this function (output).
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 */
void CytogeneticsTrioAnalysis::calculateLikelihoodH1(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownMotherGenotypeCalls, std::vector<char>& vAllegedFatherGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate)
{
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        char cSampleGenotypeCall = vSampleGenotypeCalls[uiIndex];
        char cKnownMotherGenotypeCall = vKnownMotherGenotypeCalls[uiIndex];
        char cAllegedFatherGenotypeCall = vAllegedFatherGenotypeCalls[uiIndex];
        if (cSampleGenotypeCall == -1) {cSampleGenotypeCall = 1;}
        if (cKnownMotherGenotypeCall == -1) {cKnownMotherGenotypeCall = 1;}
        if (cAllegedFatherGenotypeCall == -1) {cAllegedFatherGenotypeCall = 1;}
        vLikelihoods[uiIndex] = (float)calculateLikelihoodH1(cSampleGenotypeCall, cKnownMotherGenotypeCall, cAllegedFatherGenotypeCall, vProbabilityAAlleles[uiIndex], vProbabilityBAlleles[uiIndex], dErrorRate);
    }
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cKnownParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation.
 * @return - The calculated likelihood.
 */
double CytogeneticsTrioAnalysis::calculateLikelihoodH1(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate)
{

    return (1.0 - dErrorRate)*(1.0 - dErrorRate)*(1.0 - dErrorRate) *
        calculateTransProbability(cSampleGenotypeCall, cKnownMotherGenotypeCall, cAllegedFatherGenotypeCall) *
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        (dErrorRate * (1.0 - dErrorRate)*(1.0 - dErrorRate)) *
        (calculateTransProbability(cSampleGenotypeCall, cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateTransProbability(cSampleGenotypeCall, cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        (dErrorRate*dErrorRate * (1.0 - dErrorRate)) *
        (calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        dErrorRate*dErrorRate*dErrorRate;
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vKnownParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param vLikelihoods - std::vector<float> The likelihoods calculated by this function (output).
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 */
void CytogeneticsTrioAnalysis::calculateLikelihoodH2(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownMotherGenotypeCalls, std::vector<char>& vAllegedFatherGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate)
{
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        char cSampleGenotypeCall = vSampleGenotypeCalls[uiIndex];
        char cKnownMotherGenotypeCall = vKnownMotherGenotypeCalls[uiIndex];
        char cAllegedFatherGenotypeCall = vAllegedFatherGenotypeCalls[uiIndex];
        if (cSampleGenotypeCall == -1) {cSampleGenotypeCall = 1;}
        if (cKnownMotherGenotypeCall == -1) {cKnownMotherGenotypeCall = 1;}
        if (cAllegedFatherGenotypeCall == -1) {cAllegedFatherGenotypeCall = 1;}
        vLikelihoods[uiIndex] = (float)calculateLikelihoodH2(cSampleGenotypeCall, cKnownMotherGenotypeCall, cAllegedFatherGenotypeCall, vProbabilityAAlleles[uiIndex], vProbabilityBAlleles[uiIndex], dErrorRate);
    }
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cKnownParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation.
 * @return - The calculated likelihood.
 */
double CytogeneticsTrioAnalysis::calculateLikelihoodH2(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate)
{
    return (1.0 - dErrorRate)*(1.0 - dErrorRate)*(1.0 - dErrorRate) *
        calculateTransProbability(cSampleGenotypeCall, cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        (dErrorRate * (1.0 - dErrorRate)*(1.0 - dErrorRate)) *
        (calculateTransProbability(cSampleGenotypeCall, cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        (dErrorRate*dErrorRate * (1.0 - dErrorRate)) *
        (calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cKnownMotherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cAllegedFatherGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        dErrorRate*dErrorRate*dErrorRate;
}

/**
 * Calculate the LOD (Likelihood of odds) value for the paternity test with one parent.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param vLikelihoods - std::vector<float> The likelihoods calculated by this function (output).
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 */
void CytogeneticsTrioAnalysis::calculateLikelihoodH1(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate)
{
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        char cSampleGenotypeCall = vSampleGenotypeCalls[uiIndex];
        char cAllegedParentGenotypeCall = vAllegedParentGenotypeCalls[uiIndex];
        if (cSampleGenotypeCall == -1) {cSampleGenotypeCall = 1;}
        if (cAllegedParentGenotypeCall == -1) {cAllegedParentGenotypeCall = 1;}
        vLikelihoods[uiIndex] = (float)calculateLikelihoodH1(cSampleGenotypeCall, cAllegedParentGenotypeCall, vProbabilityAAlleles[uiIndex], vProbabilityBAlleles[uiIndex], dErrorRate);
    }
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation.
 * @return - The calculated likelihood.
 */
double CytogeneticsTrioAnalysis::calculateLikelihoodH1(char cSampleGenotypeCall, char cAllegedParentGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate)
{
    return (1.0 - dErrorRate)*(1.0 - dErrorRate) *
        calculateTransProbability(cSampleGenotypeCall, cAllegedParentGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedParentGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        (dErrorRate * (1.0 - dErrorRate)) *
        (calculateGenotypeProbability(cAllegedParentGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        dErrorRate*dErrorRate;
}

/**
 * Calculate the LOD (Likelihood of odds) value for the paternity test with one parent.
 * Note: all the vectors input to this function must be of the same length. One element per marker.
 * @param vSampleGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vAllegedParentGenotypeCalls - std::vector<char> where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param vProbabilityAAlleles - std::vector<float> The probability of the AAllele occurring in the population.
 * @param vProbabilityBAlleles - std::vector<float> The probability of the BAllele occurring in the population.
 * @param vLikelihoods - std::vector<float> The likelihoods calculated by this function (output).
 * @param dErrorRate - The error rate associated with this calculation (defaults to 0.05).
 */
void CytogeneticsTrioAnalysis::calculateLikelihoodH2(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate)
{
    for (unsigned int uiIndex = 0; (uiIndex < vSampleGenotypeCalls.size()); uiIndex++)
    {
        char cSampleGenotypeCall = vSampleGenotypeCalls[uiIndex];
        char cAllegedParentGenotypeCall = vAllegedParentGenotypeCalls[uiIndex];
        if (cSampleGenotypeCall == -1) {cSampleGenotypeCall = 1;}
        if (cAllegedParentGenotypeCall == -1) {cAllegedParentGenotypeCall = 1;}
        vLikelihoods[uiIndex] = (float)calculateLikelihoodH2(cSampleGenotypeCall, cAllegedParentGenotypeCall, vProbabilityAAlleles[uiIndex], vProbabilityBAlleles[uiIndex], dErrorRate);
    }
}

/**
 * A sub function of Calculate the LOD (Likelihood of odds) value for the paternity test with both parents.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @param dErrorRate - The error rate associated with this calculation.
 * @return - The calculated likelihood.
 */
double CytogeneticsTrioAnalysis::calculateLikelihoodH2(char cSampleGenotypeCall, char cAllegedParentGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate)
{
    return (1.0 - dErrorRate)*(1.0 - dErrorRate) *
        calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) *
        calculateGenotypeProbability(cAllegedParentGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        (dErrorRate * (1.0 - dErrorRate)) *
        (calculateGenotypeProbability(cAllegedParentGenotypeCall, fProbabilityAAllele, fProbabilityBAllele) +
        calculateGenotypeProbability(cSampleGenotypeCall, fProbabilityAAllele, fProbabilityBAllele)) +
        dErrorRate*dErrorRate;
}

/**
 * Calculate the probability of a genotype call.
 * @param cGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @return - The calculated probability.
 */
double CytogeneticsTrioAnalysis::calculateGenotypeProbability(char cGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele)
{
    switch(cGenotypeCall)
    {
    case 0: return (fProbabilityAAllele * fProbabilityAAllele); break;
    case 1: return (2.0 * fProbabilityAAllele * fProbabilityBAllele); break;
    case 2: return (fProbabilityBAllele * fProbabilityBAllele); break;
    default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
    }
    return 0;
}

/**
 * Calculate the transfer probability of a genotype call.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cKnownParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @return - The calculated probability.
 */
double CytogeneticsTrioAnalysis::calculateTransProbability(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall)
{
    switch(cSampleGenotypeCall)
    {
    case 0:
        switch(cKnownMotherGenotypeCall)
        {
        case 0:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 1; break;
            case 1: return 0.5; break;
            case 2: return 0; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 1:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0.5; break;
            case 1: return 0.25; break;
            case 2: return 0; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 2:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0; break;
            case 1: return 0; break;
            case 2: return 0; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    case 1:
        switch(cKnownMotherGenotypeCall)
        {
        case 0:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0; break;
            case 1: return 0.5; break;
            case 2: return 1; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 1:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0.5; break;
            case 1: return 0.5; break;
            case 2: return 0.5; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 2:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 1; break;
            case 1: return 0.5; break;
            case 2: return 0; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    case 2:
        switch(cKnownMotherGenotypeCall)
        {
        case 0:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0; break;
            case 1: return 0; break;
            case 2: return 0; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 1:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0; break;
            case 1: return 0.25; break;
            case 2: return 0.5; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        case 2:
            switch (cAllegedFatherGenotypeCall)
            {
            case 0: return 0; break;
            case 1: return 0.5; break;
            case 2: return 1; break;
            default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
            }
            break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
    }
    return 0;
}

/**
 * Calculate the transfer probability of a genotype call.
 * @param cSampleGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param cAllegedParentGenotypeCall - char where -1 = No Call, 0 = AA, 1 = AB, and 2 = BB.
 * @param fProbabilityAAllele - float The probability of the AAllele occurring in the population.
 * @param fProbabilityBAllele - float The probability of the BAllele occurring in the population.
 * @return - The calculated probability.
 */
double CytogeneticsTrioAnalysis::calculateTransProbability(char cSampleGenotypeCall, char cAllegedParentGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele)
{
    switch(cSampleGenotypeCall)
    {
    case 0:
        switch(cAllegedParentGenotypeCall)
        {
        case 0: return fProbabilityAAllele; break;
        case 1: return (fProbabilityAAllele / 2.0); break;
        case 2: return 0; break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    case 1:
        switch(cAllegedParentGenotypeCall)
        {
        case 0: return fProbabilityBAllele; break;
        case 1: return ((fProbabilityAAllele + fProbabilityBAllele) / 2.0); break;
        case 2: return fProbabilityAAllele; break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    case 2:
        switch(cAllegedParentGenotypeCall)
        {
        case 0: return 0; break;
        case 1: return (fProbabilityBAllele / 2.0); break;
        case 2: return fProbabilityBAllele; break;
        default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
        }
        break;
    default: Err::errAbort("Genotype Call not supported in Trio Analysis.");
    }
    return 0;
}
