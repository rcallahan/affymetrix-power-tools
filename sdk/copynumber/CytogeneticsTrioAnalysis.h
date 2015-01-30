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
 * @file CytogeneticsTrioAnalysis.h
 *
 * @brief This header contains some Trio Analysis functions for Cytogenetics.
 *        Some Trio Analysis functions for Cytogenetics.
 */

#ifndef __cytogeneticstrioanalysis_h__
#define __cytogeneticstrioanalysis_h__

#include "util/Err.h"

/**
 * @brief  A class for doing trio analysis.
 *
 */
class CytogeneticsTrioAnalysis
{
public:
    double calculateLODPaternity(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownParentGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, double dErrorRate = 0.05);
    double calculateLODPaternity(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, double dErrorRate = 0.05);

private:
    double calculateGenotypeProbability(char cObservedGenotype, float fProbabilityAlleleA, float fProbabilityAlleleB);
    double calculateTransProbability(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall);
    double calculateTransProbability(char cSampleGenotypeCall, char cAllegedParentGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele);
    double calculateLikelihoodH1(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate = 0.055);
    void calculateLikelihoodH1(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownMotherGenotypeCalls, std::vector<char>& vAllegedFatherGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate = 0.05);
    double calculateLikelihoodH2(char cSampleGenotypeCall, char cKnownMotherGenotypeCall, char cAllegedFatherGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate = 0.05);
    void calculateLikelihoodH2(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vKnownMotherGenotypeCalls, std::vector<char>& vAllegedFatherGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate = 0.05);
    double calculateLikelihoodH1(char cSampleGenotypeCall, char cAllegedParentGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate = 0.055);
    void calculateLikelihoodH1(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate = 0.05);
    double calculateLikelihoodH2(char cSampleGenotypeCall, char cAllegedParnetGenotypeCall, float fProbabilityAAllele, float fProbabilityBAllele, double dErrorRate = 0.05);
    void calculateLikelihoodH2(std::vector<char>& vSampleGenotypeCalls, std::vector<char>& vAllegedParentGenotypeCalls, std::vector<float>& vProbabilityAAlleles, std::vector<float>& vProbabilityBAlleles, std::vector<float>& vLikelihoods, double dErrorRate = 0.05);
};

#endif // __cytogeneticstrioanalysis_h__
