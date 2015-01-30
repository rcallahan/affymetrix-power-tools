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

#ifndef _CNFamilialAnalysisMethodLOD_H_
#define _CNFamilialAnalysisMethodLOD_H_
/**
 * @file CNFamilialAnalysisMethodLOD.h
 *
 * @brief This header contains the CNFamilialAnalysisMethodLOD class definition.
 */

#include "copynumber/CNFamilialAnalysisMethod.h"

/**
 * @brief Class to store allele frequency data
 */
class AlleleFrequency
{
public:
    AffxString SNPID;
    double AAlleleFrequency;
    double BAlleleFrequency;
};

/**
 * @brief Comparison predicate for AlleleFrequency pointers
 */
struct AlleleFreqComp : std::binary_function<const AlleleFrequency*, const AlleleFrequency*, bool>
{
    bool operator()(const AlleleFrequency* lhs, const AlleleFrequency* rhs) const {
        return lhs->SNPID < rhs->SNPID;
    }
};

/**
 * @brief Functor to delete objects from vectors of pointers
 */
struct DeleteAlleleFreq
{
    template<class T>
    void operator()(const T* p) const {
        delete p;
    }
};


/**
 * @brief  The LOD Familial analysis method.
 *
 */
class CEMSeed;
class CPrimeEM;

class CNFamilialAnalysisMethodLOD : public CNFamilialAnalysisMethod
{
private:
    float m_fRoleValidityThreshold;
    float m_fErrorRate;

    // helper functions
    static double estErrorRate(const std::vector<double>& mean, const std::vector<double>& sigma);
    static void setSeed(CEMSeed& seed);
    static std::pair<double, double> findIntersections(double mean1, double mean2, double sigma1, double sigma2);
    static void fit3CompModel(CPrimeEM& em, const std::vector<double>& vSCAR);

    double estMaxErrorRate(CPrimeEM& indexEM, CPrimeEM& motherEM, CPrimeEM& fatherEM);
    bool genotypeCallsAvailable(CNCychp& cychp);
    void countMIE(
                    const std::vector<char>& vIndexGenotypeCalls,
                    const std::vector<char>& vMotherGenotypeCalls,
                    const std::vector<char>& vFatherGenotypeCalls,
                    const std::vector<unsigned char>& vChromosomes
                    );
    bool isMIE(char genIndex, char genMother, char genFather);
    bool isMIE(char genIndex, char genParent);

public:
    static std::string getType() {return "paternity";}
    static std::string getDescription() {return "Familial LOD";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNFamilialAnalysisMethodLOD();
    virtual ~CNFamilialAnalysisMethodLOD() {}

    virtual bool isFamilialAnalysis() {return true;}
    virtual bool isSegmentOverlapAnalysis() {return false;}
    virtual bool isSegmentTypeAnalysis() {return false;}
    virtual AffxString getName() {return getType();}
    virtual void run();
    virtual bool readFile(const AffxString& strFileName);

    int getAlleleFreqSize() {return m_vAlleleFrequencies.size();}
    const AlleleFrequency* getAlleleFrequency(int i) const {return m_vAlleleFrequencies[i];}

private:
    std::vector<AlleleFrequency*> m_vAlleleFrequencies;
};

typedef std::vector<AlleleFrequency*>::iterator alleleIter;
typedef std::pair<alleleIter, alleleIter> aIterPair_t;


#endif


