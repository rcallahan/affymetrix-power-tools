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
 * @file CNFamilialAnalysisMethodLOD.cpp
 *
 * @brief This file contains the CNFamilialAnalysisMethodLOD class members.
 */
#include "copynumber/CNFamilialAnalysisMethodLOD.h"
//
#include "copynumber/CytogeneticsTrioAnalysis.h"
//
#include "algorithm/em/PrimeEM.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/TsvFile/TsvFile.h"
#include "stats/statfun.h"
#include "util/AffxStatistics.h"
#include "util/AptErrno.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNFamilialAnalysisMethodLOD::explainSelf()
{
    SelfDoc doc;
    doc.setDocName(CNFamilialAnalysisMethodLOD::getType());
    doc.setDocDescription(CNFamilialAnalysisMethodLOD::getDescription());
    doc.setDocOptions(CNFamilialAnalysisMethodLOD::getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNFamilialAnalysisMethodLOD::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  SelfDoc::Opt validityThreshold = {"role-validity-threshold", SelfDoc::Opt::Float, "1000.0", "1000.0", "0", "NA", " Paternity Role Validity Threshold"};
  opts.push_back(validityThreshold);

  SelfDoc::Opt errorRate = {"error-rate", SelfDoc::Opt::Float, "0.01", "0.01", "0.0", "1.0", " Paternity Error Rate"};
  opts.push_back(errorRate);

  return opts;
}

/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate* CNFamilialAnalysisMethodLOD::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNFamilialAnalysisMethodLOD* pMethod = new CNFamilialAnalysisMethodLOD();
    std::string strPrefix = getPrefix();
    pMethod->m_fRoleValidityThreshold = setupFloatParameter("role-validity-threshold", strPrefix, params, doc, opts);
    pMethod->m_fErrorRate = setupFloatParameter("error-rate", strPrefix, params, doc, opts);
    pMethod->setDocName(getType());
    pMethod->setDocDescription(getDescription());
    // specify Calvin groups and datasets this analysis needs to read from cychp file
    pMethod->addGroupDataset("ProbeSets", "CopyNumber");
    pMethod->addGroupDataset("AlgorithmData", "MarkerABSignal");
    pMethod->addGroupDataset("Genotyping", "Calls");
    return pMethod;
}

/**
 * @brief Constructor
 */
CNFamilialAnalysisMethodLOD::CNFamilialAnalysisMethodLOD()
{
    m_fRoleValidityThreshold = 0.0;
    m_fErrorRate = 0.0;
    m_strAlleleFrequencyFileName = "";
}

bool CNFamilialAnalysisMethodLOD::readFile(const AffxString& strFileName) {
  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;
  tsv.setAbortOnError(false);
  
  m_vAlleleFrequencies.reserve(affx::TsvFile::getLineCountInFile(strFileName, false));

  if ( tsv.openTable(strFileName) != affx::TSV_OK ) {
    return false;
  }
  std::string strSnpId;
  std::string strAAllele;
  std::string strBAllele;
  
  while (tsv.nextLevel(0)==affx::TSV_OK) {
    if( tsv.getColumnCount(0) == 3 ) {
      tsv.get(0,0,strSnpId);
      tsv.get(0,1,strAAllele);
      tsv.get(0,2,strBAllele);
      AlleleFrequency* p = new AlleleFrequency;
      p->SNPID = strSnpId;
      p->AAlleleFrequency = AffxByteArray(strAAllele).parseDouble();
      p->BAlleleFrequency = AffxByteArray(strBAllele).parseDouble();
      Verbose::out(4, std::string("strSnpId, strAAllele, strBAllele: ") + strSnpId + ", " +  strAAllele + ", " + strBAllele);
      m_vAlleleFrequencies.push_back(p);
    }
  }
  tsv.clear();
  std::sort(m_vAlleleFrequencies.begin(), m_vAlleleFrequencies.end(), AlleleFreqComp());

  return true;
}

bool CNFamilialAnalysisMethodLOD::genotypeCallsAvailable(CNCychp& cychp)
{
    // Elaborate as needed for future array types
    return cychp.getCychpHeader().getArrayType() == "CytoScanHD_Array";
}

/**
 * @brief Run the analysis.
 */
void CNFamilialAnalysisMethodLOD::run()
{
    Verbose::out(1, "CNFamilialAnalysisMethodLOD::run(...) start");
    Verbose::progressBegin(1, "CNFamilialAnalysisMethodLOD::run(...) ", 5, 1, 5);
    CNCychp& cychpIndex  = getCychpIndex();
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    if (!readFile(m_strAlleleFrequencyFileName))
    {
        Err::errAbort("Cannot open allele frequency file: " + m_strAlleleFrequencyFileName);
    }
    Verbose::progressStep(1);

    const bool motherPresent = cychpMother.getFileName() != "";
    const bool fatherPresent = cychpFather.getFileName() != "";

    bool useGenotypeCalls = genotypeCallsAvailable(cychpIndex);
    double errorRate;

    std::vector<char> vIndexGenotypeCalls;
    std::vector<char> vMotherGenotypeCalls;
    std::vector<char> vFatherGenotypeCalls;
    std::vector<float> vAAlleleFrequencies;
    std::vector<float> vBAlleleFrequencies;
    std::vector<unsigned char> vChromosomes;

    std::vector<char> vIndexAllGenotypeCalls;
    std::vector<char> vMotherAllGenotypeCalls;
    std::vector<char> vFatherAllGenotypeCalls;
    std::vector<unsigned char> vAllChromosomes;

    if (!useGenotypeCalls)
    {
        // Calls not available - try to deduce them from SCAR distributions
        std::vector<double> vIndexSCAR;
        std::vector<double> vMotherSCAR;
        std::vector<double> vFatherSCAR;
        for (int iMarkerIndex = 0; (iMarkerIndex < cychpIndex.getCychpProbeSetsCopyNumbers().getCount()); iMarkerIndex++)
        {
            CNCychpProbeSetsCopyNumber* pMarker = cychpIndex.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex);
            AlleleFrequency objSearch;
            objSearch.SNPID = pMarker->getProbeSetName();
            aIterPair_t range = std::equal_range(m_vAlleleFrequencies.begin(), m_vAlleleFrequencies.end(), &objSearch, AlleleFreqComp());
            if (range.first != range.second)
            {
                vIndexSCAR.push_back(pMarker->getSCAR());
                if (motherPresent)
                {
                    vMotherSCAR.push_back(cychpMother.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getSCAR());
                }
                if (fatherPresent)
                {
                    vFatherSCAR.push_back(cychpFather.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getSCAR());
                }
                vAAlleleFrequencies.push_back((*range.first)->AAlleleFrequency);
                vBAlleleFrequencies.push_back((*range.first)->BAlleleFrequency);
                vChromosomes.push_back(pMarker->getChromosome());
            }
        }
        int iSCARCount = vIndexSCAR.size();
        Verbose::progressStep(1);

        CPrimeEM indexEM;
        CPrimeEM motherEM;
        CPrimeEM fatherEM;

        // 3 component mixture model fit to the SCAR data
        fit3CompModel(indexEM, vIndexSCAR);
        vIndexSCAR.clear();
        std::vector<double>().swap(vIndexSCAR);

        if (motherPresent)
        {
            fit3CompModel(motherEM, vMotherSCAR);
            vMotherSCAR.clear();
            std::vector<double>().swap(vMotherSCAR);
        }
        if (fatherPresent)
        {
            fit3CompModel(fatherEM, vFatherSCAR);
            vFatherSCAR.clear();
            std::vector<double>().swap(vFatherSCAR);
        }

        int gLen = indexEM.getEMEstimates()->m_class.size();
        for (int iSnpIndex = 0; iSnpIndex < gLen; iSnpIndex++)
        {
            vIndexGenotypeCalls.push_back(indexEM.getEMEstimates()->m_class[iSnpIndex]);
            if (motherPresent)
            {
                vMotherGenotypeCalls.push_back(motherEM.getEMEstimates()->m_class[iSnpIndex]);
            }
            if (fatherPresent)
            {
                vFatherGenotypeCalls.push_back(fatherEM.getEMEstimates()->m_class[iSnpIndex]);
            }
        }
        errorRate = estMaxErrorRate(indexEM, motherEM, fatherEM);
    }
    else {
        // Genotype calls available - use them
        for (int iMarkerIndex = 0; iMarkerIndex < cychpIndex.getCychpProbeSetsCopyNumbers().getCount(); iMarkerIndex++)
        {
            CNCychpProbeSetsCopyNumber* pMarker = cychpIndex.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex);
            if (!pMarker->getIsGenotypeCallSet())
            {
                continue;
            }
            vIndexAllGenotypeCalls.push_back(pMarker->getGenotypeCall());
            if (motherPresent)
            {
                vMotherAllGenotypeCalls.push_back(cychpMother.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getGenotypeCall());
            }
            if (fatherPresent)
            {
                vFatherAllGenotypeCalls.push_back(cychpFather.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getGenotypeCall());
            }
            vAllChromosomes.push_back(pMarker->getChromosome());

            AlleleFrequency objSearch;
            objSearch.SNPID = pMarker->getProbeSetName();
            aIterPair_t range = std::equal_range(m_vAlleleFrequencies.begin(), m_vAlleleFrequencies.end(), &objSearch, AlleleFreqComp());
            if (range.first != range.second)
            {
                vIndexGenotypeCalls.push_back(pMarker->getGenotypeCall());
                if (motherPresent)
                {
                    vMotherGenotypeCalls.push_back(cychpMother.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getGenotypeCall());
                }
                if (fatherPresent)
                {
                    vFatherGenotypeCalls.push_back(cychpFather.getCychpProbeSetsCopyNumbers().getAt(iMarkerIndex)->getGenotypeCall());
                }
                vAAlleleFrequencies.push_back((*range.first)->AAlleleFrequency);
                vBAlleleFrequencies.push_back((*range.first)->BAlleleFrequency);
                vChromosomes.push_back(pMarker->getChromosome());
            }
        }
        errorRate = m_fErrorRate;
    }

    std::for_each(m_vAlleleFrequencies.begin(), m_vAlleleFrequencies.end(), DeleteAlleleFreq());
    Verbose::progressStep(1);
    CytogeneticsTrioAnalysis objTrioAnalysis;
    unsigned char ucMotherCall    = 0;
    unsigned char ucFatherCall    = 0;
    unsigned char ucFatherDuoCall = 0;
    double dFatherConfidence    = 0;
    double dMotherConfidence    = 0;
    double dFatherDuoConfidence = 0;

    if (motherPresent && fatherPresent)
    {
        dMotherConfidence    = objTrioAnalysis.calculateLODPaternity(vIndexGenotypeCalls, vMotherGenotypeCalls, vAAlleleFrequencies, vBAlleleFrequencies, errorRate);
        dFatherConfidence    = objTrioAnalysis.calculateLODPaternity(vIndexGenotypeCalls, vMotherGenotypeCalls, vFatherGenotypeCalls, vAAlleleFrequencies, vBAlleleFrequencies, errorRate);
        dFatherDuoConfidence = objTrioAnalysis.calculateLODPaternity(vIndexGenotypeCalls, vFatherGenotypeCalls, vAAlleleFrequencies, vBAlleleFrequencies, errorRate);
        Verbose::progressStep(1);
        if (dMotherConfidence    > m_fRoleValidityThreshold) {ucMotherCall    = (unsigned char)1;}
        if (dFatherConfidence    > m_fRoleValidityThreshold) {ucFatherCall    = (unsigned char)1;}
        if (dFatherDuoConfidence > m_fRoleValidityThreshold) {ucFatherDuoCall = (unsigned char)1;}
    }
    else if (motherPresent)
    {
        dMotherConfidence = objTrioAnalysis.calculateLODPaternity(vIndexGenotypeCalls, vMotherGenotypeCalls, vAAlleleFrequencies, vBAlleleFrequencies, errorRate);
        Verbose::progressStep(1);
        if (dMotherConfidence > m_fRoleValidityThreshold) {ucMotherCall = (unsigned char)1;}
    }
    else if (fatherPresent)
    {
        dFatherDuoConfidence = objTrioAnalysis.calculateLODPaternity(vIndexGenotypeCalls, vFatherGenotypeCalls, vAAlleleFrequencies, vBAlleleFrequencies, errorRate);
        Verbose::progressStep(1);
        if (dFatherDuoConfidence > m_fRoleValidityThreshold) {ucFatherDuoCall = (unsigned char)1;}
    }
    cychpIndex.setFamilialCall((unsigned char)1);
    cychpIndex.setFamilialConfidence((float)1);

    cychpMother.setFamilialCall(ucMotherCall);
    cychpMother.setFamilialConfidence((float)dMotherConfidence);

    cychpFather.setFamilialCall(ucFatherCall);
    cychpFather.setFamilialConfidence((float)dFatherConfidence);

    cychpFather.setFamilialDuoCall(ucFatherDuoCall);
    cychpFather.setFamilialDuoConfidence((float)dFatherDuoConfidence);

    // Mendelian inheritance errors
    if (useGenotypeCalls)
    {
        countMIE(vIndexAllGenotypeCalls, vMotherAllGenotypeCalls, vFatherAllGenotypeCalls, vAllChromosomes);
    }
    else {
        countMIE(vIndexGenotypeCalls, vMotherGenotypeCalls, vFatherGenotypeCalls, vChromosomes);
    }

    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNFamilialAnalysisMethodLOD::run(...) end");
}

void CNFamilialAnalysisMethodLOD::fit3CompModel(CPrimeEM& em, const std::vector<double>& vSCAR)
{
    CEMSeed seed;
    setSeed(seed);
    em.setData(vSCAR);
    em.setThreshold(1.1);    // force no -1 calls
    em.EMEstimate(seed);
}

double CNFamilialAnalysisMethodLOD::estMaxErrorRate(CPrimeEM& indexEM,
                                                   CPrimeEM& motherEM,
                                                   CPrimeEM& fatherEM)
{
    std::vector<double> errorRates;
    errorRates.push_back(estErrorRate(indexEM.getEMEstimates()->m_mu, indexEM.getEMEstimates()->m_sigma));
    if (getCychpMother().getFileName() != "")
    {
        errorRates.push_back(estErrorRate(motherEM.getEMEstimates()->m_mu, motherEM.getEMEstimates()->m_sigma));
    }
    if (getCychpFather().getFileName() != "")
    {
        errorRates.push_back(estErrorRate(fatherEM.getEMEstimates()->m_mu, fatherEM.getEMEstimates()->m_sigma));
    }
    return *std::max_element(errorRates.begin(), errorRates.end());
}

double CNFamilialAnalysisMethodLOD::estErrorRate(const std::vector<double>& mean, const std::vector<double>& sigma)
{
    const double sqrt2 = 1.414213562373;

    double intersect1 = findIntersections(mean[1], mean[2], sigma[1], sigma[2]).first;
    double intersect2 = findIntersections(mean[1], mean[0], sigma[1], sigma[0]).second;

    double erfArg1 = (intersect1 - mean[2])/(sigma[2]*sqrt2);
    double erfArg2 = (intersect2 - mean[0])/(sigma[0]*sqrt2);

    double errorRate1 = 0.5 + 0.5*affxstat::erf(erfArg1);
    double errorRate2 = 0.5 - 0.5*affxstat::erf(erfArg2);

    return errorRate1 > errorRate2 ? errorRate1 : errorRate2;
}

void CNFamilialAnalysisMethodLOD::setSeed(CEMSeed& seed)
{
    seed.setMu(-1.0, 0, 1.0);
    seed.setSigma(0.1f, 0.1f, 0.1f);
    seed.setWeight(0.33f, 0.34f, 0.33f);
    seed.setMinMu(-2, -.05, .25);
    seed.setMinSigma(0.02f, 0.02f, 0.02f);
    seed.setMaxMu(-.25, .05, 2);
    seed.setMaxSigma(.3, .3, .3);
}

// Find intersection points of given gaussians
std::pair<double, double> CNFamilialAnalysisMethodLOD::findIntersections(double mean1, double mean2, double sigma1, double sigma2)
{
    using std::abs;
    using std::log;
    using std::sqrt;

    const double epsilon = 1.0e-10;

    double intersect1;
    double intersect2;

    if (abs(sigma1 - sigma2) < epsilon)
    {
        intersect1 = intersect2 = 0.5*(mean1 + mean2);
    }
    else {
        double A = sigma1*sigma1 - sigma2*sigma2;
        double B = -2.0*(mean2*sigma1*sigma1 - mean1*sigma2*sigma2);
        double C = mean2*mean2*sigma1*sigma1 - mean1*mean1*sigma2*sigma2 - 2.0*sigma2*sigma2*sigma1*sigma1*(log(sigma2) - log(sigma1));
        intersect1 = (-B - sqrt(B*B - 4.0*A*C))/(2.0*A);
        intersect2 = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
    }
    return std::make_pair(intersect1, intersect2);
}

void CNFamilialAnalysisMethodLOD::countMIE(
                                    const std::vector<char>& vIndexGenotypeCalls,
                                    const std::vector<char>& vMotherGenotypeCalls,
                                    const std::vector<char>& vFatherGenotypeCalls,
                                    const std::vector<unsigned char>& vChromosomes
                                    )
{
    const bool trioCase = !vMotherGenotypeCalls.empty() && !vFatherGenotypeCalls.empty();
    const bool duoMat   = !vMotherGenotypeCalls.empty();
    const bool duoPat   = !vFatherGenotypeCalls.empty();

    if (trioCase)
    {
        for (int iIndex = 0; iIndex < vIndexGenotypeCalls.size(); iIndex++)
        {
            CNCychp::m_mMarkerCount[vChromosomes[iIndex]]++;
            if (isMIE(vIndexGenotypeCalls[iIndex], vMotherGenotypeCalls[iIndex], vFatherGenotypeCalls[iIndex]))
            {
                CNCychp::m_mMIE_Trio[vChromosomes[iIndex]]++;
            }
            if (isMIE(vIndexGenotypeCalls[iIndex], vMotherGenotypeCalls[iIndex]))
            {
                CNCychp::m_mMIE_Mat[vChromosomes[iIndex]]++;
            }
            if (isMIE(vIndexGenotypeCalls[iIndex], vFatherGenotypeCalls[iIndex]))
            {
                CNCychp::m_mMIE_Pat[vChromosomes[iIndex]]++;
            }
        }
    }
    else if (duoMat)
    {
        for (int iIndex = 0; iIndex < vIndexGenotypeCalls.size(); iIndex++)
        {
            CNCychp::m_mMarkerCount[vChromosomes[iIndex]]++;
            if (isMIE(vIndexGenotypeCalls[iIndex], vMotherGenotypeCalls[iIndex]))
            {
                CNCychp::m_mMIE_Mat[vChromosomes[iIndex]]++;
            }
        }
    }
    else if (duoPat)
    {
        for (int iIndex = 0; iIndex < vIndexGenotypeCalls.size(); iIndex++)
        {
            CNCychp::m_mMarkerCount[vChromosomes[iIndex]]++;
            if (isMIE(vIndexGenotypeCalls[iIndex], vFatherGenotypeCalls[iIndex]))
            {
                CNCychp::m_mMIE_Pat[vChromosomes[iIndex]]++;
            }
        }
    }
}

inline bool CNFamilialAnalysisMethodLOD::isMIE(char genIndex, char genMother, char genFather)
{
    // Mendelian error if either:
    //  Mother  Father  Child
    //  AA      AA      AB|BB
    //  AA      BB      AA|BB
    //  AA      AB      BB
    //  AB      AA      BB
    //  AB      BB      AA
    //  BB      AB      AA
    //  BB      AA      AA|BB
    //  BB      BB      AA|AB
    //
    return
        genMother == 0 && genFather == 0 && (genIndex == 1 || genIndex == 2)  ||
        genMother == 0 && genFather == 2 && (genIndex == 0 || genIndex == 2)  ||
        genMother == 0 && genFather == 1 && genIndex == 2                     ||
        genMother == 1 && genFather == 0 && genIndex == 2                     ||
        genMother == 1 && genFather == 2 && genIndex == 0                     ||
        genMother == 2 && genFather == 1 && genIndex == 0                     ||
        genMother == 2 && genFather == 0 && (genIndex == 0 || genIndex == 2)  ||
        genMother == 2 && genFather == 2 && (genIndex == 0 || genIndex == 1);
}

inline bool CNFamilialAnalysisMethodLOD::isMIE(char genIndex, char genParent)
{
    // Mendelian error if either:
    //  Parent  Child
    //  AA      BB
    //  BB      AA
    //
    return
        genParent == 0 && genIndex == 2  ||
        genParent == 2 && genIndex == 0;
}
