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
 * @file CNAnalysisMethodLOHCyto2.cpp
 *
 * @brief This file contains the CNAnalysisMethodLOHCyto2 class members.
 */

//
#include "copynumber/CNAnalysisMethodLOHCyto2.h"
//
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
#include "stats/stats-distributions.h"
#include "util/AffxStatistics.h"
#include "util/Fs.h"
#include "util/Guid.h"
//
#include <cmath>
#include <limits>
//
using namespace std;

using namespace affxstat;
using namespace affx;

// constructor
CNAnalysisMethodLOHCyto2::CNAnalysisMethodLOHCyto2() {
    m_bSNPAlreadyLoaded=false;
    m_iLohCNSegSeparation=0;

    m_fMuPrimeAA=0.0;
    m_fMuPrimeAB=0.0;
    m_fMuPrimeBB=0.0;

    m_fSigmaPrimeAA=0.0;
    m_fSigmaPrimeAB=0.0;
    m_fSigmaPrimeBB=0.0;

    m_fMinInformation=0.0;
    m_fMinInformationCommandLineValue=0.0;

    m_fEMSN=0.0;

    m_fLambdaCritical=0.0;
    m_fLambdaCriticalCommandLineValue=0.0;
    m_fLambdaSpacing=0.0;
    m_fLambdaStartingValue=0.0;
    m_iNumberOfLambdaPoints=0;

    m_fEMSN=0.0;
    m_fEMSNSpacing=0.0;
    m_fEMSNStartingValue=0.0;
    m_iNumberOfEMSNPoints=0;
}

CNAnalysisMethodLOHCyto2::~CNAnalysisMethodLOHCyto2() {}

std::string CNAnalysisMethodLOHCyto2::getType() {
    return "loh-cyto2";
}
std::string CNAnalysisMethodLOHCyto2::getDescription() {
    return "CopyNumber LOH Cyto2";
}
std::string CNAnalysisMethodLOHCyto2::getVersion() {
    return CN_ANALYSIS_METHOD_LOH_CYTO2_VERSION;
}

AffxString CNAnalysisMethodLOHCyto2::getName() {
    return getType();
}

bool CNAnalysisMethodLOHCyto2::isSegmentTypeAnalysis() {
    return true;
}

/*
inline double normalDistribution(double ecks, double mu, double sigma)
{
  double tau = 1.0/(sigma*sigma);
  double sqrt_tau = 1.0/sigma;
  double diff1 = ecks - mu;
  double diff2 = tau*diff1*diff1;
  return exp(-0.5*diff2);
}
*/

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodLOHCyto2::explainSelf()
{
    CNAnalysisMethodLOHCyto2 obj;
    SelfDoc doc;
    doc.setDocName(obj.getType());
    doc.setDocDescription(obj.getDescription());
    doc.setDocOptions(obj.getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNAnalysisMethodLOHCyto2::getDefaultDocOptions()
{
    std::vector<SelfDoc::Opt> opts;

    // SelfDoc::Opt(name, type, value, default, min, max, description)

    SelfDoc::Opt lohCNSegSeparation = {"lohCNSegSeparation", SelfDoc::Opt::Integer, "1000000", "1000000", "10", "10000000", "LOH CN Separation"};
    opts.push_back(lohCNSegSeparation);

    SelfDoc::Opt minInformation = {"minInformation", SelfDoc::Opt::Float, "100", "100", "NA", "NA", "A control of window size." };
    opts.push_back(minInformation);

    SelfDoc::Opt lambdaCritical = {"lambdaCritical", SelfDoc::Opt::Float, "8.0", "8.0", "NA", "NA", "A measure of required likelihood." };
    opts.push_back(lambdaCritical);

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
SelfCreate* CNAnalysisMethodLOHCyto2::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodLOHCyto2* pMethod = new CNAnalysisMethodLOHCyto2();
    std::string strPrefix = getPrefix();

    pMethod->m_iLohCNSegSeparation = setupIntParameter("lohCNSegSeparation", strPrefix, params, doc, opts);

    pMethod->m_fMinInformationCommandLineValue = setupFloatParameter("minInformation", strPrefix, params, doc, opts);
    pMethod->m_fLambdaCriticalCommandLineValue = setupFloatParameter("lambdaCritical", strPrefix, params, doc, opts);

    return pMethod;
}

void CNAnalysisMethodLOHCyto2::reset()
{
    //  There are numerical comparison on these values so they must be reset after each sample is analyzed.
    m_fMinInformation=m_fMinInformationCommandLineValue;
    m_fLambdaCritical=m_fLambdaCriticalCommandLineValue;
}


/**
 * @brief run the analysis
 */
void CNAnalysisMethodLOHCyto2::run()
{
    Verbose::out(1, "CNAnalysisMethodLOHCyto2::run(...) start");
    isSetup();
    reset();

    if(!m_bSNPAlreadyLoaded)
    {
        ///@todo  Add error messaging if no snp-reference file is input. Also exit run gracefully.
        AffxString strFileName = m_pEngine->getOpt("snp-reference-input-file");
        loadSnpReferenceFile(strFileName);
        m_bSNPAlreadyLoaded=true;

/* The following are the values presently read from the snp-reference file. They are here as a quick reference.

Full array values

        m_vVarianceInflationValues.push_back(1.423082);
        m_vVarianceInflationValues.push_back(1.348103);
        m_vVarianceInflationValues.push_back(1.273123);
        m_vVarianceInflationValues.push_back(1.198369);
        m_vVarianceInflationValues.push_back(1.148516);
        m_vVarianceInflationValues.push_back(1.098662);
        m_vVarianceInflationValues.push_back(1.0);

        m_vInformationRanges.push_back(0.4061303);
        m_vInformationRanges.push_back(0.4733475);
        m_vInformationRanges.push_back(0.5923077);
        m_vInformationRanges.push_back(0.7761194);
        m_vInformationRanges.push_back(0.9158621);
        m_vInformationRanges.push_back(1.1380471);

Focused array values
        m_vVarianceInflationValues.push_back(1.4);
        m_vVarianceInflationValues.push_back(1.2);
        m_vVarianceInflationValues.push_back(1.1);
        m_vVarianceInflationValues.push_back(1.0);


        m_vInformationRanges.push_back(2.8);
        m_vInformationRanges.push_back(3.0);
        m_vInformationRanges.push_back(3.2);
*/


        m_iNumberOfVIValues=m_vVarianceInflationValues.size();
    }


    std::vector<float> vSCAR(getProbeSets()->size());
    std::vector<float> vInformation(getProbeSets()->size());
    std::vector<int> vPosition(getProbeSets()->size());
    AffxMultiDimensionalArray<double> dHomValueArray;
    AffxMultiDimensionalArray<double> dHetValueArray;
    dHomValueArray.initialize(getProbeSets()->size(), m_iNumberOfVIValues);
    dHetValueArray.initialize(getProbeSets()->size(), m_iNumberOfVIValues);

    // Note that the vPosition is not used anywhere within the LOH calculation. It is placed here to allow comparison against external R
    // algorithms which generate beginning and ending positions of windows over the genome.  This windowing is done here in the function lohFind.
    // For comparison to these external algorithms, position information is passed to the function lohFind.

    if (!lohPreProcessing(vSCAR, vInformation, vPosition, dHetValueArray, dHomValueArray )) {throw(Except("LohPreProcessing failed."));}

    // storage for intermediate data collected for validation
    int iUseProbeSetCount = 0;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++) {
        if ( getProbeSets()->at(iIndex)->isUseForCyto2LOH() ) {
            iUseProbeSetCount++;
        }
    }
    std::vector<int> vWindowStart;
    std::vector<int> vWindowLength;
    std::vector<double> vWindowTest;
    std::vector<double> vWindowInflate;
    std::vector<char> vWindowLoh;
    vWindowStart.reserve(iUseProbeSetCount);
    vWindowLength.reserve(iUseProbeSetCount);
    vWindowTest.reserve(iUseProbeSetCount);
    vWindowInflate.reserve(iUseProbeSetCount);
    vWindowLoh.reserve(iUseProbeSetCount);
    std::vector<int> vGenomeChunkStart;
    std::vector<int> vGenomeChunkLength;

    int iLastChromosome = getProbeSets()->at(getProbeSets()->size() - 1)->getChromosome();
    Verbose::progressBegin(1, "CNAnalysisMethodLOHCyto2::run(...) ", iLastChromosome, 1, iLastChromosome);
    for (int iChromosome = 1; (iChromosome <= iLastChromosome); iChromosome++)
    {
        try
        {
            if (iChromosome > m_iYChromosome) {continue;}
            int iChromosomeProbeSetCount = 0;
            for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
                if (pobjProbeSet->getChromosome() == iChromosome && pobjProbeSet->isUseForCyto2LOH() )
                {
                    iChromosomeProbeSetCount++;
                }
            }

            if (iChromosomeProbeSetCount == 0) {continue;}
            Verbose::progressStep(1);

            std::vector<float> vChromosomeInformation(iChromosomeProbeSetCount);
            std::vector<float> vChromosomeSCAR(iChromosomeProbeSetCount);
            std::vector<int> vChromosomeProbeSetIndexes(iChromosomeProbeSetCount);
            std::vector<int> vChromosomePosition(iChromosomeProbeSetCount);
            AffxMultiDimensionalArray<double> dChromosomeHomValueArray;
            AffxMultiDimensionalArray<double> dChromosomeHetValueArray;
            dChromosomeHomValueArray.initialize(iChromosomeProbeSetCount, m_iNumberOfVIValues);
            dChromosomeHetValueArray.initialize(iChromosomeProbeSetCount, m_iNumberOfVIValues);

            int iProbeSetIndex = 0;
            for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);

                if ( pobjProbeSet->getChromosome() == iChromosome && pobjProbeSet->isUseForCyto2LOH() )
                {
                    //  This is the location that we see the necessity of having vSCARreduced and vSCAR.  We
                    //  cannot, at this point in the code have vSCAR containing only SCAR values for probes
                    //  with valid SCARs and HomHet values since we need to access them via the index of the
                    //  "getProbeSets" vector.
                    vChromosomeSCAR[iProbeSetIndex] = vSCAR[iIndex];
                    vChromosomeInformation[iProbeSetIndex] = vInformation[iIndex];
                    vChromosomeProbeSetIndexes[iProbeSetIndex] = iIndex;
                    vChromosomePosition[iProbeSetIndex] = vPosition[iIndex];

                    for(int iValueIndex=0; iValueIndex<m_iNumberOfVIValues; iValueIndex++)
                    {
                            dChromosomeHomValueArray.set(iProbeSetIndex, iValueIndex, dHomValueArray.get(iIndex, iValueIndex));
                            dChromosomeHetValueArray.set(iProbeSetIndex, iValueIndex, dHetValueArray.get(iIndex, iValueIndex));
                    }
                    iProbeSetIndex++;
                }
            }
            int iStartIndex = 0;
            for (int iIndex = 1; (iIndex <= iChromosomeProbeSetCount); iIndex++)
            {
                //  Note that iIndex may equal iChromosomeProbeSetCount but this will not cause a problem in the access
                //  to the vector vChromosomeProbeSetIndex in the next line since the first part of the  || will fail first.
                if ((iIndex == iChromosomeProbeSetCount) || ((getProbeSets()->at(vChromosomeProbeSetIndexes[iIndex])->getPosition() - getProbeSets()->at(vChromosomeProbeSetIndexes[iIndex - 1])->getPosition()) > m_iLohCNSegSeparation))
                {
                    int iSegmentProbeSetCount = ((iIndex - 1) - iStartIndex + 1);
                    std::vector<float> vSegmentSCAR(iSegmentProbeSetCount);
                    std::vector<float> vSegmentInformation(iSegmentProbeSetCount);
                    std::vector<char> vLoh(iSegmentProbeSetCount, 0);
                    std::vector<int> vSegmentPosition(iSegmentProbeSetCount);
                    AffxMultiDimensionalArray<double> dSegmentHomValueArray;
                    AffxMultiDimensionalArray<double> dSegmentHetValueArray;
                    dSegmentHomValueArray.initialize(iSegmentProbeSetCount, m_iNumberOfVIValues);
                    dSegmentHetValueArray.initialize(iSegmentProbeSetCount, m_iNumberOfVIValues);

                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentProbeSetCount); iSegmentIndex++)
                    {

                        vSegmentPosition[iSegmentIndex] = vChromosomePosition[iStartIndex + iSegmentIndex];
                        vSegmentSCAR[iSegmentIndex] = vChromosomeSCAR[iStartIndex + iSegmentIndex];
                        vSegmentInformation[iSegmentIndex] = vChromosomeInformation[iStartIndex + iSegmentIndex];
                        for(int iValueIndex=0; iValueIndex<m_iNumberOfVIValues; iValueIndex++)
                        {
                            dSegmentHomValueArray.set(iSegmentIndex, iValueIndex, dChromosomeHomValueArray.get(iStartIndex+iSegmentIndex, iValueIndex));
                            dSegmentHetValueArray.set(iSegmentIndex, iValueIndex, dChromosomeHetValueArray.get(iStartIndex+iSegmentIndex, iValueIndex));
                        }
                    }
                    std::vector<int> vSegmentWindowStart(iSegmentProbeSetCount, -1);
                    std::vector<int> vSegmentWindowLength(iSegmentProbeSetCount, -1);
                    std::vector<double> vSegmentWindowTest(iSegmentProbeSetCount, numeric_limits<float>::quiet_NaN() );
                    std::vector<double> vSegmentWindowInflate(iSegmentProbeSetCount, numeric_limits<float>::quiet_NaN() );
                    std::vector<char> vSegmentWindowLoh(iSegmentProbeSetCount, 0);
                    lohFind(vSegmentSCAR, vSegmentInformation, vSegmentPosition, vLoh, dSegmentHomValueArray, dSegmentHetValueArray, vSegmentWindowStart, vSegmentWindowLength, vSegmentWindowTest, vSegmentWindowInflate, vSegmentWindowLoh);
                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentProbeSetCount); iSegmentIndex++)
                    {
                        getProbeSets()->at(vChromosomeProbeSetIndexes[iStartIndex + iSegmentIndex])->setLoh(vLoh[iSegmentIndex]);
                        if( vSegmentWindowStart[iSegmentIndex] >= 0 ) {
                            vWindowStart.push_back( vChromosomeProbeSetIndexes[ iStartIndex + vSegmentWindowStart[iSegmentIndex] ] );
                            vWindowLength.push_back( vSegmentWindowLength[iSegmentIndex] );
                            vWindowTest.push_back( vSegmentWindowTest[iSegmentIndex] );
                            vWindowInflate.push_back( vSegmentWindowInflate[iSegmentIndex] );
                            vWindowLoh.push_back( vSegmentWindowLoh[iSegmentIndex] );
                        }
                    }
                    vGenomeChunkStart.push_back( vChromosomeProbeSetIndexes[iStartIndex] );
                    vGenomeChunkLength.push_back( iSegmentProbeSetCount );
                    iStartIndex = iIndex;
                }
            }
        }
        catch(...) {
            Verbose::out(1, "CNAnalysisMethodLOHCyto2::run(...) Working on Chromosome " + ::getInt(iChromosome));
            throw;
        }
    }

    m_vSegments.deleteAll();
    std::vector<double> vHomValues;
    std::vector<double> vHetValues;

    for(int iIndex=0; iIndex < getProbeSets()->size();iIndex++)
    {
        vHomValues.push_back(dHomValueArray.get(iIndex,m_iNumberOfVIValues-1));
        vHetValues.push_back(dHetValueArray.get(iIndex,m_iNumberOfVIValues-1));
    }

    newSegments(getSegmentType(), vSCAR, vHomValues, vHetValues);
    calculateSummaryLOH();
    newSegments(4, vSCAR, vHomValues, vHetValues);
    newSegments(5, vSCAR, vHomValues, vHetValues);

    if (m_pEngine->getOptBool("keep-intermediate-data"))
    {
        writeLOH( vSCAR, vGenomeChunkStart, vGenomeChunkLength, vWindowStart, vWindowLength, vWindowTest, vWindowInflate, vWindowLoh );
    }

    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodLOHCyto2::run(...) end");
}


void CNAnalysisMethodLOHCyto2::computeGlobalSnpParameters()
{

    // At this point a value for m_fEMSN has been computed and the value of m_fLambdaCritical gotten from command line.
    // If these values fall out of the ranges of the lookup table we set them to the boundary values of the table.
    int iEMSNIndex = 0;
    int iLambdaIndex = 0;

    float fMaxEMSNValue = m_fEMSNStartingValue + ((m_iNumberOfEMSNPoints-1) * m_fEMSNSpacing);
    if(m_fEMSN > fMaxEMSNValue)
    {
        iEMSNIndex = m_iNumberOfEMSNPoints-1;
    } else
    {
        if(m_fEMSN < m_fEMSNStartingValue)
        {
            iEMSNIndex=0;
        } else
        {
            float fEMSNFloatIndex = (m_fEMSN - m_fEMSNStartingValue) / m_fEMSNSpacing;
            double dEMSNTruncated=0.0;
            double dEMSNRemainder = modf((double)fEMSNFloatIndex, &dEMSNTruncated);
            if(dEMSNRemainder > 0.5) dEMSNTruncated = dEMSNTruncated + 1.0;
            iEMSNIndex = (int) dEMSNTruncated;
        }
    }

    float fMaxLambdaValue = m_fLambdaStartingValue + ((m_iNumberOfLambdaPoints-1) * m_fLambdaSpacing);
    if(m_fLambdaCritical > fMaxLambdaValue)
    {
        iLambdaIndex = m_iNumberOfLambdaPoints-1;
    } else
    {
        if(m_fLambdaCritical < m_fLambdaStartingValue)
        {
            iLambdaIndex=0;
        } else
        {
            float fLambdaFloatIndex = (m_fLambdaCritical - m_fLambdaStartingValue) / m_fLambdaSpacing;
            double dLambdaTruncated=0.0;
            double dLambdaRemainder = modf((double)fLambdaFloatIndex, &dLambdaTruncated);
            if(dLambdaRemainder > 0.5) dLambdaTruncated = dLambdaTruncated + 1.0;
            iLambdaIndex = (int) dLambdaTruncated;
        }
    }

    float fTmin = m_dEMSNArray.get(iEMSNIndex, iLambdaIndex);
    if(m_fMinInformation < fTmin) m_fMinInformation = fTmin;
    Verbose::out(3, "CNAnalysisMethodLOHCyto2::minInformation value is: " + AffxString::doubleToString((double)m_fMinInformation, 6, false));
}



bool CNAnalysisMethodLOHCyto2::loadSnpReferenceFile(const AffxString& strFileName)
{
    Verbose::out(3, "CNAnalysisMethodLOHCyto2::loadSnpReferenceFile");
    bool bSuccessful = false;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5;

            tsv5 = group5->openTsv("EMSNParameters", affx::FILE5_OPEN);
            double dInputDouble=0.0;
            int iInputInteger=0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &dInputDouble);  m_fLambdaSpacing=(float)dInputDouble;
                tsv5->get(0, 1, &dInputDouble);  m_fEMSNSpacing=(float)dInputDouble;
                tsv5->get(0, 2, &dInputDouble);  m_fLambdaStartingValue=(float)dInputDouble;
                tsv5->get(0, 3, &dInputDouble);  m_fEMSNStartingValue=(float)dInputDouble;
                tsv5->get(0, 4, &iInputInteger); m_iNumberOfLambdaPoints=iInputInteger;
                tsv5->get(0, 5, &iInputInteger); m_iNumberOfEMSNPoints=iInputInteger;
            }
            tsv5->close();
            delete tsv5;


            m_dEMSNArray.initialize(m_iNumberOfEMSNPoints, m_iNumberOfLambdaPoints);
            tsv5 = group5->openTsv("EMSNArray", affx::FILE5_OPEN);
            dInputDouble=0.0;
            int iRowCount=0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                for(int iIndex=0; iIndex < m_iNumberOfLambdaPoints; iIndex++)
                {
                    tsv5->get(0, iIndex, &dInputDouble);
                    m_dEMSNArray.set(iRowCount, iIndex, dInputDouble);
                }
                iRowCount++;
            }
            tsv5->close();
            delete tsv5;


            if(group5->name_exists("VIRanges") && group5->name_exists("VIValues"))
            {
                tsv5 = group5->openTsv("VIRanges", affx::FILE5_OPEN);
                if(tsv5)
                {
                    dInputDouble=0.0;
                    while (tsv5->nextLine() == affx::FILE5_OK)
                    {
                        tsv5->get(0, 0, &dInputDouble);
                        m_vInformationRanges.push_back(dInputDouble);
                    }
                }
                tsv5->close();
                delete tsv5;

                tsv5 = group5->openTsv("VIValues", affx::FILE5_OPEN);
                if(tsv5)
                {
                    dInputDouble=0.0;
                    while (tsv5->nextLine() == affx::FILE5_OK)
                    {
                        tsv5->get(0, 0, &dInputDouble);
                        m_vVarianceInflationValues.push_back(dInputDouble);
                    }
                    tsv5->close();
                    delete tsv5;
                }
            }
            else
            {
                m_vVarianceInflationValues.push_back(1.0);
            }

            group5->close();
            delete group5;
            file5.close();
            bSuccessful = true;

        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }
    return bSuccessful;
}

float CNAnalysisMethodLOHCyto2::logLikelihood(  std::vector<float>& vfEcks,
                                                float fMu,
                                                float fSigma,
                                                float fDelta,
                                                float fProbability1,
                                                float fProbability2,
                                                float fProbability3)
{
    float fLogLikelihood=0.0;
    int iNumberOfProbesets= vfEcks.size();
    for(int iIndex=0; iIndex<iNumberOfProbesets;iIndex++)
    {
        float fTermInSum=0.0;
        fTermInSum += fProbability1 * normalDistribution( (double)(vfEcks[iIndex]), (double)(fMu+fDelta), (double)fSigma);
        fTermInSum += fProbability2 * normalDistribution( (double)(vfEcks[iIndex]), (double)fMu, (double)(0.92 * fSigma));
        fTermInSum += fProbability3 * normalDistribution( (double)(vfEcks[iIndex]), (double)(fMu-fDelta), (double)fSigma);
        fLogLikelihood+=log(fTermInSum);
    }

    return fLogLikelihood;
}

bool CNAnalysisMethodLOHCyto2::computeGlobalParameters(    const std::vector<float>& vSCAR,
                                                        double & fMu,
                                                        double & fSigma,
                                                        double& fDelta,
                                                        double& fProbability1,
                                                        double& fProbability2,
                                                        double& fProbability3)
{

    // Just a bit of Canada....
    std::vector<double> vfZed1;
    std::vector<double> vfZed2;
    std::vector<double> vfZed3;

    std::vector<float> vSCARreduced;

    int iNumberOfProbeSets = vSCAR.size();
    for (int iIndex = 0; iIndex < iNumberOfProbeSets; iIndex++)
    {

        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iIndex);
        if(pobjProbeSet->getValidSCARExists() && pobjProbeSet->getUseInEMAlgorithm() == 1)
        {
             if(pobjProbeSet->getInformation() == 0)
             {
                 int i=0;
                 i++;
             }
        }


        if(pobjProbeSet->getValidSCARExists() && pobjProbeSet->getUseInEMAlgorithm() == 1 && pobjProbeSet->getInformation() !=0)
        {
            vSCARreduced.push_back(vSCAR[iIndex]);
        }
    }

    iNumberOfProbeSets=vSCARreduced.size();
    for(int iIndex=0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        vfZed1.push_back( fProbability1 * normalDistribution((double)(vSCARreduced[iIndex]), (double)(fMu+fDelta), (double) fSigma) );
        vfZed2.push_back( fProbability2 * normalDistribution((double)(vSCARreduced[iIndex]), (double) fMu, (double)(0.92*fSigma)) );
        vfZed3.push_back( fProbability3 * normalDistribution((double)(vSCARreduced[iIndex]), (double)(fMu-fDelta), (double)fSigma) );
    }

    bool bPrecisionCheckPassed=false;
    double fPrecisionOld=0.0;
    bPrecisionCheckPassed=sumLogsInVector(vfZed1, vfZed2, vfZed3,fPrecisionOld);
    if(!bPrecisionCheckPassed)
    {
        //  We expect that the failure of computeGlobalParameters will be dealt with by the calling function. i.e Error messaging.
        return false;
    }

    double fPrecisionNew = fPrecisionOld+1;;
    int iLoopIndex=1;
    double fConstant = 0.92 * 0.92;



/*      This is code to output intermediate values computed in the EM algorithm.
    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->writeTsv(getExperiment()->getExperimentName() + ".EM.txt");
    tsv->defineColumn(      0,0,"mu", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,1,"delta", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,2,"sigma", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,3,"pi1", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,4,"pi2", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,5,"pi3", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,6,"precision", affx::TSV_TYPE_DOUBLE);
    tsv->defineColumn(      0,7,"loopIndex", affx::TSV_TYPE_INT);
*/
    while(iLoopIndex < 100 && fPrecisionNew - fPrecisionOld > 0.00005)
    {
        iLoopIndex++;
        fPrecisionOld=fPrecisionNew;
        double fASum=0.0;
        double fBSum=0.0;
        double fC1Sum=0.0;
        double fC2Sum=0.0;
        double fD1Sum=0.0;

        double fZed1Term=0.0;
        double fZed2Term=0.0;
        double fZed3Term=0.0;

        double fZed1Sum=0.0;
        double fZed2Sum=0.0;
        double fZed3Sum=0.0;

        for(int iIndex=0; iIndex<iNumberOfProbeSets; iIndex++)
        {
            double fSum= vfZed1[iIndex] + vfZed2[iIndex] + vfZed3[iIndex];
            fZed1Term = vfZed1[iIndex]/fSum;
            fZed2Term = vfZed2[iIndex]/fSum;
            fZed3Term = vfZed3[iIndex]/fSum;

            vfZed1[iIndex] = fZed1Term;
            vfZed2[iIndex] = fZed2Term;
            vfZed3[iIndex] = fZed3Term;

            fZed1Sum+=fZed1Term;
            fZed2Sum+=fZed2Term;
            fZed3Sum+=fZed3Term;

            fASum += ( fZed1Term - fZed3Term );
            fBSum += ( fZed1Term + fZed3Term );
            fC2Sum += ( (fZed1Term - fZed3Term ) * (double) vSCARreduced[iIndex]);

            double fTemporaryValue= fConstant*(fZed1Term  +fZed3Term) + fZed2Term;
            fD1Sum += fTemporaryValue ;
            fC1Sum += (fTemporaryValue * (double) vSCARreduced[iIndex]);
        }

        fMu=((fBSum* fC1Sum) + (fConstant*fASum*fC2Sum)) / ( (fBSum*fD1Sum)-(fConstant*fASum*fASum) );
        fDelta =  fMu*(fASum/fBSum) + (fC2Sum/fBSum);

        fSigma=0.0;
        for(int iIndex=0; iIndex<iNumberOfProbeSets; iIndex++)
        {
            fSigma+=    (

                            ( vfZed1[iIndex]*( ((double)vSCARreduced[iIndex] - (fMu+fDelta)) * ((double) vSCARreduced[iIndex] -(fMu+fDelta)) )  )
                             +
                            ( (1.0/fConstant)*vfZed2[iIndex]*( ((double)vSCARreduced[iIndex] -fMu)*((double)vSCARreduced[iIndex] -fMu) )  )
                             +
                            ( vfZed3[iIndex]*( ((double)vSCARreduced[iIndex] - (fMu-fDelta)) * ((double)vSCARreduced[iIndex] -(fMu-fDelta)) )  )

                            );
        }

        fSigma=sqrt(fSigma/iNumberOfProbeSets);

        fProbability1=fZed1Sum/iNumberOfProbeSets;
        fProbability2=fZed2Sum/iNumberOfProbeSets;
        fProbability3=fZed3Sum/iNumberOfProbeSets;

        for(int iIndex=0; iIndex<iNumberOfProbeSets; iIndex++)
        {
            vfZed1[iIndex]=fProbability1*normalDistribution((double)(vSCARreduced[iIndex]),double(fMu+fDelta),(double)fSigma);
            vfZed2[iIndex]=fProbability2*normalDistribution((double)(vSCARreduced[iIndex]), (double)(fMu),(double)(0.92*fSigma));
            vfZed3[iIndex]=fProbability3*normalDistribution((double)(vSCARreduced[iIndex]), (double)(fMu-fDelta),(double)fSigma);
        }
        fPrecisionNew=0.0;
        bPrecisionCheckPassed=sumLogsInVector(vfZed1, vfZed2, vfZed3, fPrecisionNew);
        if(!bPrecisionCheckPassed){
            //  We expect that the failure of computeGlobalParameters will be dealt with by the calling function. i.e Error messaging.
            return false;
        }

/*              This is code to output intermediate values in the EM algorithm.
        bool testFlag = false;
        if(testFlag)
        {
            tsv->set(0,0,fMu);
            tsv->set(0,1,fDelta);
            tsv->set(0,2,fSigma);
            tsv->set(0,3,fProbability1);
            tsv->set(0,4,fProbability2);
            tsv->set(0,5,fProbability3);
            tsv->set(0,6,fPrecisionNew);
            tsv->set(0,7,iLoopIndex);
            tsv->writeLevel(0);
        }
*/
    }
/*      This is code to output intermediate values in the EM algorithm.
    delete tsv;
*/
    if(iLoopIndex==100)
        return false;
    else
        return true;
}

bool CNAnalysisMethodLOHCyto2::sumLogsInVector(            const std::vector<double>& vfZed1,
                                                        const std::vector<double>& vfZed2,
                                                        const std::vector<double>& vfZed3,
                                                        double& fSumOfLogs)
{
    int iNumberOfProbeSets = vfZed1.size();
    float fTermInSum=0.0;
    for(int iIndex=0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        fTermInSum=0.0;
        fTermInSum+=vfZed1[iIndex];
        fTermInSum+=vfZed2[iIndex];
        fTermInSum+=vfZed3[iIndex];
        if(fTermInSum<=0.0)
        {
            Verbose::warn(1, "CNAnalysisMethodLOHCyto2::Expectation-Minimization algorithm did not converge.  Invalid terms were computed which checking precision of the computed parameters.  Global snp parameters will be taken from command line or snp-reference-input-file");
            return false;
        }
        fSumOfLogs += log(fTermInSum);
    }
    return true;
}



/**
 * @brief Preprocess the log2Ratio and probeset parameters to procude the SCAR values.
 * @return bool - Was the function successful or not (true = successful)
 */
bool CNAnalysisMethodLOHCyto2::lohPreProcessing(        std::vector<float>& vSCAR,
                                                        std::vector<float>& vInformation,
                                                        std::vector<int>& vPosition,
                                                        AffxMultiDimensionalArray<double>& dHetValueArray,
                                                        AffxMultiDimensionalArray<double>& dHomValueArray )
{
    loadInformationVector(vInformation);

    loadPositionVector(vPosition);

    computeSCAR(vSCAR);

    double fMu=0.0;
    double fSigma=0.5;
    double fDelta=1.0;
    double fProbability1=0.34;
    double fProbability2=0.32;
    double fProbability3=0.34;

    bool bGP=false;
    bGP = computeGlobalParameters(  vSCAR,
                                    fMu,
                                    fSigma,
                                    fDelta,
                                    fProbability1,
                                    fProbability2,
                                    fProbability3 );
    bGP=true;
    if(bGP)
    {
        m_fMuPrimeAA=fMu+fDelta;
        m_fMuPrimeAB=fMu;
        m_fMuPrimeBB=fMu-fDelta;

        m_fSigmaPrimeAA=fSigma;
        m_fSigmaPrimeAB=0.92*fSigma;
        m_fSigmaPrimeBB=fSigma;
        m_fEMSN = (fMu + fDelta) / fSigma;
        Verbose::out(3, "CNAnalysisMethodLOHCyto2::fMu values is: " + AffxString::doubleToString((double)fMu, 6, false));
        Verbose::out(3, "CNAnalysisMethodLOHCyto2::fDelta values is: " + AffxString::doubleToString((double)fDelta, 6, false));
        Verbose::out(3, "CNAnalysisMethodLOHCyto2::fSigma values is: " + AffxString::doubleToString((double)fSigma, 6, false));
    }

    computeGlobalSnpParameters();

    computeHomHetValues(    vSCAR,
                            dHetValueArray,
                            dHomValueArray  );

    return true;
}

void CNAnalysisMethodLOHCyto2::loadInformationVector(std::vector<float> &vInformation)
{
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        vInformation[iIndex] = pobjProbeSet->getInformation();
    }
}

void CNAnalysisMethodLOHCyto2::loadPositionVector(std::vector<int> &vPosition)
{
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        vPosition[iIndex] = pobjProbeSet->getPosition();
    }
}

void CNAnalysisMethodLOHCyto2:: computeHomHetValues(    const std::vector<float>& vSCAR,
                                                        AffxMultiDimensionalArray<double>& dHetValueArray,
                                                        AffxMultiDimensionalArray<double>& dHomValueArray)
{
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
    for(int iValueIndex=0; iValueIndex < m_iNumberOfVIValues; iValueIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        pobjProbeSet->setLoh(numeric_limits<float>::quiet_NaN());

        // We attempt to calculate Hom Het values for ProbeSets with valid SCAR values.
        if( ! pobjProbeSet->processAsSNP())
        {
            //  The default value of using a probeset for hom/het calculations is true.
            //  So if it is not a SNP probeset we must set it to false.
            pobjProbeSet->setValidHomHetExists(false);
            dHomValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
            dHetValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
            continue;
        }
        if(pobjProbeSet->getValidSCARExists())
        {
            if(pobjProbeSet->getAAlleleSignal() == 0.0 || pobjProbeSet->getBAlleleSignal() == 0.0){
                Verbose::warn(1, "CNAnalysisMethodLOHCyto2:: Invalid A or B allele signals were found for the ProbeSet " + pobjProbeSet->getProbeSetName());
                dHomValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
                dHetValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
                //  The default value of using a probeset for hom/het calculations is true.
                //  So if it is no valid SCARs we calculated we must set it to false.
                pobjProbeSet->setValidHomHetExists(false);
                continue;
            }
            if(pobjProbeSet->getMuAA() == CN_INVALID_DOUBLE || pobjProbeSet->getMuAB() == CN_INVALID_DOUBLE || pobjProbeSet->getMuBB() == CN_INVALID_DOUBLE)
            {
                Verbose::warn(1, "CNAnalysisMethodLOHCyto2:: Invalid snp-reference parameters for probe set " + pobjProbeSet->getProbeSetName() + ". LOH computations will be done without this probe set.");
                dHomValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
                dHetValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
                //  The default value of using a probeset for hom/het calculations is true.
                //  So if it is no valid global parameters exist we must set it to false.
                pobjProbeSet->setValidHomHetExists(false);
                continue;
            }

           float  fSCAR = vSCAR[iIndex];
           double dVarianceInflationFactor = m_vVarianceInflationValues[iValueIndex];
           double dHomValueTemp = homDistribution( (double) fSCAR,
                                                     0.5,
                                                     0.5,
                                                     (double) m_fMuPrimeAA,
                                                     dVarianceInflationFactor * (double) m_fSigmaPrimeAA,
                                                     (double) m_fMuPrimeBB,
                                                     dVarianceInflationFactor * (double) m_fSigmaPrimeBB) ;


            double dHetValueTemp = hetDistribution( (double) fSCAR,
                                                    0.375,
                                                    0.375,
                                                    0.25,
                                                    (double) m_fMuPrimeAA,
                                                    dVarianceInflationFactor * (double) m_fSigmaPrimeAA,
                                                    (double) m_fMuPrimeBB,
                                                    dVarianceInflationFactor * (double) m_fSigmaPrimeBB,
                                                    (double) m_fMuPrimeAB,
                                                    dVarianceInflationFactor * (double) m_fSigmaPrimeAB);


           //  The default value of using a probeset for hom/het calculations is true.
           //  Thus any failure of the hom/het computations for any variance inflation will cause the
           //  probeset to be unavailable for all variance inflation computations.
           //  MG May have to modify this functionality if it doesn't match Ben's.

           if( dHomValueTemp > 10.0e-25 && dHetValueTemp > 10.0e-25)
           {
               double dInformation = pobjProbeSet->getInformation();
               dHomValueArray.set(iIndex, iValueIndex, dInformation * log(dHomValueTemp));
               dHetValueArray.set(iIndex, iValueIndex, dInformation * log(dHetValueTemp));
               //MG pobjProbeSet->setValidHomHetExists(true);
           } else
           {
               Verbose::warn(3, "CNAnalysisMethodLOHCyto2:: Invalid hom and or het values have been calculated for ProbeSet " + pobjProbeSet->getProbeSetName()+ " It will not be used in LOH computations.");
               dHomValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
               dHetValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
               pobjProbeSet->setValidHomHetExists(false);
           }
        }
        else
        {
            dHomValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
            dHetValueArray.set(iIndex, iValueIndex, CN_INVALID_DOUBLE);
            pobjProbeSet->setValidHomHetExists(false);
        }
    }
    }
}


float CNAnalysisMethodLOHCyto2::sumVector(    int iStart,
                                              int iEnd,
                                              const std::vector<float>& vInputVector)
{
    float partialSum=0.0;
    for(int i=iStart; i<=iEnd; i++)
    {
        partialSum += vInputVector[i];
    }
    return partialSum;
}

/**
 * @brief Find the runs off LOH
 * @param std::vector<char>& - The SCAR values for each probeset
 * @param std::vector<int>& - The chromosome positions
 * @param std::vector<char>& - The resulting LOH calls (1 = LOH, 0 = !LOH)
 */
bool CNAnalysisMethodLOHCyto2::lohFind( const std::vector<float>& vSCAR,
                                        const std::vector<float>& vInformation,
                                        const std::vector<int>& vPosition,
                                        std::vector<char>& vLoh,
                                        AffxMultiDimensionalArray<double>& dHomValueArray,
                                        AffxMultiDimensionalArray<double>& dHetValueArray,
                                        std::vector<int>& vWindowStart,
                                        std::vector<int>& vWindowLength,
                                        std::vector<double>& vWindowTest,
                                        std::vector<double>& vWindowInflate,
                                        std::vector<char>& vWindowLoh)

{
    int iEndOfSegment= vSCAR.size();
    int iStart=0;
    int iEnd=0;
    float fLoh=0.0;
    float fInformation=0.0;
    int iVarianceInflationIndex=0;


    while(sumVector(iStart, iEnd, vInformation) < m_fMinInformation && iEnd < iEndOfSegment-1) {
        iEnd++;
    }
    // If we don't have don't have enough probesets to generate m_fMinInformation we set LOH = 0;
    if(sumVector(iStart, iEnd, vInformation) < m_fMinInformation)
    {
        for(int i=iStart; i<=iEnd; i++)
        {
            vLoh[i] = 0;
        }
        vWindowStart[iStart] = iStart;
        vWindowLength[iStart] = iEnd - iStart + 1;
        vWindowLoh[iStart] = 0;
        return true;
    }

    iVarianceInflationIndex = determineVarianceInflationIndex(iStart, iEnd, vInformation);

    fLoh = likelihood(iStart, iEnd, dHomValueArray, dHetValueArray, iVarianceInflationIndex);
    if(fLoh > m_fLambdaCritical)
    {
        for(int i=iStart; i<=iEnd; i++)
        {
            vLoh[i] = 1;
        }
    }
    vWindowStart[iStart] = iStart;
    vWindowLength[iStart] = iEnd - iStart + 1;
    vWindowTest[iStart] = fLoh;
    vWindowInflate[iStart] = m_vVarianceInflationValues[iVarianceInflationIndex];
    if(fLoh > m_fLambdaCritical) {
        vWindowLoh[iStart] = 1;
    }
    else {
        vWindowLoh[iStart] = 0;
    }

    while(iEnd < iEndOfSegment-1)
    {
        iStart++;
        iEnd++;

        fInformation = sumVector(iStart, iEnd, vInformation);
        if( fInformation <  m_fMinInformation  )
        {
            while(fInformation < m_fMinInformation && (iEnd < iEndOfSegment-1) )
            {
                iEnd++;
                fInformation = sumVector(iStart, iEnd, vInformation);
            }
        }
        else
        {
            while(sumVector(iStart+1, iEnd, vInformation)  > m_fMinInformation )
            {
                iStart++;
            }
        }

        iVarianceInflationIndex = determineVarianceInflationIndex(iStart, iEnd, vInformation);
        fLoh = likelihood(iStart, iEnd, dHomValueArray, dHetValueArray, iVarianceInflationIndex);

        if(fLoh > m_fLambdaCritical)
        {
            for(int i=iStart; i<=iEnd; i++)
           {
               vLoh[i] = 1;
           }
        }
        vWindowStart[iStart] = iStart;
        vWindowLength[iStart] = iEnd - iStart + 1;
        vWindowTest[iStart] = fLoh;
        vWindowInflate[iStart] = m_vVarianceInflationValues[iVarianceInflationIndex];
        if(fLoh > m_fLambdaCritical) {
           vWindowLoh[iStart] = 1;
        }
        else {
           vWindowLoh[iStart] = 0;
        }
    }

    return true;
}

double CNAnalysisMethodLOHCyto2::homDistribution(double ecks,
                                                double probabilityAA,
                                                double probabilityBB,
                                                double muPrimeAA,
                                                double sigmaPrimeAA,
                                                double muPrimeBB,
                                                double sigmaPrimeBB)
{
    double temp1 = normalDistribution(ecks, muPrimeAA, sigmaPrimeAA);
    double temp2 = normalDistribution(ecks, muPrimeBB, sigmaPrimeBB);
    double temp3 = (probabilityAA * temp1) + (probabilityBB * temp2);
    return temp3;

/*  return ( (probabilityAA * normalDistribution(ecks, muPrimeAA, sigmaPrimeAA)) +
              (probabilityBB * normalDistribution(ecks, muPrimeBB, sigmaPrimeBB))
           );
*/
}


double CNAnalysisMethodLOHCyto2::hetDistribution(double ecks,
                                                double probabilityAA,
                                                double probabilityBB,
                                                double probabilityAB,
                                                double muPrimeAA,
                                                double sigmaPrimeAA,
                                                double muPrimeBB,
                                                double sigmaPrimeBB,
                                                double muPrimeAB,
                                                double sigmaPrimeAB)
{
    double temp1 = normalDistribution(ecks, muPrimeAA, sigmaPrimeAA);
    double temp2 = normalDistribution(ecks, muPrimeBB, sigmaPrimeBB);
    double temp3 = normalDistribution(ecks, muPrimeAB, sigmaPrimeAB);
    double temp4 = (probabilityAA * temp1) + ( probabilityBB * temp2) + (probabilityAB * temp3);
    return temp4;

/*  return ( (probabilityAA * normalDistribution(ecks, muPrimeAA, sigmaPrimeAA)) +
             (probabilityBB * normalDistribution(ecks, muPrimeBB, sigmaPrimeBB)) +
             (probabilityAB * normalDistribution(ecks, muPrimeAB, sigmaPrimeAB))
           );
*/
}

float CNAnalysisMethodLOHCyto2:: likelihood(    int iStart,
                                                int iEnd,
                                                AffxMultiDimensionalArray<double>& dHomValueArray,
                                                AffxMultiDimensionalArray<double>& dHetValueArray,
                                                int iVarianceInflationIndex)
{
    double total=0.0;
    for(int i=iStart; i<=iEnd; i++)
    {
        double temp1 = dHomValueArray.get(i, iVarianceInflationIndex);
        double temp2 = dHetValueArray.get(i, iVarianceInflationIndex);

        if( temp1 != CN_INVALID_DOUBLE && temp2 != CN_INVALID_DOUBLE){
            total+=temp1;
            total-=temp2;
        }
    }
    return (float) total;
}

/**
 * @brief Create new segments
 * @param int - The segment type to create
 * @return int - The number of new segments created
 */
int CNAnalysisMethodLOHCyto2::newSegments(    int iSegmentType,
                                              const std::vector<float>& vSCAR,
                                              const std::vector<double>& vHomValues,
                                              const std::vector<double>& vHetValues)
{
    AffxMultiDimensionalArray<int> vMarkerDistances(getProbeSets()->getCount());

    std::vector<float> vSCARreduced;
    std::vector<double> vHomValuesReduced;
    std::vector<double> vHetValuesReduced;
    for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iIndex);
        if(pobjProbeSet->isUseForCyto2LOH())
        {
            vSCARreduced.push_back(vSCAR[iIndex]);
            vHomValuesReduced.push_back(vHomValues[iIndex]);
            vHetValuesReduced.push_back(vHetValues[iIndex]);
        }
    }

    if (getProbeSets()->getCount() == 0) {return 0;}
    int iSegmentCount = 0;
    int iMarkerCount = 0;
    CNSegment* pobjSegment = NULL;
    CNProbeSet* pobjProbeSet = NULL;
    int iStartIndex = 0;
    if ((iSegmentType == 2) || (iSegmentType == 3) || (iSegmentType == 4) || (iSegmentType == 7))
    {
        for (;(iStartIndex < getProbeSets()->getCount()); iStartIndex++)
        {
            pobjProbeSet = getProbeSets()->getAt(iStartIndex);
            if ((double)pobjProbeSet->getLoh() == (double)pobjProbeSet->getLoh()) // Not = NaN
            {
                break;
            }
        }
    }
    char cCall;
    int iPrevMarkerPosition;
    int iPosition;
    float fConfidence=0;
    int iConfidenceStart=0;
    int iConfidenceEnd=-1;
    int iPrevCNState = -1;
    if (iStartIndex < getProbeSets()->getCount() )
    {
        pobjProbeSet = getProbeSets()->getAt(iStartIndex);
        while( (! pobjProbeSet->isUseForCyto2LOH())  && iStartIndex < getProbeSets()->getCount() )
        {
            iStartIndex++;
            if(iStartIndex==getProbeSets()->getCount()){
                Verbose::warn(1, "CNAnalysisMethodLOHCyto2::No probesets with valid information and/or Hom Het values.  Possibly a missing snp-reference file or information section in snp-reference file.");
                return 0;
            }
            pobjProbeSet = getProbeSets()->getAt(iStartIndex);
        }
        iPosition = pobjProbeSet->getPosition();
        iPrevMarkerPosition = iPosition;
        iPrevCNState = pobjProbeSet->getCNState();
        cCall = pobjProbeSet->getSegmentTypeCall(iSegmentType, getExperiment()->getCNCallGenderAsInt(), m_iXChromosome, m_iYChromosome);

        pobjSegment = new CNSegment;
        pobjSegment->setSegmentType(iSegmentType);
        pobjSegment->setChromosome(pobjProbeSet->getChromosome());
        pobjSegment->setStartPosition(pobjProbeSet->getPosition());
        pobjSegment->setEndPosition(pobjProbeSet->getPosition());
        pobjSegment->setCall(cCall);
        iMarkerCount = 0;

        for (int iIndex = iStartIndex; (iIndex < getProbeSets()->getCount()); iIndex++)
        {
            pobjProbeSet = getProbeSets()->getAt(iIndex);
            if(pobjProbeSet->isUseForCyto2LOH()) {
                iPosition = pobjProbeSet->getPosition();
                cCall = pobjProbeSet->getSegmentTypeCall(       iSegmentType,
                                                                getExperiment()->getCNCallGenderAsInt(),
                                                                m_iXChromosome,
                                                                m_iYChromosome);

                // For CNNeutralLOH segments, keep cCall = 1 if current probe set's CN = -1
                // and last seen probe set with a valid CN had CN = 2.
                if (iSegmentType == 4 && cCall == 0  &&
                    pobjProbeSet->getCNState() == -1 && pobjProbeSet->getLoh() == 1 && iPrevCNState == 2)
                {
                    cCall = (char)1;
                }
                // Same for normal diploid segments
                if (iSegmentType == 5 && cCall == 0  &&
                    pobjProbeSet->getCNState() == -1 && pobjProbeSet->getLoh() == 0 && iPrevCNState == 2)
                {
                    cCall = (char)1;
                }
                if ((iSegmentType == 2) || (iSegmentType == 3) || (iSegmentType == 4) || (iSegmentType == 7))
                {
                    if ((double)pobjProbeSet->getLoh() != (double)pobjProbeSet->getLoh()) // NaN
                    {
                        continue;
                    }
                }
                if ((pobjProbeSet->getChromosome() != pobjSegment->getChromosome()) || (cCall != pobjSegment->getCall())
                       || (iPosition - iPrevMarkerPosition) > m_iLohCNSegSeparation )
                {
                    pobjSegment->setMeanMarkerDistance(vMarkerDistances.mean(iMarkerCount - 1));
                    pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
                    pobjSegment->setMarkerCount(iMarkerCount);
                    if(iMarkerCount==1)
                    {
                        fConfidence = 0.0;
                        pobjSegment->setCall( (char) 0 );
                    } else
                    {
                        fConfidence = computeConfidence(iConfidenceStart, iConfidenceEnd, vSCARreduced, vHomValuesReduced, vHetValuesReduced);
                    }
                    iConfidenceStart = iConfidenceEnd+1;
                    pobjSegment->setConfidence(fConfidence);
                    m_vSegments.add(pobjSegment);
                    iSegmentCount++;
                    pobjSegment = new CNSegment;
                    pobjSegment->setSegmentType(iSegmentType);
                    pobjSegment->setCall(cCall);
                    pobjSegment->setChromosome(pobjProbeSet->getChromosome());
                    pobjSegment->setStartPosition(pobjProbeSet->getPosition());
                    iMarkerCount = 0;
                }
                iConfidenceEnd++;
                iMarkerCount++;
                pobjSegment->setEndPosition(pobjProbeSet->getPosition());
                if (iMarkerCount > 1)
                {
                    vMarkerDistances.set((iMarkerCount - 2), (pobjProbeSet->getPosition() - iPrevMarkerPosition));
                }
                iPrevMarkerPosition = iPosition;
                if (pobjProbeSet->getCNState() != -1) {
                    iPrevCNState = pobjProbeSet->getCNState();
                }
            }
        }
    }
    if (pobjSegment != NULL)
    {
        pobjSegment->setMarkerCount(iMarkerCount);
        pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
        fConfidence = computeConfidence(iConfidenceStart, iConfidenceEnd, vSCARreduced, vHomValuesReduced, vHetValuesReduced);
        pobjSegment->setConfidence(fConfidence);
        pobjSegment->setMeanMarkerDistance(vMarkerDistances.mean(iMarkerCount - 1));
        m_vSegments.add(pobjSegment); iSegmentCount++;
}
    return iSegmentCount;
}


float CNAnalysisMethodLOHCyto2::computeConfidence(      int iStart,
                                                        int iEnd,
                                                        const std::vector<float>& vSCAR,
                                                        const std::vector<double>& vHomValues,
                                                        const std::vector<double>& vHetValues)

{
    double homTotal=0.0;
    double hetTotal=0.0;
    double temp1=0.0;
    double temp2=0.0;
    for(int i=iStart; i<=iEnd; i++)
    {
        temp1 = vHomValues[i];
        temp2 = vHetValues[i];

        if( temp1 != CN_INVALID_DOUBLE && temp2 != CN_INVALID_DOUBLE)
        {
            homTotal+=temp1;
            hetTotal+=temp2;
        }
    }

    double fExpHom = exp(homTotal);
    double fExpHet = exp(hetTotal);

    double fConfidence = fExpHom / (fExpHom + fExpHet);
    if(fConfidence == fConfidence) // !isnan
    {
        return (float) fConfidence;
    } else
    {
        if(homTotal > hetTotal)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
     }
}

int CNAnalysisMethodLOHCyto2::determineVarianceInflationIndex(     int iStart,
                                                                int iEnd,
                                                                const std::vector<float>& vInformation)
{
    if(m_iNumberOfVIValues == 1)
    {
        return 0;
    }

    std::vector<float> vSCARStoSort;
    for(int iIndex=iStart; iIndex<=iEnd; iIndex++)
    {
        vSCARStoSort.push_back(vInformation[iIndex]);
    }
    sort(vSCARStoSort.begin(), vSCARStoSort.end());
    int iQuantileIndex = Util::round( (iEnd-iStart+1)*0.25) - 1;
    //int iQuantileIndex = (int)(nearbyint(round( (iEnd-iStart+1)*0.25) ) - 1);

    float fVarianceInflationInput = vSCARStoSort[iQuantileIndex];
    int iPositionIndex=0;
    for(; iPositionIndex < m_vInformationRanges.size(); iPositionIndex++)
    {
        if(fVarianceInflationInput <= m_vInformationRanges[iPositionIndex])
        {
            break;
        }
    }
    if(fVarianceInflationInput > m_vInformationRanges[m_vInformationRanges.size()-1])
    {
        return m_vInformationRanges.size();
    }
    else
    {
        return iPositionIndex;
    }
}


void CNAnalysisMethodLOHCyto2::writeLOH(    std::vector<float>& vSCAR,
                                            std::vector<int>& vGenomeChunkStart,
                                            std::vector<int>& vGenomeChunkLength,
                                            std::vector<int>& vWindowStart,
                                            std::vector<int>& vWindowLength,
                                            std::vector<double>& vWindowTest,
                                            std::vector<double>& vWindowInflate,
                                            std::vector<char>& vWindowLoh )
{
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv_global = NULL;
    affx::File5_Tsv* tsv_model = NULL;
    affx::File5_Tsv* tsv_scar= NULL;
    affx::File5_Tsv* tsv_chunk = NULL;
    affx::File5_Tsv* tsv_window = NULL;
    affx::File5_Tsv* tsv_segment = NULL;
    affx::File5_Tsv* tsv_summary = NULL;

    std::string analysisString = "analysis";
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),analysisString,
                                  strExperimentName+".loh.a5");

    Verbose::out(3, "Writing LOH data for " + fileName);

    try {
        file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
        group5 = file5.openGroup("LOH", affx::FILE5_REPLACE);

        tsv_global = group5->openTsv("GlobalParams", affx::FILE5_REPLACE);
        tsv_global->defineColumn(0,0,"t_min", affx::FILE5_DTYPE_FLOAT);
        tsv_global->defineColumn(0,1,"lambda_crit", affx::FILE5_DTYPE_FLOAT);

        tsv_global->set_f(0,0, (float) m_fMinInformation );
        tsv_global->set_f(0,1, (float) m_fLambdaCritical );
        tsv_global->writeLevel(0);
        tsv_global->close();
        delete tsv_global;

        tsv_model = group5->openTsv("ModelParams", affx::FILE5_REPLACE);
        tsv_model->defineColumn(0,0,"muAA", affx::FILE5_DTYPE_FLOAT);
        tsv_model->defineColumn(0,1,"muAB", affx::FILE5_DTYPE_FLOAT);
        tsv_model->defineColumn(0,2,"muBB", affx::FILE5_DTYPE_FLOAT);
        tsv_model->defineColumn(0,3,"sigmaAA", affx::FILE5_DTYPE_FLOAT);
        tsv_model->defineColumn(0,4,"sigmaAB", affx::FILE5_DTYPE_FLOAT);
        tsv_model->defineColumn(0,5,"sigmaBB", affx::FILE5_DTYPE_FLOAT);

        tsv_model->set_f(0,0, (float) m_fMuPrimeAA );
        tsv_model->set_f(0,1, (float) m_fMuPrimeAB );
        tsv_model->set_f(0,2, (float) m_fMuPrimeBB );
        tsv_model->set_f(0,3, (float) m_fSigmaPrimeAA );
        tsv_model->set_f(0,4, (float) m_fSigmaPrimeAB );
        tsv_model->set_f(0,5, (float) m_fSigmaPrimeBB );
        tsv_model->writeLevel(0);
        tsv_model->close();
        delete tsv_model;

        tsv_scar = group5->openTsv("SCAR", affx::FILE5_REPLACE);
        tsv_scar->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING, 20);
        tsv_scar->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_INT);
        tsv_scar->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
        tsv_scar->defineColumn(0,3,"SCAR", affx::FILE5_DTYPE_FLOAT);

        for(int iProbeSetIndex=0; iProbeSetIndex < m_pvProbeSets->getCount(); iProbeSetIndex++) {
            CNProbeSet* pProbeSet = m_pvProbeSets->getAt(iProbeSetIndex);
            if ( pProbeSet->isUseForCyto2LOH() ) {
                tsv_scar->set_string(0,0, pProbeSet->getProbeSetName());
                tsv_scar->set_i(0,1, pProbeSet->getChromosome());
                tsv_scar->set_i(0,2, pProbeSet->getPosition());
                tsv_scar->set_f(0,3, (float) vSCAR[iProbeSetIndex] );
                tsv_scar->writeLevel(0);
            }
        }
        tsv_scar->close();
        delete tsv_scar;

        tsv_chunk = group5->openTsv("Chunks", affx::FILE5_REPLACE);
        tsv_chunk->defineColumn(0,0,"StartChromosome", affx::FILE5_DTYPE_INT);
        tsv_chunk->defineColumn(0,1,"StartPosition", affx::FILE5_DTYPE_INT);
        tsv_chunk->defineColumn(0,2,"Length", affx::FILE5_DTYPE_INT);

        for(unsigned int iChunkIndex=0; iChunkIndex < vGenomeChunkStart.size(); iChunkIndex++) {
            CNProbeSet* pProbeSet = m_pvProbeSets->getAt(vGenomeChunkStart[iChunkIndex]);
            tsv_chunk->set_i(0,0, pProbeSet->getChromosome());
            tsv_chunk->set_i(0,1, pProbeSet->getPosition());
            tsv_chunk->set_i(0,2, vGenomeChunkLength[iChunkIndex] );
            tsv_chunk->writeLevel(0);
        }
        tsv_chunk->close();
        delete tsv_chunk;

        tsv_window = group5->openTsv("Windows", affx::FILE5_REPLACE);
        tsv_window->defineColumn(0,0,"StartChromosome", affx::FILE5_DTYPE_INT);
        tsv_window->defineColumn(0,1,"StartPosition", affx::FILE5_DTYPE_INT);
        tsv_window->defineColumn(0,2,"Length", affx::FILE5_DTYPE_INT);
        tsv_window->defineColumn(0,3,"Test", affx::FILE5_DTYPE_FLOAT);
        tsv_window->defineColumn(0,4,"Inflate", affx::FILE5_DTYPE_FLOAT);
        tsv_window->defineColumn(0,5,"Loh", affx::FILE5_DTYPE_CHAR);

        for(unsigned int iWindowIndex=0; iWindowIndex < vWindowStart.size(); iWindowIndex++) {
            CNProbeSet* pProbeSet = m_pvProbeSets->getAt(vWindowStart[iWindowIndex]);
            tsv_window->set_i(0,0, pProbeSet->getChromosome());
            tsv_window->set_i(0,1, pProbeSet->getPosition());
            tsv_window->set_i(0,2, vWindowLength[iWindowIndex] );
            tsv_window->set_f(0,3, (float) vWindowTest[iWindowIndex] );
            tsv_window->set_f(0,4, (float) vWindowInflate[iWindowIndex] );
            tsv_window->set_c(0,5, vWindowLoh[iWindowIndex] );
            tsv_window->writeLevel(0);
        }
        tsv_window->close();
        delete tsv_window;

        tsv_segment = group5->openTsv("Segments", affx::FILE5_REPLACE);
        tsv_segment->defineColumn(0,0,"StartChromosome", affx::FILE5_DTYPE_INT);
        tsv_segment->defineColumn(0,1,"StartPosition", affx::FILE5_DTYPE_INT);
        tsv_segment->defineColumn(0,2,"EndPosition", affx::FILE5_DTYPE_INT);
        tsv_segment->defineColumn(0,3,"NumberMarkers", affx::FILE5_DTYPE_INT);
        tsv_segment->defineColumn(0,4,"MeanMarkerDistance", affx::FILE5_DTYPE_INT);
        tsv_segment->defineColumn(0,5,"LOHState", affx::FILE5_DTYPE_CHAR);
        tsv_segment->defineColumn(0,6,"Confidence", affx::FILE5_DTYPE_FLOAT);

        for(unsigned int iSegmentIndex=0; iSegmentIndex < m_vSegments.size(); iSegmentIndex++) {
            CNSegment *pSegment = m_vSegments.getAt(iSegmentIndex);
            if( pSegment->getSegmentType() == 3 ) {
                tsv_segment->set_i(0,0, pSegment->getChromosome());
                tsv_segment->set_i(0,1, pSegment->getStartPosition() );
                tsv_segment->set_i(0,2, pSegment->getEndPosition() );
                tsv_segment->set_i(0,3, pSegment->getMarkerCount() );
                tsv_segment->set_i(0,4, (int) pSegment->getMeanMarkerDistance() );
                tsv_segment->set_c(0,5, pSegment->getCall() );
                tsv_segment->set_f(0,6, (float) pSegment->getConfidence() );
                tsv_segment->writeLevel(0);
            }
        }
        tsv_segment->close();
        delete tsv_segment;

        std::pair<char, float> data;
        tsv_summary = group5->openTsv("Summary", affx::FILE5_REPLACE);
        tsv_summary->defineColumn(0,0,"Metric", affx::FILE5_DTYPE_STRING, 10);
        tsv_summary->defineColumn(0,1,"Value", affx::FILE5_DTYPE_FLOAT);
        for(unsigned int iChrIndex=0; iChrIndex < getExperiment()->getNumberOfChromosomesToReport(); iChrIndex++) {
            data = getExperiment()->getChromosomeLOH(iChrIndex);
            if( data.first != 24 ) {
                tsv_summary->set_string(0,0, "Chr"+ToStr((int)data.first) );
            }
            else {
                tsv_summary->set_string(0,0, "ChrX" );
            }
            tsv_summary->set_f(0,1, (float) data.second );
            tsv_summary->writeLevel(0);
        }

        tsv_summary->set_string(0,0, "Overall" );
        tsv_summary->set_f(0,1, (float) getExperiment()->getGenomeLOH() );
        tsv_summary->writeLevel(0);

        tsv_summary->close();
        delete tsv_summary;

        delete group5;
        file5.close();
    }
    catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump intensities to file."));}
}
