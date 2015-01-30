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

#ifndef _CNAnalysisMethodReference_H_
#define _CNAnalysisMethodReference_H_
/**
 * @file CNAnalysisMethodReference.h
 *
 * @brief This header contains the CNAnalysisMethodReference class definition.
 */

#include "copynumber/CNAnalysisMethodChipstream.h"
//
#include "util/AffxBinaryFile.h"
#include "util/AffxMultiDimensionalArray.h"
//

// @todo this needs to be renamed to "CN_INVALID_DOUBLE and put in the correct
//       header file. -- jhg
#ifndef INVALID_DOUBLE
#define INVALID_DOUBLE  (-9999.9999)
#endif

/**
 * @brief  The Reference Analysis Method.
 *
 */
class CNAnalysisMethodReference : public CNAnalysisMethodChipstream
{
public:
    static std::string getType() {return "reference";}
    static std::string getDescription() {return "Copynumber Reference";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodReference();
    virtual ~CNAnalysisMethodReference() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run() {}

    void runPart1();
    void runPart2();

    AffxMultiDimensionalArray<char,long int>& getGenotypeCalls() {return m_mxGenotypeCalls;}

    AffxString getTempFileName() {return m_strTempFileName;}
    void setTempFileName(const AffxString& str) {m_strTempFileName = str;}

    /**
    * Utility function for reading a ChipLayout based on the options passed
    * to the base engine (i.e. cdf-file or spf-file)
     */
    static void loadLayout(ChipLayout *pLayout, BaseEngine *engine);

protected:
    int m_iWarningCount;
    AffxString m_strTempFileName;
    AffxMultiDimensionalArray<char,long int> m_mxGenotypeCalls;
    double m_dAASum;
    double m_dABSum;
    double m_dBBSum;
    double m_dAASumSquares;
    double m_dABSumSquares;
    double m_dBBSumSquares;
    int m_iAAIndex;
    int m_iABIndex;
    int m_iBBIndex;
    float m_fLambdaSpacing;
    float m_fEMSNSpacing;
    float m_fLambdaStartingValue;
    float m_fEMSNStartingValue;
    int   m_iNumberOfLambdaPoints;
    int   m_iNumberOfEMSNPoints;
    AffxMultiDimensionalArray<double> m_dEMSNArray;

    std::vector<float> m_vInformationRanges;
    std::vector<float> m_vVarianceInflationValues;

    std::vector<float> m_vThreePeakFLD_X;
    std::vector<float> m_vThreePeakFLD_Y;
    std::vector<float> m_vFourPeakFLD_X;
    std::vector<float> m_vFourPeakFLD_Y;

    std::vector<float> m_vThreePeakShrink_X;
    std::vector<float> m_vThreePeakShrink_Y;
    std::vector<float> m_vFourPeakShrink_X;
    std::vector<float> m_vFourPeakShrink_Y;


protected:
    void calculateReferenceSketch();
    bool writeReferenceSketch(const std::string& strName, std::vector<double>& vReferenceSketch);
    virtual void newProbes();
    void processIntensities();
    void callPlier(int iProbeSetIndex, char cAllele, std::vector<int>& vProbeIDs, int iProbeCount, affx::File5_Tsv* tsv5, std::vector<float>& vMedianIntensities);
    void setSignals(    char cAllele,
                                double* pdAAlleleSignal,
                                double* pdBAlleleSignals,
                                affx::File5_Tsv** pptsv5Signals,
                                AffxString strProbeSetName);
    void callBrlmmp(int iProbeSetIndex, int iProbeCount, double* pdAAlleleSignals, double* pdBAlleleSignals, std::vector<affx::GType>& vCalls, affx::File5_Tsv* tsv5);
    void calculateMedians(int iProbeSetIndex, double* pdAAlleleSignals, double* pdBAlleleSignals, std::vector<affx::GType>& vCalls);
    void calculateXYMedians(int iProbeSetIndex, double* pdAAlleleSignals, double* pdBAlleleSignals);
    void writeCopyNumberReference();
    void writeAdapterTypeMatrix(AffxMultiDimensionalArray<char>& mx);
        void writeSnpReference();

        bool loadSnpReferenceFile(const AffxString& strFileName);
        void calculateSnpReference(    int iProbeSetIndex,
                                        double* pdAAlleleSignals,
                                        double* pdBAlleleSignals,
                                        std::vector<affx::GType>& vCalls);

        float computeSCAR(CNProbeSet* pobjProbeSet);

        void warn(int iVerbosity, const std::string& strMessage);


};

#endif


