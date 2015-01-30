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

#ifndef _CNReferenceMethodWaveCorrection_h_
#define _CNReferenceMethodWaveCorrection_h_

/**
 * @file CNWaveCorrection.h
 *
 * @brief This header contains the CNWaveCorrection class definition.
 */

#define NEWMAT_USING_FLOAT
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNCychp.h"
//
#include "util/AffxArray.h"
//
#include "../external/newmat/newmatap.h"
//
#include <vector>
//

class CNWaveCorrectionProbeSet
{
public:
	int Order;
    AffxString ProbeSetName;
    unsigned char Chromosome;
    unsigned int Position;
    float Signal;
    float MedianSignal;
    float XXMedianSignal;
    float YMedianSignal;
    double Log2Ratio;
    std::vector<double> Waves;
    std::vector<float> AllCovariates;
    unsigned char ProcessFlag;

    float getCovariateValue(int i);
    unsigned char getProcessFlag();
    unsigned char getChromosome();
    double getLog2Ratio();
    void setLog2Ratio(double f);

    CNWaveCorrectionProbeSet()
    {
		Order = 0;
        Chromosome = 0;
        Position = 0;
        Signal = 0;
        MedianSignal = 0;
        XXMedianSignal = 0;
        YMedianSignal = 0;
        Log2Ratio = 0;
        ProcessFlag = (unsigned char)0;
    }

    int compareTo(CNWaveCorrectionProbeSet& that, int iCompareCode)
    {
        int iCompareResult = 0;
        switch (iCompareCode)
        {
        case 0:
            iCompareResult = strcmp(ProbeSetName.c_str(), that.ProbeSetName.c_str());
            break;
        case 1:
            iCompareResult = AffxArray<unsigned char>::compare(Chromosome, that.Chromosome);
            if (iCompareResult == 0) {iCompareResult = AffxArray<unsigned int>::compare(Position, that.Position);}
            if (iCompareResult == 0) {iCompareResult = strcmp(ProbeSetName.c_str(), that.ProbeSetName.c_str());}
            break;
        }
        return iCompareResult;
    }

    template<int k> struct ComparePred {
        bool operator()(const CNWaveCorrectionProbeSet* lhs, const CNWaveCorrectionProbeSet* rhs) const {
            Err::errAbort("CNWaveCorrectionProbeSet: ComparePred instantiated with an invalid compare code = " + ToStr(k));
            return false;
        }
    };
};

template<> struct CNWaveCorrectionProbeSet::ComparePred<0> {
    bool operator()(const CNWaveCorrectionProbeSet* lhs, const CNWaveCorrectionProbeSet* rhs) const {
        return strcmp(lhs->ProbeSetName.c_str(), rhs->ProbeSetName.c_str()) < 0;
    }
};
template<> struct CNWaveCorrectionProbeSet::ComparePred<1> {
    bool operator()(const CNWaveCorrectionProbeSet* lhs, const CNWaveCorrectionProbeSet* rhs) const {
        int iCompareResult = AffxArray<unsigned char>::compare(lhs->Chromosome, rhs->Chromosome);
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<unsigned int>::compare(lhs->Position, rhs->Position);
        }
        if (iCompareResult == 0) {
            iCompareResult = strcmp(lhs->ProbeSetName.c_str(), rhs->ProbeSetName.c_str());
        }
        return iCompareResult < 0;
    }
};

/**
 * @brief  A class for calulating WaveCorrection.
 *
 */
class CNReferenceMethodWaveCorrection : public CNAnalysisMethod
{
protected:
    AffxString m_strTempFileName;
    AffxString m_strReferenceFileName;
    AffxString m_strGroup;
    AffxString m_strTsv;
    AffxString m_strSelectedQC;

    int m_iLastAutosomeChromosome;
    int m_iAutosomeCount;
    int m_iXChromosomeCount;
    int m_iYChromosomeCount;
    double m_dTrim;
    double m_dPercentile;
    int m_iWaveCount;
    bool m_bDemean;
    double m_dWaveCutoff;
    bool m_bKeepTempData;
    bool m_bForce;
    int m_iWavinessSegCountCutoff;
    bool m_bUseHighWavinessSegCount;
    std::vector<double> m_vAutosomeCV;
    std::vector<bool> m_vHasXX;
    std::vector<bool> m_vHasY;
    std::vector<double> m_vYTargets;
    std::vector<bool> m_vPassedQC;
    int m_iRefWaveCount;
    double m_dCNQCCutoff;
    double m_dSNPQCCutoff;
    unsigned int m_uiPassedQCCount;
    bool m_bCyto2;

public:
    CNReferenceMethodWaveCorrection();
    virtual ~CNReferenceMethodWaveCorrection();

    static std::string getType() { return "wave-correction-reference-method"; }
    static std::string getDescription() { return "CopyNumber WaveCorrection"; }
    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string> & param);

    virtual AffxString getName() { return getType(); }
    virtual void run();
    virtual bool isSegmentTypeAnalysis() { return false; }

private:
    void calculateXY(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iProbeSetCount, bool& bXX, bool& bY);
    void adjustLog2Ratios(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iProbeSetCount, int iType, unsigned int uiCelIndex);
    void loadPCALikeVector(const std::string& strTempFileName, AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iWaveIndex);
    void calculateChromosomeStatistics(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets);
    void calibrateMedianAutosomeMedian(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets);
    void loadMatrix(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, std::vector<int>& vProbeSetsIndexes, unsigned int uiCelIndex, Matrix& mx);
    double corr(Matrix& mx1, Matrix& mx2);
//    double norm(Matrix& v);
    void updateReference(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iProbeSetCount);
    void covariateAdjustLog2Ratio(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, unsigned int uiCelIndex);
    void writeWaveCorrectionLog2Ratios(std::string l2rFileInfix, AffxArray<CNWaveCorrectionProbeSet>& vProbeSets);
};

#endif // _CNReferenceModuleWaveCorrection_h_


