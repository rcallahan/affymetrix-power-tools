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

#ifndef _CNAnalysisMethodChipstream_H_
#define _CNAnalysisMethodChipstream_H_
/**
 * @file CNAnalysisMethodChipstream.h
 *
 * @brief This header contains the CNAnalysisMethodChipstream class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNProbe.h"
//
#include "chipstream/AdapterTypeNormTran.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "chipstream/QuantLabelZ.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "plier/affyplier.h"
//

/**
 * @brief  The Chipstream Analysis Method.
 *
 */
class CNAnalysisMethodChipstream : public CNAnalysisMethod
{
public:
    static std::string getType() {return "chipstream";}
    static std::string getDescription() {return "Copynumber Chipstream";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodChipstream();
    virtual ~CNAnalysisMethodChipstream();

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run();

    void calculateXY();

    double getXXRatio() {return m_dXXRatio;}
    double getYRatio() {return m_dYRatio;}

protected:
    virtual void determineSketchProbeSets();
        void loadSnpProbeSets(std::vector<CNProbeSet*> &vSnpProbeSets);
    void normalizeIntensities();
    void extractSketch(std::vector<float>& vData, std::vector<float>& vSketch);
    float interpolate(float fIntensity, AffxMultiDimensionalArray<float>& vReferenceSketch, AffxMultiDimensionalArray<float>& vChipSketch, double dRatio);

    void adjustIntensitiesPDNN();
    void adjustIntensitiesHighPassFilter();
    void adjustIntensitiesUsingCovariateSignalAdjustment();

    void calculateSignalEstimates();
    void calculateGenotypes();

    virtual void newProbes();
    void deleteProbes();

    void newPlier(int iChipCount);
    void deletePlier();

    void newSketch();
    void deleteSketch();

    void loadIntensities();
    void adapterTypeNormalization(affymetrix_fusion_io::FusionCELData& cel);
    void setProbeSetSignal(CNProbe* p, double dSignal, double dMedianIntensity);

    void initializeSnpParameters(snp_param& sp, int iCM);
    void setBrlmmpParameters(snp_param& sp);
    void loadSnpDistributions();

    bool loadGenotypeCallsFromFile();
    AffxString getSampleId(const AffxString& strCelFileName);

    void calculateRawSNPQCs(ChipLayout* pLayout);

    void adjustAlternateIntensitiesUsingCovariateSignalAdjustment(affymetrix_fusion_io::FusionCELData& cel, 
                                                                  std::vector<float>& adjustedIntensities);

protected:
    int m_iPrevGenderCode;
    int m_iGenderCode;

    unsigned int m_uiMaxProbeCount;
    double* m_pdProbeEffects;
    double* m_pdChipEffects;
    double** m_ppdPM;
    double** m_ppdMM;
    double** m_ppdResiduals;
    caffyplier* m_pPlier;
    int m_iPlierChipCount;

    double m_dXXRatio;
    double m_dYRatio;

    AdapterTypeNormTran* m_pobjAdapterTypeNormTran;
    CNProbe* m_pProbes;

    int m_iInstanceCount;
    CnProbeGenderCelListener* m_pobjGenderCaller;
    std::vector<float>* m_pvSNPReferenceSketch;
    std::vector<float>* m_pvCNReferenceSketch;
    std::vector<float>* m_pvChipSketch;
    std::vector<float>* m_pvSNPData;
    std::vector<float>* m_pvCNData;

        double m_fNDBandwidth;


};

#endif


