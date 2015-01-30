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

#ifndef _CNLog2RatioData_H_
#define _CNLog2RatioData_H_
/**
 * @file CNLog2RatioData.h
 *
 * @brief This header contains the ProbeSet and Log2RatioData class definitions.
 */

#include "copynumber/CNExperiment.h"
#include "copynumber/CNProbeSet.h"
#include "copynumber/CNSegment.h"
//
#include "util/AffxArray.h"
#include "util/AffxSplitArray.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

/**
 * @brief  A class to the log2 ratio data as used by the log2 ratio engine.
 */
class CNLog2RatioData
{
private:
    bool m_bLog2RatioEngine;
    BaseEngine* m_pEngine;
    CNExperimentArray m_arExperiments;
    AffxSplitArray<CNProbeSet> m_vProbeSets;
    CNProbeSetArray m_arProbeSets;
    CNProbeSetArray m_arProbeSetsAlternateSort;

    AffxMultiDimensionalArray<float> m_mxAAlleleEstimates;
    AffxMultiDimensionalArray<float> m_mxBAlleleEstimates;
    AffxMultiDimensionalArray<char> m_mxGenotypeCalls;
    AffxMultiDimensionalArray<float> m_mxGenotypeConfidences;

    AffxMultiDimensionalArray<float> m_v1;

        double m_dMuPrimeAA;
        double m_dMuPrimeAB;
        double m_dMuPrimeBB;

        double m_dSigmaPrimeAA;
        double m_dSigmaPrimeAB;
        double m_dSigmaPrimeBB;

public:
    CNLog2RatioData();
    virtual ~CNLog2RatioData();

    void clear();
    void clearProbeSets();

    void setLog2RatioEngine(bool b) {m_bLog2RatioEngine = b;}
    void setEngine(BaseEngine* pEngine) {m_pEngine = pEngine;}
    CNExperimentArray* getExperiments() {return &m_arExperiments;}
    void setExperiments(CNExperimentArray& v) {m_arExperiments = v;}
    CNProbeSetArray* getProbeSets() {return &m_arProbeSets;}
    CNProbeSetArray* getProbeSetsAlternateSort() {return &m_arProbeSetsAlternateSort;}
    AffxMultiDimensionalArray<float>* getAAlleleEstimates() {return &m_mxAAlleleEstimates;}
    AffxMultiDimensionalArray<float>* getBAlleleEstimates() {return &m_mxBAlleleEstimates;}
    AffxMultiDimensionalArray<char>* getGenotypeCalls() {return &m_mxGenotypeCalls;}
    AffxMultiDimensionalArray<float>* getGenotypeConfidences() {return &m_mxGenotypeConfidences;}

    AffxMultiDimensionalArray<float>& getV1() {initialize(m_v1); return m_v1;}

    void loadAnnotation(bool bNewProbeSets);

    bool loadModelFile(const AffxString& strFileName);
    bool writeModelFile(const AffxString& strFileName);
    void loadCyto2ModelFile(const AffxString& strFileName);
    void writeCyto2ModelFile(const AffxString& strFileName);

    bool loadSnpReferenceFile(const AffxString& strFileName);

    void loadGenotypeReport(const AffxString& strFileName);
    void dumpGenotypeReport(const AffxString& strFileName);
    void dumpCyto2Report(const AffxString& strFileName);

    int loadExperiments(const AffxString& strFileName);
    void loadQCMetrics(const AffxString& strFileName);

    bool loadAlleleEstimatesByProbeSet(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess);
    bool loadAlleleEstimatesByProbeSetHdf5(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess);
    bool loadAlleleEstimatesByExperiment(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess);
    bool loadAlleleEstimatesByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess);

    bool loadGenotypeCallsByProbeSet(const AffxString& strCallsFileName, int iProbeSetIndex, int iProbeSetsToProcess);
    bool loadGenotypeCallsByProbeSetHdf5(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess);
    bool loadGenotypeCallsByExperiment(const AffxString& strCallsFileName, int iExperimentIndex, int iExperimentsToProcess);
    bool loadGenotypeCallsByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess);

    bool loadGenotypeConfidencesByExperiment(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess);
    bool loadGenotypeConfidencesByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess);

        bool loadGenderOverrideFile(const AffxString& strFileName);

    void calculateMedianSignal(int iProbeSetIndex, int iProbeSetsToProcess);
    void calculateMedianSignalsPerGenotype(int iProbeSetIndex, int iProbeSetsToProcess);
    void calculateXYMedianSignal(int iProbeSetIndex, int iProbeSetsToProcess);

    void loadProbeSets();


private:
    /**
     * @brief Initialize the vector
     * @param AffxMultiDimensionalArray<float>& - The vector to initialize
     */
    void initialize(AffxMultiDimensionalArray<float>& v)
    {
        if (v.getXDimension() == m_arProbeSets.getCount()) {v.initialize();}
        else {v.initialize(m_arProbeSets.getCount());}
    }

    void loadProbeSetNamesFromRestrictList(AffxArray<AffxString>& ar);

    bool isNaN(const double& d) {return d != d;}

    double getSnpSignalEstimate(int iCol, int iRow, CNProbeSet* p);
    double getAllelicDifference(int iCol, int iExperiment, int iRow);

    void overrideGcContent();

    void readCovariatesFile();
};


#endif


