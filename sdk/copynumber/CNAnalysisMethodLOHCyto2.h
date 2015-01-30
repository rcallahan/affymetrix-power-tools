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

#ifndef _CNAnalysisMethodLOHCyto2_H_
#define _CNAnalysisMethodLOHCyto2_H_
/**
 * @file CNAnalysisMethodLOHCyto2.h
 *
 * @brief This header contains the CNAnalysisMethodLOHCyto2 class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

#define CN_ANALYSIS_METHOD_LOH_CYTO2_VERSION "1.0"

/**
 * @brief  The Cyto2LOH analysis method.
 *
 */
class CNAnalysisMethodLOHCyto2 : public CNAnalysisMethod
{
private:
    int m_iLohCNSegSeparation;

        float m_fMuPrimeAA;
        float m_fMuPrimeAB;
        float m_fMuPrimeBB;

        float m_fSigmaPrimeAA;
        float m_fSigmaPrimeAB;
        float m_fSigmaPrimeBB;

        float m_fMinInformation;
        float m_fMinInformationCommandLineValue;

        float m_fLambdaCritical;
        float m_fLambdaCriticalCommandLineValue;
        float m_fLambdaSpacing;
        float m_fLambdaStartingValue;
        int m_iNumberOfLambdaPoints;

        float m_fEMSN;
        float m_fEMSNSpacing;
        float m_fEMSNStartingValue;
        int m_iNumberOfEMSNPoints;

        AffxMultiDimensionalArray<double> m_dEMSNArray;


        // MG add these to constructor, copy
        std::vector<float> m_vInformationRanges;
        std::vector<float> m_vVarianceInflationValues;
        int m_iNumberOfVIValues;

        bool m_bSNPAlreadyLoaded;

public:
        // constructor
    CNAnalysisMethodLOHCyto2();

    virtual ~CNAnalysisMethodLOHCyto2();

    static std::string getType();
    static std::string getDescription();
    static std::string getVersion();

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual AffxString getName();
    virtual void run();
    virtual bool isSegmentTypeAnalysis();

protected:
        void reset();

        bool loadSnpReferenceFile(const AffxString& strFileName);

        void loadInformationVector(std::vector<float>& vInformation);

        void loadPositionVector(std::vector<int> &vPosition);

        void computeGlobalSnpParameters();

        bool lohPreProcessing(  std::vector<float>& vSCAR,
                                std::vector<float>& vInformation,
                                std::vector<int>& vPosition,
                                AffxMultiDimensionalArray<double>& dHomValueArray,
                                AffxMultiDimensionalArray<double>& dHetValueArray );

        bool lohFind(   const std::vector<float>& vSCAR,
                       const std::vector<float>& vInformation,
                       const std::vector<int>& vPosition,
                       std::vector<char>& vLoh,
                       AffxMultiDimensionalArray<double>& dHomValueArray,
                       AffxMultiDimensionalArray<double>& dHetValueArray,
                       std::vector<int>& vWindowStart,
                       std::vector<int>& vWindowLength,
                       std::vector<double>& vWindowTest,
                       std::vector<double>& vWindowInflate,
                       std::vector<char>& vWindowLoh);

        inline double homDistribution(   double ecks,
                                double probabilityAA,
                                double probabilityBB,
                                double muPrimeAA,
                                double sigmaPrimeAA,
                                double muPrimeBB,
                                double sigmaPrimeBB);

        inline double hetDistribution(   double ecks,
                                double probabilityAA,
                                double probabilityBB,
                                double probabilityAB,
                                double muPrimeAA,
                                double sigmaPrimeAA,
                                double muPrimeBB,
                                double sigmaPrimeBB,
                                double muPrimeAB,
                                double sigmaPrimeAB);

        inline float likelihood(int iStart,
                                int iEnd,
                                AffxMultiDimensionalArray<double>& dHomValueArray,
                                AffxMultiDimensionalArray<double>& dHetValueArray,
                                int iVarianceInflationIndex);

        inline float sumVector( int iStart,
                                int iEnd,
                                const std::vector<float>& vInputVector);

        int newSegments(    int iSegmentType,
                                const std::vector<float>& vSCAR,
                                const std::vector<double>& vHomValues,
                                const std::vector<double>& vHetValues);

        inline float computeConfidence(int IStart,
                                       int iEnd,
                                       const std::vector<float>& vSCAR,
                                       const std::vector<double>& vHomValues,
                                       const std::vector<double>& vHetValues);

        float logLikelihood(    std::vector<float>& vfEcks,
                                float fMu,
                                float fSigma,
                                float fDelta,
                                float fProbability1,
                                float fProbability2,
                                float fProbability3);


        bool computeGlobalParameters(   const std::vector<float>& vSCAR,
                                        double& fMu,
                                        double& fSigma,
                                        double& fDelta,
                                        double& fProbability1,
                                        double& fProbability2,
                                        double& fProbability3);

        bool sumLogsInVector(  const std::vector<double>& vfZed1,
                                const std::vector<double>& vfZed2,
                                const std::vector<double>& vfZed3,
                                double& fSumOfLogs);

        void computeHomHetValues(        const std::vector<float>& vSCAR,
                                         AffxMultiDimensionalArray<double>& vHetValues,
                                         AffxMultiDimensionalArray<double>& vHomValues );

        int determineVarianceInflationIndex(    int iStart,
                                                int iEnd,
                                                const std::vector<float>& vInformation );

        void writeLOH(    std::vector<float>& vSCAR,
                        std::vector<int>& vGenomeChunkStart,
                        std::vector<int>& vGenomeChunkLength,
                        std::vector<int>& vWindowStart,
                        std::vector<int>& vWindowLength,
                        std::vector<double>& vWindowTest,
                        std::vector<double>& vWindowInflate,
                        std::vector<char>& vWindowLoh );

};

#endif
