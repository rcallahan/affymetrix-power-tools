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

#ifndef _pdnn_h_
#define _pdnn_h_
/**
 * @file pdnn.h
 *
 * @brief This header contains the pdnn class definition.
 */

#define NEWMAT_USING_FLOAT
#include "copynumber/CNAnalysisMethod.h"
//
#include "util/AffxArray.h"
//
#include <vector>
//

class PDNNProbe
{
public:
    int ProbeID;
    double ProbeEffect;
    int ProbeSetIndex;
    char Allele;
    AffxString Sequence;
    float MedianIntensity;
    float PredictedIntensity;

    PDNNProbe()
    {
        ProbeID = 0;
        ProbeEffect = 0;
        ProbeSetIndex = -1;
        Allele = 0;
        MedianIntensity = 0;
        PredictedIntensity = 0;
    }

    int compareTo(PDNNProbe& that, int iCompareCode)
    {
        int iCompareResult = 0;
        switch (iCompareCode)
        {
        case 0:
            iCompareResult = AffxArray<int>::compare(ProbeID, that.ProbeID);
            break;
        }
        return iCompareResult;
    }

    template<int k> struct ComparePred {
        bool operator()(const PDNNProbe* lhs, const PDNNProbe* rhs) const {
            Err::errAbort("PDNNProbe: ComparePred instantiated with an invalid compare code = " + ToStr(k));
            return false;
        }
    };
};

template<> struct PDNNProbe::ComparePred<0> {
    bool operator()(const PDNNProbe* lhs, const PDNNProbe* rhs) const {
        return AffxArray<int>::compare(lhs->ProbeID, rhs->ProbeID) < 0;
    }
};

/**
 * @brief  A class for calulating pdnn.
 *
 */
class CNReferenceMethodPDNN : CNAnalysisMethod
{
public:
    CNReferenceMethodPDNN();
    virtual ~CNReferenceMethodPDNN();

    static std::string getType() { return "pdnn-reference-method"; }
    static std::string getDescription() { return "CopyNumber PDNN"; }
    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string> & param);

    virtual AffxString getName() { return getType(); }
    virtual void run();
    virtual bool isSegmentTypeAnalysis() { return false; }

private:
    void loadProbes(const std::string& strReferenceFileName, const std::string& str1lqFileName, AffxArray<PDNNProbe>& vPDNNProbes);
    void getRandomSample(int iPopulationSize, int iSampleSize, std::vector<int>& vRandomSample);
    char getDinucleotideType(char c1, char c2);
    void svd(Matrix& mxX, DiagonalMatrix& mxDiagonal, Matrix& mxU, Matrix& mxV, Matrix& mxData, ColumnVector& vResult);
    void invertDiagonal(DiagonalMatrix& mxDiagonal);
    void flipSignIfNeeded(ColumnVector& vWeights, ColumnVector& vEnergy, Matrix& mxX);
    bool hasConverged(unsigned int uiProbeLength, int iIteration, ColumnVector& vSampledIntensities, Matrix& X, ColumnVector& vEnergy, ColumnVector& vPrevEnergy, ColumnVector& vWeights, ColumnVector& vPrevWeights, double& dPrevFit, double dTOL);
    void updateReference(const std::string& strReferenceFileName, AffxArray<PDNNProbe>& vPDNNProbes);
};

#endif // _pdnn_h_


