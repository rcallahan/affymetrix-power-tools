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

#ifndef _CNAnalysisMethod_H_
#define _CNAnalysisMethod_H_
/**
 * @file CNAnalysisMethod.h
 *
 * @brief This header contains the CNAnalysisMethod class definition.
 */

#include "copynumber/CNExperiment.h"
#include "copynumber/CNProbe.h"
#include "copynumber/CNProbeSet.h"
#include "copynumber/CNSegment.h"
//
#include "calvin_files/data/src/CHPMultiDataData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include "../external/newmat/newmatap.h"
//

#define CN_INVALID_DOUBLE  (-9999.9999)

using namespace std;

/**
 * @brief  A base class for copy number analysis methods.
 *
 */
class CNAnalysisMethod : public SelfDoc, public SelfCreate
{
private:
  static int m_iInstanceCount;
  bool m_boundsFilled;

  void fillChrBounds();

protected:
  BaseEngine* m_pEngine;
  int m_iXChromosome;
  int m_iYChromosome;
  CNExperiment* m_pobjExperiment;
  CNExperimentArray* m_pvExperiments;     // introduced for raw SNPQC calculation
  CNProbeSetArray* m_pvProbeSets;
  CNProbeArray* m_pvProbes;
  map<int, pair<int, int> > m_chrBounds;
  vector<double> m_vCoarseAllelePeaks;
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType> m_vCelFileParams;
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType> m_vParams;

  static const int allelePeakCount = 6;
  static const int CytoScanHDAllelePeakCount = 12;
  static const int iMaxPeak = 3;
  static const int dPeakMultiplier = 170;

protected:
  static AffxString getPrefix();
  static SelfDoc::Opt* getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName);
  static bool setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static int setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static float setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static double setupDoubleParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static std::string setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);

  static void doBitReverseOrder(vector<double>& vec);
  static int reverseBits(int ii, int nbits);
  static void doFourierTransform(vector<double>& re, vector<double>& im, int expSign);

  struct binCompare : std::binary_function<std::pair<int, int>, std::pair<int, int>, bool>
  {
      bool operator()(const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) const {
          return lhs.first < rhs.first;
      }
  };

  struct comparex : std::binary_function<std::vector<double>, std::vector<double>, bool>
  {
      bool operator()(const std::vector<double> &a, const std::vector<double> &b) const {
          return (a[0] < b[0]);
      }
  };

protected:
  unsigned int m_uiDataSetOffset;
  CNSegmentArray m_vSegments;

  // for peak shrinking
  int m_iStep;
  int m_iWindow;
  int m_iPointCount;
  double m_dBandwidth;
  double m_dCutoff;
  double m_dCleanThreshold;
  bool m_bSymmetry;

  // alternate parameters for custom peak shrinking
  struct PeakShrinkOverride {
      int m_iStep_override;
      int m_iWindow_override;
      int m_iPointCount_override;
      double m_dBandwidth_override;
      double m_dCutoff_override;
      double m_dCleanThreshold_override;
      bool m_bSymmetry_override;
      int m_iCovariateIndex;
  };

  vector<float> m_vThreePeakFLD_X;
  vector<float> m_vThreePeakFLD_Y;
  vector<float> m_vFourPeakFLD_X;
  vector<float> m_vFourPeakFLD_Y;

  vector<float> m_vThreePeakShrink_X;
  vector<float> m_vThreePeakShrink_Y;
  vector<float> m_vFourPeakShrink_X;
  vector<float> m_vFourPeakShrink_Y;

private:
    // for deep regression
    struct apParameters {
        int window;
        int step;
        double bandwidthFactor;
        int densityPointCount;
        double peaksCutoff;
        double cleanThreshold;
        bool symmetry;
        int covariateIndex;
    };
    struct cutoffShrinkage {
        double FLD3threshold;
        double FLD4threshold;
        double shrinkage3;
        double shrinkage4;
        double SNPQC;
    };
    //

public:
  CNAnalysisMethod();
  virtual ~CNAnalysisMethod();

  virtual AffxString getName() = 0;
  int getSegmentType() ;
  AffxString getSegmentTypeString();

  virtual bool isSegmentTypeAnalysis() = 0;
  CNSegmentArray& getSegments();
  CNSegmentArray getSegments(int iType);
  int getSegmentCount(int iType);

  virtual void setEngine(BaseEngine* p);
  BaseEngine* getEngine();
  CNExperiment* getExperiment();
  CNExperimentArray* getExperiments();
  CNProbeSetArray* getProbeSets() { return m_pvProbeSets; }
  void setProbes(CNProbeArray& vProbes);
  CNProbeArray* getProbes();
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType>* getCelFileParams();
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType>* getParams();

  static float getConfidenceThreshold(const std::string& brlmmpStr);

  void setup(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets, CNProbeArray* pvProbes=NULL);
  void setup(CNExperimentArray& vExperiments, int experimentIndex, CNProbeSetArray& vProbeSets);
  void setup(CNExperimentArray& vExperiments, CNProbeSetArray& vProbeSets);
  void isSetup();

  virtual void run() = 0;

  void setDataSetOffset(unsigned int ui);
  unsigned int getDataSetOffset();

  pair<int, int> getChrBounds(const int chr, CNProbeSetArray* pvProbeSets);
  std::pair<int, int> getChrBounds(int chr);
  std::vector<int>getChromosomes(CNProbeSetArray*);

protected:
  bool isNaN(const double& d) { return d != d; }
  int getProbeSetCount(int iChromosome, CNProbeSetArray* );
  void bin(std::vector<float>& vValues, std::vector<int>& vBinIndexes, int iBinCount);
  void writeIntensities(std::string intensityFileInfix);
  void writeSignals(std::string signalFileInfix);
  void writeLog2Ratios(std::string l2rFileInfix);
  void writeSmoothedLog2Ratios(std::string l2rFileInfix, CNProbeSetArray* vProbeSets);
  void writeLog2RatiosTsv(std::string l2rFileInfix, CNProbeSetArray* pvLocalProbeSets);
  void writeSignalsTsv(std::string signalFileInfix, CNProbeSetArray* pvLocalProbeSets);
  void writeMediansVectorAndMedian(std::string infix, const std::vector<float>& medians, float grandMedian = 0.0);
  void writeMediansVector(std::string infix, const std::vector<float>& medians);
  void computeSCAR(std::vector<float>& vSCAR);

  void calculateWindows( vector<pair<int, int> >& windowBds, vector<pair<int, int> >& stepWindowBds, const pair<int, int>& chrBound, int iStep, int iWindow);
  void calculateWeights(CNProbeSetArray* probeSets, vector<double>& vValsToProcess, vector<double>& weights, const pair<int, int>& windowBds, bool symmetry);
  int findpeaks(vector<int> &valleys, vector<int> &peaks, vector<double> &y, double delta, vector<double> &x);
  void findMaxPeaks(CNProbeSetArray* probeSets, const vector<pair<int, int> >& windowBds, const vector<vector<double> >& peaks);
  void shrinkToPeaks(CNProbeSetArray* probeSets, PeakShrinkOverride *paramOverride = NULL);
  void resetAllelePeakInitialValues();
  void fitDensityCurve( const std::vector<double>& values,
                        const std::vector<double>& weights,
                        std::vector<double>& xi,
                        std::vector<double>& density,
                        int numPoints,
                        double bandwidth,
                        bool symmetry
                        );

  // for deep regression
  void writeIntermediateData(
                    CNProbeSetArray* probeSets,
                    const vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& peakWindows,
                    const vector<pair<int, int> >& stepWindows,
                    const vector<pair<int, double> >& closestPeak,
                    const vector<pair<pair<int, int>, vector<double> > >& peakSets,
                    const vector<pair<int, double> >& processedValues,
                    const apParameters& apParametersValues,
                    const cutoffShrinkage& cutoffShrinkageValues
                    );
  //

  static float getPercentile(const vector<float>& vec, double dPercentile);
  static void binEqualNumber(
                    const std::vector<float>& vValues,
                    std::vector<int>& vBinIndexes,
                    int iBinCount);
  static void binEqualSpacing(
                    const std::vector<float>& vValues,
                    std::vector<int>& vBinIndexes,
                    int iBinCount);
  static int covariateIsBinAssignment(const std::vector<float>& vValues, std::vector<int>& vBinIndexes);

  virtual void filterCutOff(
                    CNProbeSetArray* probeSets,
                    int chrStart,
                    int chrEnd,
                    float cutoff3,
                    float cutoff4,
                    vector<bool>& filteredProbeSets
                    )
  {}

  int findcutoffInd(const vector<float>& cutoff);


  void fillChrBoundsImpl(CNProbeSetArray* probeSets,  map<int, pair<int, int> >& chrBounds) {
      // Assume that the genome is not being dynamically modified.
      if (m_boundsFilled) return;

      // Clear out the map just in case.
      chrBounds.clear();

      int chr = -1;
      for (int i = 0; i < probeSets->size(); i++) {
        // find out the chromosome number of this probe set
        //int this_chr = probeSets->operator[](i)->getChromosome();
        int this_chr = (*probeSets)[i]->getChromosome();

        // if it is the same as before increment the position of the end
        // and this iteration is done.
        if (this_chr == chr) {
          chrBounds[chr].second++;
          continue;
        }

        // If it is not the same as before, make sure the order is not
        // mangled.  The order of chromosomes is not important only
        // that they be contiguous.  If the order is a concern, check for
        // the order where the concern is motivated.
        if (chrBounds.find(this_chr) != chrBounds.end()) {
          printf("chr = %d\n", chr);
          printf("this_chr = %d\n", this_chr);
          throw(Except("CNAnalysisMethod::fillChrBounds() : "
                       "chromosomes are not contiguous\n"));
        }

        // Things look good (contiguous) so start a new chromosome.
        chrBounds[this_chr] = make_pair(i, i + 1);
        chr = this_chr;
      }

      m_boundsFilled = true;  // No need to do this again.
  }

    virtual void filterNoMansLand(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks
                    )
    {}

    virtual void filterShrinkTowardPeaks(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks,
                    float shrinkFactor3,
                    float shrinkFactor4,
                    vector<pair<int, double> >& closestPeak,
                    vector<pair<int, double> >& processedValues,
                    bool saveAllelePeaksFlag
                    )
    {}

    // These templates perform the actual filtering work and collect the code that's almost
    // identical for both allele peaks and allelic differences processing except for the calls
    // to getSCAR() (alelle peaks) and getAllelicDifference() (allelic differences). Since these
    // calls reside deep inside tight loops, it was decided not to use the virtual function
    // mechanism since it turns off inlining, but instead hide them under a statically dispatched
    // function T::getRequiredValue().
    //
    template <class T>
    void filterNoMansLandImpl(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks
                    )
    {
        for (int stepWinInd = 0; stepWinInd < stepWindowBds.size(); stepWinInd++) {
            for (int i = stepWindowBds[stepWinInd].first; i <= stepWindowBds[stepWinInd].second; i++) {
                if (!filteredProbeSets[i]) {
                    continue;
                }
                if (probeSets->getAt(i)->getMaxPeaks() != peaks[stepWinInd].size()) {
                    continue;
                }
                for (int j = 1; j < peaks[stepWinInd].size(); j++) {
                    double delta = (1.0 - m_dCleanThreshold)/2.0 * (peaks[stepWinInd][j] - peaks[stepWinInd][j-1]);
                    if (T::getRequiredValue(probeSets->getAt(i)) > peaks[stepWinInd][j-1] + delta  &&
                        T::getRequiredValue(probeSets->getAt(i)) < peaks[stepWinInd][j]   - delta)
                    {
                        filteredProbeSets[i] = false;
                    }
                }
            }
        }
    }

    template <class T>
    void filterShrinkTowardPeaksImpl(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks,
                    float shrinkFactor3,
                    float shrinkFactor4,
                    vector<pair<int, double> >& closestPeak,
                    vector<pair<int, double> >& processedValues,
                    bool saveAllelePeaksFlag
                    )
    {
        const bool keepIntermediateData = m_pEngine->getOptBool("keep-intermediate-data");
        for (int stepWinInd = 0; stepWinInd < stepWindowBds.size(); stepWinInd++) {
            if (peaks[stepWinInd].empty()) {
                continue;
            }
            for (int i = stepWindowBds[stepWinInd].first; i <= stepWindowBds[stepWinInd].second; i++) {
                if (!filteredProbeSets[i]) {
                    if (saveAllelePeaksFlag) {
                        // mark as invalid, do not use this value for density estimation later
                        m_vCoarseAllelePeaks[i] = std::numeric_limits<double>::quiet_NaN();
                    }
                    continue;
                }
                double minDist = std::fabs(T::getRequiredValue(probeSets->getAt(i)) - peaks[stepWinInd][0]);
                int minIndex = 0;
                for (int j = 1; j < peaks[stepWinInd].size(); j++) {
                    float dist = std::fabs(T::getRequiredValue(probeSets->getAt(i)) - peaks[stepWinInd][j]);
                    if (minDist > dist) {
                        minDist = dist;
                        minIndex = j;
                    }
                }

                if (keepIntermediateData) {
                    closestPeak.push_back(make_pair(i, peaks[stepWinInd][minIndex]));
                }

                float shrinkFactor;

                if (probeSets->getAt(i)->getMaxPeaks() <= 3) {
                    shrinkFactor = shrinkFactor3;
                } else {
                    shrinkFactor = shrinkFactor4;
                }

                float origValue = T::getRequiredValue(probeSets->getAt(i));

                float shrinkValue = origValue + shrinkFactor*(peaks[stepWinInd][minIndex] - origValue);

                if (keepIntermediateData) {
                    processedValues.push_back(make_pair(i, shrinkValue));
                }

                if (saveAllelePeaksFlag)
                {
                    m_vCoarseAllelePeaks[i] = shrinkValue;
                }
                else {
                    // Encode: restrict shrinkValue to [-iMaxPeak, iMaxPeak] for encoding
                    if (shrinkValue < -iMaxPeak) {
                        shrinkValue = -iMaxPeak;
                    }
                    else if (shrinkValue > iMaxPeak) {
                        shrinkValue = iMaxPeak;
                    }
                    unsigned int uiCompositeValue = (unsigned int)(::roundDouble(((shrinkValue + iMaxPeak) * (double)dPeakMultiplier), 0) + 1);
                    probeSets->getAt(i)->setAllelePeaks1(uiCompositeValue);
                }
            }
        }
    }

public:
  static void memory(const AffxString& str);

protected:

  void runmean(double *In, double *Out, const int *nIn, const int *nWin);

  double corr(Matrix& mx1, Matrix& mx2);

  double norm(Matrix& mx);

  double var(Matrix& mx);

  /**
   * Get the last autosome chromosome. As defined by being the last chromosome before the X chromosome value.
   * @return int - The chromosome number
   */
  int getLastAutosomeChromosome();

#define ELEM_TYPE double

  ELEM_TYPE klowest_select(ELEM_TYPE data[], int n, int klowest) ;

#undef ELEM_TYPE

  double percentile(double prctile, double* data, int m);

public:

  double bwnrd(const vector<double> &x,
               double fac);

  void kdensity(vector<double> &dat,
                vector<double> &density,
                vector<double> &xOut,
                vector<double> &weight,
                double bandWidth);

  double trapzoid(vector<double> &x,
                  vector<double> &y);

  double vectorMax(vector <double> &dat);

  double vectorMin(vector <double> &dat);

  double phi(double x);

  void calculateSummaryLOH();

  double assignConfidenceSigmoidal(     int iHetsInSegment,
                                        int iNumberInSegment,
                                        double dHetRate,
                                        double dPError,
                                        int iShift,
                                        int iMinSegmentSize);

  double shiftValue(    double ecks,
                        int iShift);
};

#endif


