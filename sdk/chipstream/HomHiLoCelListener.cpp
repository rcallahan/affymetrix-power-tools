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

//
#include "chipstream/HomHiLoCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
#include "stats/stats.h"
#include "util/Fs.h"


using namespace std;
using namespace affx;
using namespace affymetrix_fusion_io;

/** 
 * Process another cel files worth of data.
 */
double HomHiLoCelListener::computeStatistic(std::vector<double>& vData, double dEmThreshold, double dBinSize) {

	dBinSize = dBinSize / 2.0;
  // Find the 3 peaks in contrast values (AA, AB, BB)
  double peak1, peak2, peak3;
  computeClusterPeaks(vData,peak1,peak2,peak3,dEmThreshold, 0.10);
  
  Verbose::out(4,"HomHiLo - Peaks: " + ToStr(peak1) + ", " + ToStr(peak2) + ", " + ToStr(peak3));

  // Find the valleys
  double valley1 = computeClusterValleyFixedEndPoint(vData,peak1,peak2,dBinSize);
  double valley2 = computeClusterValleyFixedEndPoint(vData,peak2,peak3,dBinSize);
  
  Verbose::out(4,"HomHiLo -  Valleys: " + ToStr(valley1) + ", " + ToStr(valley2));

  // How many contrast values under each peak
  double densityPeak1 = getContrastDensity(vData,peak1,dBinSize);
  double densityPeak3 = getContrastDensity(vData,peak3,dBinSize);
  
  Verbose::out(4,"HomHiLo -  densityPeaks: " + ToStr(densityPeak1) + ", " + ToStr(densityPeak3));

  // How many contrast values under each valley
  double densityValley1 = getContrastDensity(vData,valley1,dBinSize);
  double densityValley2 = getContrastDensity(vData,valley2,dBinSize);
  
  Verbose::out(4,"HomHiLo - densityValleys: " + ToStr(densityValley1) + ", " + ToStr(densityValley2));

  // Compute Hom Hi Low
  double minHomHiLo = 0.1;
  double mn1 = densityPeak1 - densityValley1;
  double mn2 = densityPeak3 - densityValley2;

  // Get Min
  if (mn1 < mn2)
   minHomHiLo = mn1;
  else
   minHomHiLo = mn2;

  return minHomHiLo;
}

/** 
 * Process another cel files worth of data.
 */
void HomHiLoCelListener::newChip(FusionCELData *cel) {
  std::vector<ChipSummary::Metric> metrics;
  m_SummaryStats.push_back(metrics);
  m_CelNames.push_back(cel->GetFileName());

  // Get Contrast Values
  vector<double> contrastValues;
  fillInContrastValues(cel, m_ProbeSets, contrastValues, m_K);

  // Find the 3 peaks in contrast values (AA, AB, BB)
  double peak1, peak2, peak3;
  computeClusterPeaks(contrastValues,peak1,peak2,peak3,m_EmThresh);
  Verbose::out(4,"HomHiLo - " + Fs::basename(cel->GetFileName()) + " -  " + m_Label + 
          " - Peaks: " + ToStr(peak1) + ", " + ToStr(peak2) + ", " + ToStr(peak3));

  // Find the valleys
  double valley1 = computeClusterValley(contrastValues,peak1,peak2,m_BinSize);
  double valley2 = computeClusterValley(contrastValues,peak2,peak3,m_BinSize);
  Verbose::out(4,"HomHiLo - " + Fs::basename(cel->GetFileName()) + " -  " + m_Label + 
          " - Valleys: " + ToStr(valley1) + ", " + ToStr(valley2));

  // How many contrast values under each peak
  double densityPeak1 = getContrastDensity(contrastValues,peak1,m_BinSize);
  double densityPeak3 = getContrastDensity(contrastValues,peak3,m_BinSize);
  Verbose::out(4,"HomHiLo - " + Fs::basename(cel->GetFileName()) + " -  " + m_Label + 
          " - Peak Density: " + ToStr(densityPeak1) + ", " + ToStr(densityPeak3));

  // How many contrast values under each valley
  double densityValley1 = getContrastDensity(contrastValues,valley1,m_BinSize);
  double densityValley2 = getContrastDensity(contrastValues,valley2,m_BinSize);
  Verbose::out(4,"HomHiLo - " + Fs::basename(cel->GetFileName()) + " -  " + m_Label + 
          " - Valley Density: " + ToStr(densityValley1) + ", " + ToStr(densityValley2));

  // Compute Hom Hi Low
  double minHomHiLo = 0.1;
  double mn1 = densityPeak1 - densityValley1;
  double mn2 = densityPeak3 - densityValley2;

  // Get Min
  if (mn1 < mn2)
   minHomHiLo = mn1;
  else
   minHomHiLo = mn2;

  Verbose::out(4,"HomHiLo - " + Fs::basename(cel->GetFileName()) + " -  " + m_Label + 
          " - mn1,mn2,min: " + ToStr(mn1) + ", " + ToStr(mn2) + ", " + ToStr(minHomHiLo));

  // Populate the Summary Metrics
  m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_Label, minHomHiLo));

  ///@todo option to report other metrics?

  setValid(true);
}

/** 
 * Loop through the probesets provided and calculate a contrast value
 * for each one using the median of PM probes for A allele and B
 * allele.
 */
///@todo refactor? simmilar code in EmGenderCelListener
void HomHiLoCelListener::fillInContrastValues(FusionCELData *cel,
                                              std::vector<ProbeListPacked> &probeSets, 
                                              std::vector<double> &contrastValues,
                                              double k)
{
  contrastValues.clear();
  contrastValues.reserve(probeSets.size());
  for(int psIx = 0; psIx < probeSets.size(); psIx++) {
    const ProbeSet *ps =  ProbeListFactory::asProbeSet(probeSets[psIx]);
    contrastValues.push_back(CalculateEMContrast(cel, ps, k));
    delete ps;
  }
}

/**
 * Calculate contrast values
 */
///@todo refactor? simmilar code in EmGenderCelListener
double HomHiLoCelListener::CalculateEMContrast(FusionCELData* cel, const ProbeSet* ps, const double k) {
    float Amedian = -1, Bmedian = -1;
    bool medianOk = alleleMedians(cel, ps, Amedian, Bmedian);
    if(!medianOk || Amedian <= 0 || Bmedian <= 0) {
      Err::errAbort("HomHiLoCelListener::fillInContrastValues() - Warning. alleleMedian() failed for probeset: " + 
              ToStr(ps->name));
    }
    double contrast = -2, strength = -2;
    ContrastExtremesStretch(k, Amedian, Bmedian, contrast, strength);

    if(fabs(contrast) > 1) {
      Err::errAbort("HomHiLoCelListener::fillInContrastValues() - Can't have abs(contrast) > 1 for probeset: " + 
              ToStr(ps->name));
    }
    return contrast;
}

/**
 * Compute the cluster peaks using EM
 */
void HomHiLoCelListener::computeClusterPeaks(vector<double> &contrastValues, double &peak1, double &peak2, double &peak3, double emThresh, double dMinMu) {
    CEMSeed seed;
    seed.setMu(-0.66f, 0.0f, 0.66f);
    seed.setSigma(0.1f, 0.1f, 0.1f);
    seed.setWeight(0.33f, 0.34f, 0.33f);
    seed.setMinMu(-2.0f, -.05f, dMinMu);
    seed.setMinSigma(0.02f, 0.02f, 0.02f);
    seed.setMaxMu(-dMinMu, .05f, 2.0f);
    seed.setMaxSigma(.3f, .3f, .3f);

    CPrimeEM em;
    em.setData(contrastValues);
    em.setThreshold(emThresh);
    em.EMEstimate(seed);
    CEMEst* pEst = em.getEMEstimates();

    peak1 = pEst->m_mu[0];
    peak2 = pEst->m_mu[1];
    peak3 = pEst->m_mu[2];
}

/**
 * Compute the Valley Between Two Points (The corrected version.)
 */ 
double HomHiLoCelListener::computeClusterValleyFixedEndPoint(const std::vector<double>& contrast, double x1, double x2, double bin){
    double valleyX = x1;
    double valleyDensity = getContrastDensity(contrast, valleyX, bin);
    double x = x1;
    while (x < x2) {
        x += 2.0 * bin;
		if (x > x2) {x = x2;}
        double density = getContrastDensity(contrast, x, bin);
        if (density < valleyDensity){
            valleyX = x;
            valleyDensity = density;
        }
    }
    return valleyX;
}

/**
 * Compute the Valley Between Two Points
 */ 
double HomHiLoCelListener::computeClusterValley(const std::vector<double>& contrast, double x1, double x2, double bin){
    double valleyX = x1;
    double valleyDensity = getContrastDensity(contrast, valleyX, bin);
    double x = x1 + 2.0 * bin;
    while (x < x2) {
        double density = getContrastDensity(contrast, x, bin);
        if (density < valleyDensity){
            valleyX = x;
            valleyDensity = density;
        }
        x += 2.0 * bin;
    }
    return valleyX;
}

/**
 * Compute the density of values for a given bin
 */
double HomHiLoCelListener::getContrastDensity(const std::vector<double>& contrast, double x, double bin){
    double lb = x - bin;
    double ub = x + bin;
    size_t n = contrast.size();
    int count = 0;
    for (size_t i=0; i < n; i++) {
        if ( (contrast[i]> lb) && (contrast[i] <= ub))
            count++;
    }
    return 100.0*(double(count)/double(n));
}


// Need to undef min macro to enable the STL numeric_limits min().
// Remember this will remain in effect until the end of the file.
#undef min


void HomHiLoCelListener::computeMixtureModel(
                                  const vector<double>& xdata,
                                  double& sigma,
                                  double& delta,
                                  double& P,
                                  double& K,
                                  bool symmetrize,
                                  double sigmaInit,
                                  double deltaInit,
                                  double PInit,
                                  double KInit)
{
    const int maxIter = 200;
    const double epsilon = 5.0e-7;
    const double nearZero = std::numeric_limits<double>::min();

    vector<double> z1(xdata.size());
    vector<double> z2(xdata.size());
    vector<double> z3(xdata.size());

    int dataSize = xdata.size();
    if (symmetrize) {
        dataSize *= 2;
    }

    sigma = sigmaInit;
    delta = deltaInit;
    P = PInit;
    K = KInit;
    double logL = logLikelihood(xdata, sigma, delta, P, K);

    double logLprev;
    int iCount = 0;

    do {
        for (int i = 0; i < xdata.size(); i++) {
            double dd = mixSum(xdata[i], sigma, delta, P, K);
            if (dd < nearZero) {
                z2[i] = 0.0;
                if (xdata[i] > 0.0) {
                    z1[i] = 1.0;
                    z3[i] = 0.0;
                } else {
                    z1[i] = 0.0;
                    z3[i] = 1.0;
                }
            } else {
                z1[i] = P/2.0 * dnorm(xdata[i], sigma, delta)/dd;
                z2[i] = (1.0 - P) * dnorm(xdata[i], K*sigma, 0.0)/dd;
                z3[i] = P/2.0 * dnorm(xdata[i], sigma, -delta)/dd;
            }
        }
        double A1 = 0.0;
        double A2 = 0.0;
        double A3 = 0.0;
        double M1 = 0.0;
        double M3 = 0.0;
        for (int i = 0; i < xdata.size(); i++) {
            A1 += z1[i];
            A2 += z2[i];
            A3 += z3[i];
            M1 += z1[i]*xdata[i];
            M3 += z3[i]*xdata[i];
        }
        if (symmetrize) {
            // The following amounts to doing the loop above with xdata[i] replaced by -xdata[i].
            A1 = A1 + A3;
            A2 += A2;
            A3 = A1;
            M1 = M1 - M3;
            M3 = -M1;
        }

        delta = (M1 - M3)/(A1 + A3);

        double M11 = 0.0;
        double M21 = 0.0;
        double M31 = 0.0;
        for (int i = 0; i < xdata.size(); i++) {
            M11 += z1[i]*(xdata[i] - delta)*(xdata[i] - delta);
            M21 += z2[i]*xdata[i]*xdata[i];
            M31 += z3[i]*(xdata[i] + delta)*(xdata[i] + delta);
        }
        if (symmetrize) {
            // The following amounts to doing the loop above with xdata[i] replaced by -xdata[i].
            M11 = M11 + M31;
            M21 += M21;
            M31 = M11;
        }

        // maxima (critical points)
        P = (A1 + A3)/dataSize;
        sigma = sqrt((M11 + M21/(K*K) + M31)/dataSize);
        K = sqrt(M21/A2)/sigma;

        if (K > 1.0) {
            K = 1.0;
        }
        logLprev = logL;
        logL = logLikelihood(xdata, sigma, delta, P, K);

        iCount++;
    } while (iCount < maxIter && logL - logLprev > epsilon);
}

float HomHiLoCelListener::computeSNPQC(const vector<double>& xdata, bool symmetrize)
{
    double sigma;
    double delta;
    double P;
    double K;
    if (xdata.empty()) {
        return numeric_limits<float>::quiet_NaN();
    }
    computeMixtureModel(xdata, sigma, delta, P, K, symmetrize);

    return (delta*delta)/(sigma*sigma*(K*K + 1));
}

double HomHiLoCelListener::logLikelihood(
                              const std::vector<double>& xdata,
                              double sigma,
                              double delta,
                              double P,
                              double K)
{
    double result = 0.0;
    for (int i = 0; i < xdata.size(); i++) {
        result += log(mixSum(xdata[i], sigma, delta, P, K));
    }
    // This works for both symmetrized and unsymmetrized case
    // since it's either result/size or (2*result)/(2*size).
    result /= xdata.size();

    return result;
}

double HomHiLoCelListener::mixSum(
                              double x,
                              double sigma,
                              double delta,
                              double P,
                              double K)
{
    double result = P/2.0 * dnorm(x, sigma, delta) +
                    (1.0 - P) * dnorm(x, K*sigma, 0.0) +
                    P/2.0 * dnorm(x, sigma, -delta);
    return result;
}

double HomHiLoCelListener::dnorm(
                              double x,
                              double sigma,
                              double delta)
{
    const double pi = 3.14159265358979323846;
    double result = exp(-0.5*(((x - delta)/sigma)*((x - delta)/sigma)));
    result /= sqrt(2.0 * pi);
    result /= sigma;

    return result;
}
