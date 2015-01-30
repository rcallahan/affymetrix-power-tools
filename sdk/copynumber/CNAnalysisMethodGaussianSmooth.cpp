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
 * @file CNAnalysisMethodGaussianSmooth.cpp
 *
 * @brief This file contains the CNAnalysisMethodGaussianSmooth class members.
 */
#include "copynumber/CNAnalysisMethodGaussianSmooth.h"

CNAnalysisMethodGaussianSmooth::CNAnalysisMethodGaussianSmooth() { }

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodGaussianSmooth::explainSelf()
    {
    CNAnalysisMethodGaussianSmooth obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodGaussianSmooth::
getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt expSmoothSignal = {"expSmoothSignal", SelfDoc::Opt::Integer,
        "1", "1", "NA", "NA", "Exp Smooth Signal"};
  opts.push_back(expSmoothSignal);

  SelfDoc::Opt smooth_sigma_multiplier = {"smooth_sigma_multiplier",
        SelfDoc::Opt::Double, "2", "2", "0.0", "NA", "Smooth Sigma Multiplier"};
  opts.push_back(smooth_sigma_multiplier);

  SelfDoc::Opt smooth_window = {"smooth_bw", SelfDoc::Opt::Integer,
        "50000", "50000", "1", "NA", "Smooth Gaussian BW"};
  opts.push_back(smooth_window);

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
SelfCreate * CNAnalysisMethodGaussianSmooth::
newObject(std::map<std::string,std::string> & params)
    {
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();

    CNAnalysisMethodGaussianSmooth * pMethod;
    pMethod = new CNAnalysisMethodGaussianSmooth();

    std::string strPrefix = getPrefix();

    pMethod->m_iExpSmoothSignal = setupIntParameter("expSmoothSignal",
        strPrefix, params, doc, opts);

    pMethod->m_fSmoothSigmaMultiplier = setupFloatParameter(
        "smooth_sigma_multiplier", strPrefix, params, doc, opts);

    pMethod->m_smooth_bw = setupIntParameter("smooth_bw",
        strPrefix, params, doc, opts);

    return pMethod;
    }

/**
 * @brief Run the analysis
 */
void CNAnalysisMethodGaussianSmooth::run()
    {
    Verbose::out(1, "CNAnalysisMethodGaussianSmooth::run(...) start");
    isSetup();
    m_pEngine->setOpt("gaussian-smooth-exp", ::getInt(m_iExpSmoothSignal));

    int last_chr = getProbeSets()->at(getProbeSets()->size() - 1)->getChromosome();

//    Verbose::progressBegin(1, "CNAnalysisMethodGaussianSmooth::run(...) ",
//        last_chr, 1, last_chr);

    for (int chr=1; (chr<=last_chr); chr++)
        {
        int iProbeSetCount = getProbeSetCount(chr, getProbeSets());
        if (iProbeSetCount == 0) continue;
//        Verbose::progressStep(1);
        try
            {
            std::vector<int> vPositions(iProbeSetCount);
            std::vector<float> vLog2Ratios(iProbeSetCount);
            int iIndex = 0;
            for (int i=0; i<getProbeSets()->size(); i++)
                {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(i);
                if (pobjProbeSet->getChromosome() == chr)
                    {
                    vPositions[iIndex] = pobjProbeSet->getPosition();
                    vLog2Ratios[iIndex] = pobjProbeSet->getLog2Ratio();
                    iIndex++;
                    }
                }

            if (m_iExpSmoothSignal >= 0)
                {
                calculateGaussianSmooth(chr, vPositions, vLog2Ratios);
                }
            }

        catch(...)
            {
            Verbose::out(1, "CNAnalysisMethodGaussianSmooth::run(...) Working on Chromosome " + ::getInt(chr));
            throw;
            }

        }
    normalizeGaussianSmooth();
//    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodGaussianSmooth::run(...) end");
    }

/**
 * @brief Smooth the log2 ratios
 * @param int - Chromosome to process
 * @param bool - Is this an chromosome an autosome
 * @param std::vector<int>& - The vector of chromosome positions
 * @param std::vector<float>& - The vector of log2 ratios to smooth
 */
void CNAnalysisMethodGaussianSmooth::calculateGaussianSmooth(
        int iChromosome,
        std::vector<int>& vPositions,
        std::vector<float>& vLog2Ratios)
    {
    try
        {
        smooth(vLog2Ratios, vPositions, m_smooth_bw,
                m_fSmoothSigmaMultiplier, 0, (int)vLog2Ratios.size());
        }

    catch (Except e) { throw; }
    int iIndex = 0;

    for (int i=0; i<getProbeSets()->size(); i++)
        {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(i);
        if (pobjProbeSet->getChromosome() == iChromosome)
            {
            pobjProbeSet->setSmoothedLog2Ratio(vLog2Ratios[iIndex]);
            iIndex++;
            }
        }
    }

/**
 * @brief Given a list of log2 ratios and the corresponding positions on
 * a given chromosome this method will do gaussian smoothing on the
 * data.  This will smooth all positions with the unsmoothed data to
 * avoid biasing the smoothing process.  It would be biased because
 * the left side would be smoothed and the right side wouldn't be.
 *
 * @param lrvals   - Log2 ratios to smooth for a given sample on a given
 * chromosome
 * @param position - The SNP positions corresponding to the log2 ratios ("lrvals").
 * @param bw       - The bandwidth to use for the smoothing process.
 * @param startidx - The start position in the position vector.
 * @param sz       - The number of positions that should be processed.
 */
void CNAnalysisMethodGaussianSmooth::smooth( std::vector<float>& lrvals, const std::vector<int>& position,
                                             int smooth_bw,
                                             float sigmaMultiplier, int startidx, int sz)
{
  if ( !(lrvals.size()==position.size()) || !(smooth_bw > 0) || !(sz<=(int)lrvals.size())
    || !(startidx<sz) ) {
    GSERR_CHECK(lrvals.size()==position.size(), "Mismatch between the number of log2 ratios and the number of positions: #log2 ratios = " + ToStr(lrvals.size()) + "  # positions = " + ToStr(position.size()));
    GSERR_CHECK(smooth_bw>0, "The gaussian-smooth parameter  smooth_bw must be > 0");
    GSERR_CHECK(sz<=(int)lrvals.size(),"end index past end of array  end = "
               + ToStr(sz-1) + "  array length = " + ToStr(lrvals.size()));
    GSERR_CHECK(startidx<(int)sz, "start index " + ToStr(startidx)
              + " is larger than the end index " + ToStr(sz-1));
  }

  // Need more than 1 SNP to do smoothing
  if (lrvals.size() == 1) return;

  float smoothDist = smooth_bw * sigmaMultiplier;
  // The original log2 ratios must be kept track of since they are needed
  // for smoothing.  The smoothed log2 ratios will be kept in "tmplr"
  // until smoothing is done, then they will be copied to "lrvals".
  float* tmplr = new float[sz];         // Holds the smoothed log2 ratios.
  int iIndex = 0;
  int     lbidx     = 0; // Lower bound
  int     ubidx     = 0; // Upper bound
  int     curridx   = 0; // Index of current SNP being Smoothed.
  long    dist      = 0; // Distance between SNPs.
  double numerator = 0.0, denominator = 0.0;
  int idx = 0;
  for (curridx = startidx; curridx < sz; ++curridx) {
    // Compute lower bandwidth index
    while ( (dist = position[curridx] - position[lbidx]) > smoothDist ) {
      ++lbidx;
    }
    // Compute upper bandwidth index
    ubidx = curridx + 1; // Set to the position after the current pos.
    while (    ubidx < sz
           && (dist = position[ubidx] - position[curridx]) <= smoothDist) {
      ++ubidx;
    }
    /////////////////////////////////////////////////////////
    // Do smoothing starting from the "lbidx" to the "ubidx".
    /////////////////////////////////////////////////////////
    // Compute numerator and denominator.
    numerator   = 0.0;
    denominator = 0.0;
    // Loop from the lower bound ("lbidx") to the upper bound ("ubidx")
    // updating the numerator and denominator.
    for (idx = lbidx; idx < ubidx; ++idx) {
      double posdiff = position[idx]-position[curridx];
      double t       = posdiff/smooth_bw;
      double expon   = exp(-0.5*(t*t));
      numerator   += lrvals[idx] * expon;
      denominator += expon;
    }
    tmplr[iIndex] = (float)(numerator/denominator);
    iIndex++;
  }
  // Copy smoothed values to original log2 ratio vector.
  for (int idx = 0; idx < iIndex; ++idx) {
    lrvals[idx+startidx] = tmplr[idx];
  }
  delete[] tmplr;
} // smooth(...)



/**
 * @brief Normalize the smoothed values
 */
void CNAnalysisMethodGaussianSmooth::normalizeGaussianSmooth()
    {
        // This is now being done ni the cnreporter classes.
    if (m_iExpSmoothSignal != 1) return;
    /*
    double scale = 1.0/m_pEngine->getOptDouble("alpha-cn-calibrate");
    double shift = m_pEngine->getOptDouble("beta-cn-calibrate");

    for (int i=0; i<getProbeSets()->size(); i++)
        {
        CNProbeSet * pobjProbeSet = getProbeSets()->at(i);
        double log_intensity = pobjProbeSet->getSmoothedLog2Ratio();
        log_intensity *= scale;
        log_intensity += shift;

        pobjProbeSet->setSmoothedLog2Ratio((float)exp(log_intensity*log(2.0)));
        }
        */
    }
