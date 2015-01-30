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
 * @file CNAnalysisMethodMosaicism.cpp
 *
 * @brief This file contains the CNAnalysisMethodMosaicism class members.
 */
#include "copynumber/CNAnalysisMethodMosaicism.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//
#include "file/TsvFile/TsvFile.h"

// uncomment to dump data for mosaicism debugging
//#define MOSAICISM_DEBUG 1

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodMosaicism::explainSelf()
{
    CNAnalysisMethodMosaicism obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodMosaicism::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
    SelfDoc::Opt gains_boundries = {"gains-boundries", SelfDoc::Opt::String,
        "0.08764945,0.15380349,0.21465931,0.27100300", "0.08764945,0.15380349,0.21465931,0.27100300", "", "", "Mosaicism Gains Boundries"};
    opts.push_back(gains_boundries);

    SelfDoc::Opt losses_boundries = {"losses-boundries", SelfDoc::Opt::String,
        "-0.08293345,-0.17551812,-0.28048196,-0.40165383", "-0.08293345,-0.17551812,-0.28048196,-0.40165383", "", "", "Mosaicism Losses Boundries"};
    opts.push_back(losses_boundries);

    SelfDoc::Opt marker_bandwidth = {"marker-bandwidth", SelfDoc::Opt::Integer,
        "6000", "6000", "1", "NA",
        "Mosaicism Marker Bandwitdth"};
    opts.push_back(marker_bandwidth);

    SelfDoc::Opt confidence_window = {"confidence-window", SelfDoc::Opt::Integer,
        "251", "251", "31", "NA",
        "Mosaicism confidence running median window size"};
    opts.push_back(confidence_window);

    SelfDoc::Opt run_y_chromosome = {"run-y-chromosome", SelfDoc::Opt::Boolean,
        "true", "true", "NA", "NA",
        "Run Mosaicism analysis on Y Chromosome"};
    opts.push_back(run_y_chromosome);

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
SelfCreate* CNAnalysisMethodMosaicism::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodMosaicism* pMethod = new CNAnalysisMethodMosaicism();
    std::string strPrefix = getPrefix();

    AffxByteArray ba;
    AffxString str;

    // m_vGainsBoundries
    str = setupStringParameter("gains-boundries", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vGainsBoundries.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
    {
        pMethod->m_vGainsBoundries[k] = ba.getParameter(k+1).parseDouble();
    }

    // m_vLossesBoundries
    str = setupStringParameter("losses-boundries", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vLossesBoundries.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
    {
        pMethod->m_vLossesBoundries[k] = ba.getParameter(k+1).parseDouble();
    }

    // marker bandwidth
    pMethod->m_iMarkerBandwidth = setupIntParameter("marker-bandwidth", strPrefix, params, doc, opts);
    pMethod->m_iConfidenceWindow = setupIntParameter("confidence-window", strPrefix, params, doc, opts);

    if (pMethod->m_vGainsBoundries.size() != pMethod->m_vLossesBoundries.size())
    {
        throw(Except("losses-boundries must conform to gains-boundries."));
    }

    if ((pMethod->m_vGainsBoundries.size() != 4) || (pMethod->m_vLossesBoundries.size() != 4)) {throw(Except("gains-boundries and losses-boundries must contain four floating point numbers."));}

    for (int k=0; k<pMethod->m_vGainsBoundries.size(); k++)
    {
        if ((k > 0) && (pMethod->m_vGainsBoundries[k-1] > pMethod->m_vGainsBoundries[k])) {throw(Except("gains-boundries parameters must ordered increasing."));}
    }

    for (int k=0; k<pMethod->m_vLossesBoundries.size(); k++)
    {
        if ((k > 0) && (pMethod->m_vLossesBoundries[k-1] < pMethod->m_vLossesBoundries[k])) {throw(Except("losses-boundries parameters must ordered decreasing."));}
    }

    // run on Y chromosome
    pMethod->m_bRunYChromosome = setupBoolParameter("run-y-chromosome", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodMosaicism::CNAnalysisMethodMosaicism()
{
    m_iMarkerBandwidth = 0;
    m_iConfidenceWindow = 0;
    m_iStartIndex = 0;
    m_iEndIndex = 0;
    m_bLocalProbeSetsDetermined=false;
    m_bRunYChromosome=true;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodMosaicism::run()
{
    Verbose::out(1, "CNAnalysisMethodMosaicism::run(...) start");
    isSetup();
        determineLocalProbeSets();

    m_vSegments.deleteAll();
    vector<int> vChromosomes = getChromosomes(getProbeSets());
    int iLastChromosome = vChromosomes[vChromosomes.size() - 1];
    Verbose::progressBegin(1,"CNAnalysisMethodMosaicism::run(...) ", iLastChromosome, 1, iLastChromosome);

    for (int i = 0; (i < vChromosomes.size()); i++)
    {
        int iChromosome = vChromosomes[i];
        if (iChromosome == 255) {continue;}
        if (iChromosome == m_iYChromosome && !m_bRunYChromosome) {
          continue;
        }
        int iProbeSetCount = getProbeSetCount(iChromosome, getProbeSets());
        if (iProbeSetCount == 0) continue; // probably 23 or Y
        Verbose::progressStep(1);

        // debug: for later use
        m_current_experiment=getExperiment()->getExperimentName();
        m_current_chr=ToStr(vChromosomes[i]);

        std::vector<int> vPositions(iProbeSetCount);
        double* pLog2Ratios = new double[iProbeSetCount];
        int iStartChromosome = getChrBounds(iChromosome, getProbeSets()).first;
        m_iStartIndex = 0;
        m_iEndIndex = (iProbeSetCount - 1);
        CNProbeSetArray* psets = getProbeSets();
        for (int j = 0; (j < iProbeSetCount); j++)
        {
            vPositions[j] = psets->at(iStartChromosome + j)->getPosition();
            pLog2Ratios[j] = psets->at(iStartChromosome + j)->getLog2Ratio();
        }
        double* pRunMean = new double[iProbeSetCount];
        int iMarkerBandwidth = Min(iProbeSetCount, m_iMarkerBandwidth);
        runmean(pLog2Ratios, pRunMean, &iProbeSetCount, &iMarkerBandwidth);

        // debug - this has been ruled out.
#ifdef MOSAICISM_DEBUG_MEANS
        writeRunningMeans(m_current_experiment+".chr"+m_current_chr+".means.tsv",
                          iProbeSetCount,
                          pLog2Ratios,
                          pRunMean);
#endif

        CNSegmentArray vSegments;
        newSegments((unsigned char)iChromosome, pRunMean, iProbeSetCount, vPositions, vSegments);
        // dont need this anymore.
        delete[] pRunMean;

        // debug
#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".newseg1.tsv");
#endif
        //
        cleanSegments(vSegments);

        // Set marker count and median marker distance.
        for (int iSegmentIndex = 0; (iSegmentIndex < vSegments.getCount()); iSegmentIndex++)
        {
            CNSegment* p = vSegments.getAt(iSegmentIndex);
            p->setStartPosition(vPositions[p->getStartIndex()]);
            p->setEndPosition(vPositions[p->getEndIndex()]);
            calculateSegmentConfidence(pLog2Ratios, p);
            int iMarkerCount = 0;
            for (unsigned int iProbeSetIndex = 0; (iProbeSetIndex < vPositions.size()); iProbeSetIndex++)
            {
                if ((vPositions[iProbeSetIndex] >= p->getStartPosition()) && (vPositions[iProbeSetIndex] <= p->getEndPosition()))
                {
                    iMarkerCount++;
                }
            }
            if (iMarkerCount < 2)
            {
                p->setMarkerCount(iMarkerCount);
                p->setMeanMarkerDistance(0);
            }
            else
            {
                AffxMultiDimensionalArray<int> vMarkerDistances(iMarkerCount - 1);
                int iIndex = 0;
                for (unsigned int iProbeSetIndex = 1; (iProbeSetIndex < vPositions.size()); iProbeSetIndex++)
                {
                    int iPrevPosition = vPositions[iProbeSetIndex - 1];
                    int iNextPosition = vPositions[iProbeSetIndex];
                    if ((iPrevPosition >= p->getStartPosition()) && (iNextPosition <= p->getEndPosition()))
                    {
                        getProbeSets()->getAt(iStartChromosome + iProbeSetIndex - 1)->setMosaicismMixture(p->getMixtureAsDouble());
                        getProbeSets()->getAt(iStartChromosome + iProbeSetIndex)->setMosaicismMixture(p->getMixtureAsDouble());
                        vMarkerDistances.set(iIndex, (iNextPosition - iPrevPosition));
                        iIndex++;
                    }
                }
                p->setMarkerCount(iMarkerCount);
                p->setMeanMarkerDistance(vMarkerDistances.mean());
            }
            // transfer of ownership from vSegments to m_vSegments
            m_vSegments.add(p);
        }

#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".seg2.tsv");
#endif
        // why "nullAll()" instead of "deleteAll()"?
        // because the pointers were copied to m_vSegments above and
        // AffxArray will delete the pointers in vSegments when it goes out of scope.
        vSegments.nullAll();
        delete[] pLog2Ratios;
    }
    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodMosaicism::run(...) end");
}


// This function produces an initial segmentation.
void CNAnalysisMethodMosaicism::newSegments(unsigned char cChromosome, 
                                            double* pRunMean, 
                                            int iProbeSetCount, 
                                            std::vector<int>& vPositions, 
                                            CNSegmentArray& vSegments)
{
    double dNeutralZoneLow = m_vLossesBoundries[0];
    double dNeutralZoneHigh = m_vGainsBoundries[0];
    int iStartIndex = 0;
    int iEndIndex = 0;
    int iIndex = 0;
    double dCurrentMax = 0;

//    int last_vSegments_size=vSegments.size();
    while (iIndex < (iProbeSetCount - 1))
    {
        iStartIndex = iIndex;

        if (((dNeutralZoneLow < pRunMean[iIndex]) && (pRunMean[iIndex] <= dNeutralZoneHigh)) ||
            (pRunMean[iIndex] != pRunMean[iIndex]))
        {
            // it's a no-call segment
            while ((((dNeutralZoneLow < pRunMean[iIndex]) && (pRunMean[iIndex] <= dNeutralZoneHigh)) ||
                    (pRunMean[iIndex] != pRunMean[iIndex])) && (iIndex < (iProbeSetCount - 1)))
            {
              iIndex++;
            }
            iEndIndex = iIndex;
            if (iIndex < (iProbeSetCount - 1)) {iEndIndex--;}
            if (iEndIndex < iStartIndex) {iEndIndex = iStartIndex;}
            CNSegment* pSegment = new CNSegment();
            pSegment->setCall(0);
            pSegment->setConfidence(0.5);
            pSegment->setSegmentType(getSegmentType());
            pSegment->setChromosome(cChromosome);
            pSegment->setStartIndex(iStartIndex);
            pSegment->setEndIndex(iEndIndex);
            pSegment->setStartPosition(vPositions[iStartIndex]);
            pSegment->setEndPosition(vPositions[iEndIndex]);
            pSegment->setMixture(MCLASS_00);
            pSegment->setCalibratedCN(getCalibratedCN(0, pSegment->getChromosome()));
            vSegments.add(pSegment);
        }
        else if (pRunMean[iIndex] <= dNeutralZoneLow)
        {
            // it's a loss region
            dCurrentMax = pRunMean[iIndex];
            while ((pRunMean[iIndex] < dNeutralZoneLow) && (iIndex < (iProbeSetCount - 1)))
            {
                iIndex++;
                if (pRunMean[iIndex] < dCurrentMax) {dCurrentMax = pRunMean[iIndex];}
            }
            iEndIndex = iIndex;
            if (iIndex < (iProbeSetCount - 1)) {iEndIndex--;}
            if (iEndIndex < iStartIndex) {iEndIndex = iStartIndex;}

            CNSegment* pSegment = new CNSegment();
            pSegment->setCall(1);
            pSegment->setConfidence(0.5);
            pSegment->setSegmentType(getSegmentType());
            pSegment->setChromosome(cChromosome);
            pSegment->setStartIndex(iStartIndex);
            pSegment->setEndIndex(iEndIndex);
            pSegment->setStartPosition(vPositions[iStartIndex]);
            pSegment->setEndPosition(vPositions[iEndIndex]);
            pSegment->setMixture(MosaicClassNeg(getMosaicClass(dCurrentMax, m_vLossesBoundries)));
            pSegment->setCalibratedCN(getCalibratedCN(dCurrentMax, pSegment->getChromosome()));
            double tmp=pSegment->getMixtureAsDouble();
            if (((int)tmp)==tmp) {
              pSegment->setCall(0);
            }
            vSegments.add(pSegment);
        }
        else
        {
            // it's a gain region       
            dCurrentMax = pRunMean[iIndex];
            while ((pRunMean[iIndex] > dNeutralZoneHigh) && (iIndex < (iProbeSetCount - 1)))
            {
                iIndex++;
                if (pRunMean[iIndex] > dCurrentMax) {dCurrentMax = pRunMean[iIndex];}
            }
            iEndIndex = iIndex;
            if (iIndex < (iProbeSetCount - 1)) {iEndIndex--;}
            if (iEndIndex < iStartIndex) {iEndIndex = iStartIndex;}

            CNSegment* pSegment = new CNSegment();
            pSegment->setCall(1);
            pSegment->setConfidence(0.5);
            pSegment->setSegmentType(getSegmentType());
            pSegment->setChromosome(cChromosome);
            pSegment->setStartIndex(iStartIndex);
            pSegment->setEndIndex(iEndIndex);
            pSegment->setStartPosition(vPositions[iStartIndex]);
            pSegment->setEndPosition(vPositions[iEndIndex]);
            pSegment->setMixture(getMosaicClass(dCurrentMax, m_vGainsBoundries));
            pSegment->setCalibratedCN(getCalibratedCN(dCurrentMax, pSegment->getChromosome()));
            //float fMixture = pSegment->getMixture();
            double tmp=pSegment->getMixtureAsDouble();
            if (((int)tmp)==tmp) {
              pSegment->setCall(0);
            }
            vSegments.add(pSegment);
        }

#ifdef MOSAICISM_DEBUG_PSEG
        CNSegment* ts=vSegments[vSegments.size()-1];
        printf("### PSeg: new seg: idx:%d-%d pos:%d-%d mix=%d\n",
               ts->getStartIndex(),ts->getEndIndex(),
               ts->getStartPosition(),ts->getEndPosition(),
               ts->getMixture());
#endif
    }
}


void CNAnalysisMethodMosaicism::cleanSegments(CNSegmentArray& vSegments)
{
#ifdef MOSAICISM_DEBUG
  vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-1.tsv");
#endif

  // 20 is marker index count, not bp
  int iMinSegmentSize = 20;

  // R: set the mixture to MCLASS_00 for small segments.
  //   region.size <- x[,2] - x[,1]
  //   x[region.size < min.size,3] <- 0
  for (int idx1=0;idx1<vSegments.getCount();idx1++) {
        CNSegment* p = vSegments.getAt(idx1);
        if (p->getIndexLen()<iMinSegmentSize) {
            p->setCall(0);
            p->setConfidence(0.5);
            p->setMixture(MCLASS_00);
            p->setCalibratedCN(getCalibratedCN(0, p->getChromosome()));
        }
  }

  //
  CNSegmentArray tmp_segments;
  CNSegment* p1;
  int idx1=0;
  while (idx1<vSegments.getCount()) {
    p1=vSegments.getAt(idx1);
    tmp_segments.push_back(p1);
    vSegments.setAt(idx1,NULL);
    idx1++;

#ifdef MOSAICISM_DEBUG
    printf("%3d: p1->getMixtureAsDouble()=%f\n",idx1,p1->getMixtureAsDouble());
#endif

    // keep non-zeros as is.
    if (p1->getMixture()!=MCLASS_00) {
      continue;
    }

    // a zero, gobble up the other zeros
    while ((idx1<vSegments.getCount()) &&
           (vSegments.getAt(idx1)->getMixture()==MCLASS_00)) {
       p1->setEndIndex(vSegments.getAt(idx1)->getEndIndex());
       delete vSegments.getAt(idx1);
       vSegments.setAt(idx1,NULL);
       idx1++;
    }
  }
  vSegments.nullAll();
  vSegments.clear();
  vSegments=tmp_segments;
  tmp_segments.nullAll();

#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-2.tsv");
#endif

    // Now the hard part, fixing the boundries.
    // First we enlarge or shrink the regions where a CN change is detected.
    if (vSegments.getCount() <= 1)
    {
        return;
    }

    int iSegmentIndex = 0;
    CNSegment* p2;

    for (iSegmentIndex = 0; iSegmentIndex < (vSegments.getCount() - 1); iSegmentIndex++) {
        p1 = vSegments.getAt(iSegmentIndex);
        p2 = vSegments.getAt(iSegmentIndex + 1);
        // first case *p1 has no CN change and *p2 has cn change
        if (p1->getMixture() == MCLASS_00 && p2->getMixture() != MCLASS_00) {
            p2->setStartIndex(Max(1, (p2->getStartIndex() + getBoundaryAdjustment(p2->getMixture()))));
            p1->setEndIndex(p2->getStartIndex() - 1);

        }
        // second case *p1 has CN change and *p2 has no cn change
        else if (p1->getMixture() != MCLASS_00 && p2->getMixture() == MCLASS_00) {
            p1->setEndIndex(Max(0, (p1->getEndIndex() - getBoundaryAdjustment(p1->getMixture()))));
            p2->setStartIndex((p1->getEndIndex() + 1));
        }
        // otherwise p1 and p2 have CN change - don't adjust
    }
#ifdef MOSAICISM_DEBUG
    vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-3.tsv");
#endif

    // Check for segments where mixture = 0 and start > end, and remove them.
    for (int iSegmentIndex = 0; (iSegmentIndex < vSegments.getCount()); iSegmentIndex++)
    {
        CNSegment* p = vSegments.getAt(iSegmentIndex);
        if ((p->getMixture() == MCLASS_00) &&
            (p->getStartIndex() > p->getEndIndex()))
        {
            CNSegment* p2 = vSegments.getAt(iSegmentIndex);
            vSegments.setAt(iSegmentIndex,NULL);
            vSegments.removeAt(iSegmentIndex);
            delete p2;
            iSegmentIndex--;
        }
    }

#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-4.tsv");
#endif

    // Check for segments out of bounds and trim as needed.
    for (int iSegmentIndex = 0; (iSegmentIndex < vSegments.getCount()); iSegmentIndex++)
    {
        CNSegment* p = vSegments.getAt(iSegmentIndex);
        if (p->getStartIndex() < m_iStartIndex) {p->setStartIndex(m_iStartIndex);}
        if (p->getStartIndex() > m_iEndIndex) {p->setStartIndex(m_iEndIndex);}
        if (p->getEndIndex() < m_iStartIndex) {p->setEndIndex(m_iStartIndex);}
        if (p->getEndIndex() > m_iEndIndex) {p->setEndIndex(m_iEndIndex);}
    }

    // Join adjacent segments with same Mixture value. What to do about calibratedCN? Take max of min depending on sign of Mixture.
    // Adjust end points of overlapping segments with different mixture value
    for (int iSegmentIndex1 = 0; (iSegmentIndex1 < (vSegments.getCount() - 1)); iSegmentIndex1++)
    {
        p1 = vSegments.getAt(iSegmentIndex1);
        p2 = vSegments.getAt(iSegmentIndex1 + 1);
        while (p1->getMixture() == p2->getMixture()) {
          // if segments are ajacent then join and go on to next segment
          if (p2->getStartIndex() <= p1->getEndIndex() + 1) {
             // Overlap exists.
             p1->setStartIndex(Min(p1->getStartIndex(), p2->getStartIndex()));
             p1->setEndIndex(Max(p1->getEndIndex(), p2->getEndIndex()));
             if (p1->getMixtureAsDouble() > 0.0)
             {
                 p1->setCalibratedCN(Max(p1->getCalibratedCN(), p2->getCalibratedCN()));
             }
             else
             {
                 p1->setCalibratedCN(Min(p1->getCalibratedCN(), p2->getCalibratedCN()));
             }
             vSegments.setAt(iSegmentIndex1 + 1, NULL);
             vSegments.removeAt(iSegmentIndex1 + 1);
             delete p2;
             if (iSegmentIndex1 < vSegments.getCount() - 1) {
               p2 = vSegments.getAt(iSegmentIndex1 + 1);
             }
             else {
               // last segment
               break;
             }
          }
          else {
            // not adjacent
            break;
          }
        }
    }

#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-5.tsv");
#endif

    // Remove any overlaps by adjusting the start of each segment to be after the end point of the previous segment.
    // If any segments have negative size delete them
    //
    int iLastEnd = 0;
    for (int iSegmentIndex1 = 0; iSegmentIndex1 < vSegments.getCount(); iSegmentIndex1++)
    {
        p1 = vSegments.getAt(iSegmentIndex1);
        // make sure start is after last end point
        if (p1->getStartIndex() < iLastEnd) {
          p1->setStartIndex(iLastEnd);
        }
        // check for end < start
        if (p1->getEndIndex() < p1->getStartIndex()) {
          // delete p1
          delete p1;
          vSegments.removeAt(iSegmentIndex1);
          iSegmentIndex1--;
          continue;
        }
        // set next end point + 1
        iLastEnd = p1->getEndIndex() + 1;
    }
#ifdef MOSAICISM_DEBUG
        vSegments.writeToTsv(m_current_experiment+"-chr"+m_current_chr+".segclean-6.tsv");
#endif
}


// This function determines which class the mosaic falls into
// R-func: get.mos.prop
MosaicClass_t CNAnalysisMethodMosaicism::getMosaicClass(double dMostExtreme, std::vector<double>& vBoundries)
{
    if (fabs(dMostExtreme) >= fabs(vBoundries[3])) {
      return MCLASS_p10;
    }
    if (fabs(dMostExtreme) >= fabs(vBoundries[2])) {
      return MCLASS_p07;
    }
    if (fabs(dMostExtreme) >= fabs(vBoundries[1])) {
      return MCLASS_p05;
    }
    if (fabs(dMostExtreme) >= fabs(vBoundries[0])) {
      return MCLASS_p03;
    }
    return MCLASS_00;
}

float CNAnalysisMethodMosaicism::getCalibratedCN(double dCurrentMax)
{
    dCurrentMax /= m_pEngine->getOptDouble("alpha-cn-calibrate");
    dCurrentMax += m_pEngine->getOptDouble("beta-cn-calibrate");
    return (float)exp(dCurrentMax*log(2.0));
}

float CNAnalysisMethodMosaicism::getCalibratedCN(double dCurrentMax, int chromosome)
{
    double alphaCalibrate;
    double betaCalibrate;
    if (chromosome < m_iXChromosome) {
        alphaCalibrate = m_pEngine->getOptDouble("alpha-cn-calibrate");
        betaCalibrate = m_pEngine->getOptDouble("beta-cn-calibrate");
    }
    else if (chromosome == m_iXChromosome) {
        alphaCalibrate = m_pEngine->getOptDouble("alpha-X-cn-calibrate");
        betaCalibrate = m_pEngine->getOptDouble("beta-X-cn-calibrate");
    }
    else if (chromosome == m_iYChromosome) {
        alphaCalibrate = m_pEngine->getOptDouble("alpha-Y-cn-calibrate");
        betaCalibrate = m_pEngine->getOptDouble("beta-Y-cn-calibrate");
    }
    else {
        alphaCalibrate = 1.0;
        betaCalibrate = 0.0;
    }
    dCurrentMax /= alphaCalibrate;
    dCurrentMax += betaCalibrate;
    return (float)exp(dCurrentMax*log(2.0));
}

// This function determines how to adjust the segment boundaries
// R-func: GetBoundaryAdjustment
//int CNAnalysisMethodMosaicism::getBoundaryAdjustment(double dMosaicClass)
int getBoundaryAdjustment_arr[]=
  { 0,     // Bad
    2165,  // n10
    1777,  // n07
    1054,  // n05
    -1780, // n03
    0,     // _00 - no change
    -2128, // p03
    209,   // p05
    909,   // p07
    1342   // p10
  };
int CNAnalysisMethodMosaicism::getBoundaryAdjustment(MosaicClass_t val)
{
  assert((MCLASS_n10<=val) && (val<=MCLASS_p10));
  return getBoundaryAdjustment_arr[val];
}

void CNAnalysisMethodMosaicism::runmed(double* p, int iCount, int iWindowSize)
{
    std::vector<float>vOut(iCount);
    if ((iWindowSize % 2) == 0) {iWindowSize++;} // Window size should be odd.
    if (iCount <= iWindowSize)
    {
        iWindowSize = iCount - 2;
        if ((iWindowSize % 2) == 0) {iWindowSize--;} // Window size should be odd.
        if (iWindowSize < 0) {return;}
    }
    int iHalfWindow = (int)((double)iWindowSize/2.0);
    for (int iIndex = iHalfWindow; (iIndex < (iCount - iHalfWindow)); iIndex++)
    {
        vOut[iIndex] = percentile(50, (p + iIndex - iHalfWindow), iWindowSize);
    }
    for (int iIndex = (iHalfWindow - 1); (iIndex >= 0); iIndex--)
    {
        vOut[iIndex] = vOut[iIndex + 1];
    }
    for (int iIndex = (iCount - iHalfWindow); (iIndex < iCount); iIndex++)
    {
        vOut[iIndex] = vOut[iIndex - 1];
    }
    for (int iIndex = 0; (iIndex < iCount); iIndex++)
    {
        p[iIndex] = vOut[iIndex];
    }
}

void CNAnalysisMethodMosaicism::calculateSegmentConfidence(double* pRunMean, CNSegment* p)
{
    // iCount should be positive, but was "-312", causing an abort.
    int iCount = p->getEndIndex() - p->getStartIndex() + 1;
    if (iCount<=1) {
#ifdef MOSAICISM_DEBUG
      printf("calculateSegmentConfidence: iCount=%d start=%d end=%d\n",
             iCount,
             p->getStartIndex(),
             p->getEndIndex());
#endif
      return;
    }
    double* pData = new double[iCount];
    int i1 = 0;
    for (int i = p->getStartIndex(); (i <= p->getEndIndex()); i++)
    {
        pData[i1] = pRunMean[i]; i1++;
    }
    runmed(pData, iCount, m_iConfidenceWindow);
    double dSum = 0;
    for (int i = 0; (i < iCount); i++)
    {
        if ((p->getMixtureAsDouble() >= 0.0) && (pData[i] > m_vGainsBoundries[0])) {dSum++;}
        if ((p->getMixtureAsDouble() <= 0.0) && (pData[i] < m_vLossesBoundries[0])) {dSum++;}
    }
    if (iCount == 0) {p->setConfidence(0);}
    else
    {
        double d = fabs(dSum / (double)iCount);
        if (d > 1) {d = 1;}
        p->setConfidence((float)d);
    }
    delete[] pData;
}

void CNAnalysisMethodMosaicism::determineLocalProbeSets()
{
        if(m_bLocalProbeSetsDetermined)
        {
                return;
        }
        int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
        for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
        {
                if( CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsCN() )
                {
                        getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
                }
        }
        m_bLocalProbeSetsDetermined=true;
}

void CNAnalysisMethodMosaicism::writeRunningMeans(const std::string& path,
                                                  int cnt,
                                                  double* pLog2Ratios,
                                                  double* pRunMean)
{
  Verbose::out(1,std::string("CNAnalysisMethodMosaicism::writeRunningMeans: ")+path);
  
  affx::TsvFile tsv;
  tsv.setPrecision(15); // match R digits
  tsv.defineFile("i\tlog2\tmean");
  tsv.addHeader("m_current_chr",m_current_chr);
  tsv.addHeader("m_current_experiment",m_current_experiment);
  tsv.writeTsv(path);
  
  for (int i=0;i<cnt;i++) {
    tsv.set(0,0,i);
    tsv.set(0,1,pLog2Ratios[i]);
    tsv.set(0,2,pRunMean[i]);
    tsv.writeLevel(0);
  }
  tsv.close();
}
