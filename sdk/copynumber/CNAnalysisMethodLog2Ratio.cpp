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
 * @file CNAnalysisMethodLog2Ratio.cpp
 *
 * @brief This file contains the CNAnalysisMethodLog2Ratio class members.
 */
#include "copynumber/CNAnalysisMethodLog2Ratio.h"
//
#include "copynumber/CNAnalysisMethodFactory.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

#define IQR_LIMIT 0.01
/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodLog2Ratio::explainSelf()
{
    CNAnalysisMethodLog2Ratio obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodLog2Ratio::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt opt2 = {"gc-correction", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Log2Ratio GC Correction"};
  opts.push_back(opt2);

  SelfDoc::Opt opt3 = {"median-autosome-median-normalization", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Log2Ratio Median Autosmome Median Normalization"};
  opts.push_back(opt3);

  SelfDoc::Opt opt4 = {"median-smooth-marker-count", SelfDoc::Opt::Integer, "5", "5", "3", "NA", "Log2Ratio Median Smooth marker count."};
  opts.push_back(opt4);

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
SelfCreate* CNAnalysisMethodLog2Ratio::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodLog2Ratio* pMethod = new CNAnalysisMethodLog2Ratio();
    std::string strPrefix = getPrefix();

    pMethod->m_bGCCorrection = setupBoolParameter("gc-correction", strPrefix, params, doc, opts);
    pMethod->m_bMedianAutosomeMedianNormalization = setupBoolParameter("median-autosome-median-normalization", strPrefix, params, doc, opts);
    pMethod->m_iMedianSmoothWindowSize = setupIntParameter("median-smooth-marker-count", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodLog2Ratio::CNAnalysisMethodLog2Ratio()
{
    m_dYTarget = 0;
    m_bGCCorrection = true;
    m_bMedianAutosomeMedianNormalization = true;
    m_iMedianSmoothWindowSize = 0;
    m_dTrimLow = 0;
    m_dTrimHigh = 0;
    m_bLocalProbeSetsDetermined=false;
}

void CNAnalysisMethodLog2Ratio::determineLocalProbeSets()
{
    if(m_bLocalProbeSetsDetermined)
    {
        return;
    }
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
    for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        if( CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAll() )
        {
            if (isCytoScanHD && (!CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsCN())) {continue;}
            getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
        }
    }
    m_bLocalProbeSetsDetermined=true;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodLog2Ratio::run()
{
//    Verbose::out(1, "CNAnalysisMethodLog2Ratio::run(...) start");
    isSetup();
    m_dYTarget = m_pEngine->getOptDouble("yTarget");
//    Verbose::progressBegin(1, "CNAnalysisMethodLog2Ratio::run(...) ", 5, 1, 5);

    determineLocalProbeSets();

    bool cytoScanHD = m_pEngine->getOptBool("cytoscan-hd");

    // if cytoScanHD == true => adjustYTarget == false.  yTarget adjustment is done later by applyYTargetAdjustment().
    calculateLog2Ratios(!cytoScanHD);

	highPassFilterLog2Ratios();

    adjustUsingCovariateLRAdjustment();

    applyWaveCorrection();
//    Verbose::progressStep(1);
    if (m_bGCCorrection)
    {
        GCCorrection(getExperiment());
    }

    if (cytoScanHD) {applyYTargetAdjustment();}

//    Verbose::progressStep(1);
    if (m_bMedianAutosomeMedianNormalization)
    {
        calculateMedianAutosomeMedian(getExperiment());
//        Verbose::progressStep(1);
        normalizeLog2Ratios(getExperiment());
    }
//    Verbose::progressStep(1);
    calculateQCMetrics(getExperiment());
    getProbeSets()->calculateLog2RatioMedianSmooth(m_iMedianSmoothWindowSize);
//    Verbose::progressEnd(1, "Done");
//    Verbose::out(1, "CNAnalysisMethodLog2Ratio::run(...) end");
}

/**
 * Calculate the log2 ratios.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * Note that the signal is A + B.
 * @param int - The col index.
 * @param int - The experiment index.
 */
void CNAnalysisMethodLog2Ratio::calculateLog2Ratios(bool adjustYTarget)
{
//    Verbose::out(3, "CNAnalysisMethodLog2Ratio::calculateLog2Ratios");
    double yTarget = ((adjustYTarget == true) ? m_dYTarget : 0.0);
    
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        double dSignal = (pobjProbeSet->getAAlleleSignal() + pobjProbeSet->getBAlleleSignal());
        double dLog2Ratio = 0;
        if (pobjProbeSet->getChromosome() == m_iXChromosome)
        {
            dLog2Ratio = log2(dSignal / pobjProbeSet->getXXMedianSignal());
        }
        else if (pobjProbeSet->getChromosome() == m_iYChromosome)
        {
            dLog2Ratio = log2(dSignal / pobjProbeSet->getYMedianSignal()) + yTarget;
        }
        else // Autosome
        {
            dLog2Ratio = log2(dSignal / pobjProbeSet->getMedianSignal());
        }
        pobjProbeSet->setLog2Ratio((float)dLog2Ratio);
    }
}

void CNAnalysisMethodLog2Ratio::applyYTargetAdjustment()
{
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        if (pobjProbeSet->getChromosome() == m_iYChromosome)
        {
            double dLog2Ratio = pobjProbeSet->getLog2Ratio() + m_dYTarget;
            pobjProbeSet->setLog2Ratio((float)dLog2Ratio);
        }
    }
}

/**
 * Apply the GC content correction. This is used to reduce wavyness in the smoothed log2 ratios.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param int - The col index.
 * @param int - The experiment index.
 */
void CNAnalysisMethodLog2Ratio::GCCorrection(CNExperiment* pobjExperiment)
{
//    Verbose::out(3, "CNLog2RatioData::GCCorrection");


    // This functionality was moved here from the original location in the Annotations file read.
    // Since we do GC correction on a subset of the probesets we allow the module to determine which
    // probe sets to use via the local choice in getProbeSets.
    getProbeSets()->setupGCCorrectionBins(m_pEngine->getOptInt("gc-correction-bin-count"));


    int iEnzymeCount = 0;
    AffxMultiDimensionalArray<float>& vLog2RatiosBeforeCorrection = getV2();
    std::vector<int> vBinIndexes(getProbeSets()->getCount());
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        if ((pobjProbeSet->isSty()) || (pobjProbeSet->isNsp())) {iEnzymeCount++;}
        vLog2RatiosBeforeCorrection.set(iRowIndex, pobjProbeSet->getLog2Ratio());
        vBinIndexes[iRowIndex] = pobjProbeSet->getGCBinIndex();
    }
    // calculateIqr uses V1
    double dOldIqr = calculateIqr();

    if (iEnzymeCount == 0)
    {
        processMarkerGroupGCCorrection(true, false, false, vBinIndexes);
        processMarkerGroupGCCorrection(false, false, false, vBinIndexes);
    }
    else
    {
        processMarkerGroupGCCorrection(true, true, false, vBinIndexes);
        processMarkerGroupGCCorrection(true, false, true, vBinIndexes);
        processMarkerGroupGCCorrection(true, true, true, vBinIndexes);
        processMarkerGroupGCCorrection(true, false, false, vBinIndexes);
        processMarkerGroupGCCorrection(false, true, false, vBinIndexes);
        processMarkerGroupGCCorrection(false, false, true, vBinIndexes);
        processMarkerGroupGCCorrection(false, true, true, vBinIndexes);
        processMarkerGroupGCCorrection(false, false, false, vBinIndexes);
    }

    double dNewIqr = calculateIqr();
    float fLog2Ratio = 0;
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        fLog2Ratio = pobjProbeSet->getLog2Ratio();
        if (dNewIqr != 0) {fLog2Ratio *= (dOldIqr / dNewIqr);}
        pobjProbeSet->setLog2Ratio(fLog2Ratio);
    }
    AffxMultiDimensionalArray<float>& vLog2RatioCorrectionDifferences = getV3();
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        fLog2Ratio = pobjProbeSet->getLog2Ratio();
        vLog2RatioCorrectionDifferences.set(iRowIndex, fabs(vLog2RatiosBeforeCorrection.get(iRowIndex) - fLog2Ratio));
    }
    pobjExperiment->setGCCorrectionMetric((float)vLog2RatioCorrectionDifferences.median());
}

/**
 * @brief Process the marker group for GC correction
 * @param bool - Is the marker group SNPs
 * @param bool - Is the marker group Sty
 * @param bool - Is the marker group Nsp
 */
void CNAnalysisMethodLog2Ratio::processMarkerGroupGCCorrection(bool bSnp, bool bSty, bool bNsp, std::vector<int>& vBinIndexes)
{
    AffxMultiDimensionalArray<float>& vLog2Ratios = getV4();
    float fLog2Ratio = 0;
    int iMarkerCount = 0;
    double dPrevMedian = 0;
    std::vector<int> vProbeSetIndexes(getProbeSets()->getCount());
    for (int iBinIndex = 0; (iBinIndex < m_pEngine->getOptInt("gc-correction-bin-count")); iBinIndex++)
    {
        int iIndex = 0;
        iMarkerCount = 0;
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            if (vBinIndexes[iRowIndex] == iBinIndex)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
                if (isInGCCorrectionGroup(*pobjProbeSet, (m_iXChromosome - 1), bSnp, bSty, bNsp, iBinIndex))
                {
                    fLog2Ratio = pobjProbeSet->getLog2Ratio();
                    vLog2Ratios.set(iMarkerCount, fLog2Ratio);
                    vProbeSetIndexes[iMarkerCount] = iRowIndex;
                    iMarkerCount++;
                }
                if (isInGCCorrectionGroup(*pobjProbeSet, m_iYChromosome, bSnp, bSty, bNsp, iBinIndex))  // apply to autosomes and sex chroms.
                {
                    vProbeSetIndexes[iIndex] = iRowIndex;
                    iIndex++;
                }
            }
        }
        if (iIndex < getProbeSets()->getCount()) {vProbeSetIndexes[iIndex] = -1;}
        double dMedian = vLog2Ratios.median(iMarkerCount);
        if (dMedian != dMedian) {dMedian = dPrevMedian;} // Check for NaN
        if (dMedian != dMedian) {dMedian = 0;}
        dPrevMedian = dMedian;

        iMarkerCount = 0;
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            if (vProbeSetIndexes[iRowIndex] == -1) {break;}
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(vProbeSetIndexes[iRowIndex]);
            pobjProbeSet->setGcAdjustment((float)dMedian);
            fLog2Ratio = pobjProbeSet->getLog2Ratio();
            fLog2Ratio -= dMedian;
            pobjProbeSet->setLog2Ratio(fLog2Ratio);
            vLog2Ratios.set(iMarkerCount, fLog2Ratio);
            iMarkerCount++;
        }
        // Took out bug where if medain == 0 the continue.
        double dIqr = (vLog2Ratios.percentile(0.75, iMarkerCount) - vLog2Ratios.percentile(0.25, iMarkerCount));
        if (((-IQR_LIMIT <= dIqr) && (dIqr <= IQR_LIMIT)) || (dIqr != dIqr)) {continue;}
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            if (vProbeSetIndexes[iRowIndex] == -1) {break;}
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(vProbeSetIndexes[iRowIndex]);
            fLog2Ratio = pobjProbeSet->getLog2Ratio();
            if (dIqr != 0) {fLog2Ratio /= dIqr;}
            pobjProbeSet->setLog2Ratio(fLog2Ratio);
        }
    }
}

/**
 * @brief Is the probe set belong to the marker group
 * @param CNProbeSet& - The specified probe set
 * @param int - The last autosome chromosome
 * @param bool - Is the marker group SNPs
 * @param bool - Is the marker group Sty
 * @param bool - Is the marker group Nsp
 * @param int - The marker group GC content bin index
 */
bool CNAnalysisMethodLog2Ratio::isInGCCorrectionGroup(CNProbeSet& objProbeSet, int iLastAutosomeChromosome, bool bSnp, bool bSty, bool bNsp, int iBinIndex)
{
    return ((objProbeSet.getChromosome() <= iLastAutosomeChromosome) &&
        (objProbeSet.processAsSNP() == bSnp) && ((objProbeSet.isSty() == bSty) && (objProbeSet.isNsp() == bNsp)) &&
        (objProbeSet.getGCBinIndex() == iBinIndex));
}

/**
 * Calcuate the median autosome median. This is used to normalize the log2 ratios.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param CNExperiment* - The experiment to process
 */
void CNAnalysisMethodLog2Ratio::calculateMedianAutosomeMedian(CNExperiment* pobjExperiment)
{
//    Verbose::out(3, "CNLog2RatioData::calculateMedianAutosomeMedian in");
    int iLastAutosomeChromosome = getLastAutosomeChromosome();
    AffxMultiDimensionalArray<int> vChrCounts(iLastAutosomeChromosome);
    AffxMultiDimensionalArray<int> vChrIndexes(iLastAutosomeChromosome);
    AffxMultiDimensionalArray<float> vChrMedians(iLastAutosomeChromosome);
    AffxMultiDimensionalArray<float>* pvar = new AffxMultiDimensionalArray<float>[iLastAutosomeChromosome];
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        if (pobjProbeSet->getChromosome() <= 0) {continue;}
        if (pobjProbeSet->getChromosome() <= iLastAutosomeChromosome)
        {
            vChrCounts.increment(pobjProbeSet->getChromosome() - 1);
        }
    }
    for (int iIndex = 0; (iIndex < iLastAutosomeChromosome); iIndex++)
    {
        pvar[iIndex].initialize(vChrCounts.get(iIndex));
    }
    vChrIndexes.initialize();
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        if (pobjProbeSet->getChromosome() <= 0) {continue;}
        if (pobjProbeSet->getChromosome() <= iLastAutosomeChromosome)
        {
            pvar[pobjProbeSet->getChromosome() - 1].set(vChrIndexes.get(pobjProbeSet->getChromosome() - 1), pobjProbeSet->getLog2Ratio());
            vChrIndexes.increment(pobjProbeSet->getChromosome() - 1);
        }
    }
    for (int iIndex = 0; (iIndex < iLastAutosomeChromosome); iIndex++)
    {
        vChrMedians.set(iIndex, (float)pvar[iIndex].median());
    }
    pobjExperiment->setMedianAutosomeMedian((float)vChrMedians.finiteMedian());
    delete[] pvar;
//    Verbose::out(3, "CNLog2RatioData::calculateMedianAutosomeMedian out");
}

/**
 * Normalize the log2 ratios by subtracting out the median autosome median.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param CNExperiment* - The experiment to process
 */
void CNAnalysisMethodLog2Ratio::normalizeLog2Ratios(CNExperiment* pobjExperiment)
{
//    Verbose::out(3, "CNLog2RatioData::normalizeLog2Ratios");
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        pobjProbeSet->setLog2Ratio((float)(pobjProbeSet->getLog2Ratio() - pobjExperiment->getMedianAutosomeMedian()));
    }
}

/**
 * Trim the log2 ratios.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param CNExperiment* - The experiment to process
 */
void CNAnalysisMethodLog2Ratio::trim()
{
//    Verbose::out(3, "CNLog2RatioData::normalizeLog2Ratios");
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        float fLog2Ratio = pobjProbeSet->getLog2Ratio();
        if (fLog2Ratio > m_dTrimHigh) {fLog2Ratio = m_dTrimHigh;}
        else if (fLog2Ratio < m_dTrimLow) {fLog2Ratio = m_dTrimLow;}
        pobjProbeSet->setLog2Ratio(fLog2Ratio);
    }
}

/**
 * Calcuate the only the iqr QC Metric.
 */
double CNAnalysisMethodLog2Ratio::calculateIqr()
{
    AffxMultiDimensionalArray<float>& v = getV1();
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        v.set(iRowIndex, pobjProbeSet->getLog2Ratio());
    }
    return (v.percentile(0.75, getProbeSets()->getCount()) - v.percentile(0.25, getProbeSets()->getCount()));
}

/**
 * Calcuate the QC Metrics.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param CNExperiment* - The experiment to process
 */
void CNAnalysisMethodLog2Ratio::calculateQCMetrics(CNExperiment* pobjExperiment)
{
//    Verbose::out(3, "CNLog2RatioData::calculateQCMetrics");
    AffxMultiDimensionalArray<float>& vAbsDiff = getV1();
    AffxMultiDimensionalArray<float>& v = getV2();
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        v.set(iRowIndex, pobjProbeSet->getLog2Ratio());
    }
    // m_arProbeSetsAlternateSort order is Chromosome, Position
    int iCount = 0;
    for (int iRowIndex = 1; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjPrevProbeSet = getProbeSets()->getAt(iRowIndex - 1);
        CNProbeSet* pobjNextProbeSet = getProbeSets()->getAt(iRowIndex);
        if (pobjPrevProbeSet->getChromosome() == pobjNextProbeSet->getChromosome())
        {
            vAbsDiff.set(iCount, fabs(pobjNextProbeSet->getLog2Ratio() - pobjPrevProbeSet->getLog2Ratio()));
            iCount++;
        }
    }
    pobjExperiment->setMadDiffCN((float)vAbsDiff.median(iCount));
    pobjExperiment->setIqr((float)(v.percentile(0.75, getProbeSets()->getCount()) - v.percentile(0.25, getProbeSets()->getCount())));

    AffxMultiDimensionalArray<float>& vAbs = getV1();
    iCount = 0;
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        float fLog2Ratio = pobjProbeSet->getLog2Ratio();
        vAbs.set(iCount, fabs(fLog2Ratio));
        iCount++;
    }
    pobjExperiment->setMeanAbsRle(vAbs.mean(iCount));
}

void CNAnalysisMethodLog2Ratio::applyWaveCorrection()
{
    if (!m_pEngine->isOptDefined("wave-correction-log2ratio-adjustment-method")) {return;}
    if (m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") == "none") {return;}
    if (m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") == "") {return;}

    CNAnalysisMethodFactory amFactory;
    CNAnalysisMethod* am = NULL;
    am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method"));
    am->setEngine(m_pEngine);
    am->setup(*getExperiment(), *CNAnalysisMethod::getProbeSets());
    am->run();
    delete am;
}

void CNAnalysisMethodLog2Ratio::adjustUsingCovariateLRAdjustment()
{
    if (m_pEngine->getEngineName() != "CNCytoEngine")
    {
        return;
    }
    CNAnalysisMethodFactory amFactory;
    CNAnalysisMethod* am = NULL;
    am = amFactory.CNAnalysisMethodForString("covariate-lr-adjuster");
    am->setEngine(m_pEngine);
    am->setup(*getExperiment(), *getProbeSets());
    am->run();
    delete am;
}

void CNAnalysisMethodLog2Ratio::highPassFilterLog2Ratios()
{
    if (!m_pEngine->isOptDefined("log2ratio-adjustment-method")) {return;}
    if (m_pEngine->getOpt("log2ratio-adjustment-method") == "none") {return;}
    if (m_pEngine->getOpt("log2ratio-adjustment-method") == "") {return;}

    CNAnalysisMethodFactory amFactory;
    CNAnalysisMethod* am = NULL;
    am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("log2ratio-adjustment-method"));
    am->setEngine(m_pEngine);
    am->setup(*getExperiment(),*CNAnalysisMethod::getProbeSets(),getProbes());
    am->run();

    delete am;
}

