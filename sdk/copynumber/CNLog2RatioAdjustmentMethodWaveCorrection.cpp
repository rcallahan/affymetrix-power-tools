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
 * @file CNLog2RatioAdjustmentMethodWaveCorrection.cpp
 *
 * @brief This file contains the CNLog2RatioAdjustmentMethodWaveCorrection class members.
 */
#include "copynumber/CNLog2RatioAdjustmentMethodWaveCorrection.h"
//
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
//
#include "util/Err.h"
#include "util/ErrHandler.h"

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNLog2RatioAdjustmentMethodWaveCorrection::explainSelf()
{
    CNLog2RatioAdjustmentMethodWaveCorrection obj;
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
std::vector<SelfDoc::Opt> CNLog2RatioAdjustmentMethodWaveCorrection::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt Bandwidth = {"bandwidth", SelfDoc::Opt::Integer, "101", "101", "NA", "NA", "Wave Correction bandwidth."};
  opts.push_back(Bandwidth);

  SelfDoc::Opt BinCount = {"bin-count", SelfDoc::Opt::Integer, "25", "25", "NA", "NA", "Wave Correction bin count."};
  opts.push_back(BinCount);

  SelfDoc::Opt WaveCount = {"wave-count", SelfDoc::Opt::Integer, "-1", "-1", "-1", "NA", "Wave Correction count."};
  opts.push_back(WaveCount);

  SelfDoc::Opt WaveSmooth = {"wave-smooth", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                             "Turn the non-parametric smoothing on or off."};
  opts.push_back(WaveSmooth);

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
SelfCreate* CNLog2RatioAdjustmentMethodWaveCorrection::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNLog2RatioAdjustmentMethodWaveCorrection* pMethod = new CNLog2RatioAdjustmentMethodWaveCorrection();
    std::string strPrefix = getPrefix();

    pMethod->m_iBandwidth = setupIntParameter("bandwidth", strPrefix, params, doc, opts);
    pMethod->m_iBinCount = setupIntParameter("bin-count", strPrefix, params, doc, opts);
    pMethod->m_iWaveCount = setupIntParameter("wave-count", strPrefix, params, doc, opts);
    pMethod->m_bWaveSmooth = setupBoolParameter("wave-smooth", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNLog2RatioAdjustmentMethodWaveCorrection::CNLog2RatioAdjustmentMethodWaveCorrection()
{
    m_iBandwidth = 0;
    m_iBinCount = 0;
    m_iWaveCount = -1;
}

/**
 * @brief Run the analysis.
 */
void CNLog2RatioAdjustmentMethodWaveCorrection::run()
{
    Verbose::out(1, "CNLog2RatioAdjustmentMethodWaveCorrection::run(...) start");
    isSetup();

    determineLocalProbeSets();
    if (m_pEngine->getOptBool("keep-intermediate-data"))
    {
        writeLog2Ratios("beforeWaveCorrection");
    }

    // The number of waves that we use will be set in the variable iNumberOfWavesToUse.  If no wave-count value is set on the
    // command line, at this point in the code we have m_iWaveCount equal to -1.  The plan is then to default
    // to the number of waves found in the CNReference input file.  If the wave-count value has been set on the
    // the command line we use this value only if it does not exceed the number of waves defined in the CNReference
    // file. If it does exceed the number of waves in the reference file, we default to the number found in the
    // reference and issue a warning.

    // One more complexity is that any instantiation of amFactory between samples will reset the value of m_iWaveCount
    // to whatever is on the command line, or the default of -1. Thus we need to calculate the value of iNumberOfWaveToUse
    // for each sample.  (This will add an open, read of number of columns  and close of the reference file.) Another option
    // which was previously implemented was to find the number of waves on probeset object.  This seemed to rely too much on
    // the algorithm not modifying the probesets in a particular way.  I didn't want to rely on that.

    AffxString strReferenceFileName = m_pEngine->getOpt("reference-input");
    if (strReferenceFileName == "") {strReferenceFileName = m_pEngine->getOpt("reference-output");}
    int iProbeSetCount = 0;
    if (m_pEngine->getOptBool("cyto2")) {iProbeSetCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "WaveCorrection");}
    if (iProbeSetCount <= 0)
    {
        m_strGroup = "CopyNumber";
        iProbeSetCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CopyNumber", "WaveCorrection");
        if (iProbeSetCount <= 0)
        {
            Verbose::out(1, "WARNING: Cannot apply Wave Correction as the supporting data does not exist in the reference.");
            return;
        }
    } else {m_strGroup = "Cyto2";}

    // Here we find the number of waves in the reference file.
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
    group5 = file5.openGroup(m_strGroup, affx::FILE5_OPEN_RO);
    double d = 0;
    tsv5 = group5->openTsv("WaveCorrection", affx::FILE5_OPEN_RO);
    int iIndex = 0;
    int iColumnCount = tsv5->getColumnCount(0);

    // We report the number of waves found in the reference file in the "affymetrix-algorithm-state" section of the cychp file.
    if (getExperiment()->getIndex() == 0)
    {
        getEngine()->setOpt("reference-wave-count", ToStr(iColumnCount));
    }
	
    int iNumberOfWavesToUse;
    if (m_iWaveCount <= -1)
    {
        Verbose::out(2, "NOTE: The number of waves used in the wave-correction-log2ratio-adjustment-method will be the number of waves found in the CNReference file. This is because a separate wave-count parameter was not specified in command line input.");
        iNumberOfWavesToUse =  iColumnCount;
    } else
    {
        if (m_iWaveCount > iColumnCount)
        {
            Err::errAbort("Analysis wave correction wave-count parameter value " + ToStr(m_iWaveCount) +
                " exceeds wave-count " + ToStr(iColumnCount) +
                " available in reference file " + strReferenceFileName +
                ". The number of waves used will be the number found in the reference file.");
        }
        iNumberOfWavesToUse = Min(m_iWaveCount, iColumnCount);
    }
    // We report the number of waves actually used in the "affymetrix-chipsummary-" section of the cychp file.
    getExperiment()->setNumberOfWavesUsed(iNumberOfWavesToUse);

    // Here we load the waves information into the probesets.  This is only done for the first experiment.
    if (getExperiment()->getIndex() == 0)
    {
        CNAnalysisMethod::getProbeSets()->quickSort(0);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
			if (iIndex >= CNAnalysisMethod::getProbeSets()->size()) {Err::errAbort("CN Ref file and CNAnalysisMethod::ProbeSets vector are out of sync.\t" + ToStr(CNAnalysisMethod::getProbeSets()->size()));}
            CNProbeSet* p = CNAnalysisMethod::getProbeSets()->getAt(iIndex);
            p->getWaves().resize(iNumberOfWavesToUse);
            for (int iColumnIndex = 0; (iColumnIndex < iNumberOfWavesToUse); iColumnIndex++)
            {
                tsv5->get(0, iColumnIndex, &d); p->getWaves()[iColumnIndex] = (float)d;
            }
            iIndex++;
        }
        CNAnalysisMethod::getProbeSets()->quickSort(1);
    }
    tsv5->close();
    delete tsv5;
    group5->close();
    delete group5;
    file5.close();

    Verbose::out(2,"The number of waves used for wave-correction is: " + ToStr(iNumberOfWavesToUse));

    int iLastAutosomeChromosome = getLastAutosomeChromosome();
    int iAutosomeCount = 0;
    int iLocalProbeSetCount=getProbeSets()->getCount();
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iLocalProbeSetCount); iProbeSetIndex++)
    {
        CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
        if (p->getChromosome() <= iLastAutosomeChromosome) {iAutosomeCount++;}
    }
    if (m_bMedianAutosomeMedianNormalization)
    {
        calculateMedianAutosomeMedian(getExperiment());
        normalizeLog2Ratios(getExperiment());
    }

    for (int iWaveIndex = 0; iWaveIndex < iNumberOfWavesToUse; iWaveIndex++)
    {
        adjustLog2Ratios(iAutosomeCount, iLastAutosomeChromosome, iWaveIndex);
        if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
            writeLog2Ratios( "afterWaveSubtract" + getInt(iWaveIndex+1) );
        }
    }

    if (m_bWaveSmooth) {
        for (int iWaveIndex = 0; iWaveIndex < iNumberOfWavesToUse; iWaveIndex++)
        {
            applyRunningMean(iWaveIndex);
            if (m_pEngine->getOptBool("keep-intermediate-data"))
            {
                writeLog2Ratios( "afterWaveSmooth" + getInt(iWaveIndex+1) );
            }
        }
    }

    Verbose::out(1, "CNLog2RatioAdjustmentMethodWaveCorrection::run(...) end");
}

void CNLog2RatioAdjustmentMethodWaveCorrection::adjustLog2Ratios(       int iAutosomeCount,
                                                                        int iLastAutosomeChromosome,
                                                                        int iWaveIndex)
{
    int iProbeSetCount = getProbeSets()->getCount();
    Matrix mxTemp(iAutosomeCount, 1);
    int iTempIndex = 0;
    double dSum = 0;
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
        float fValue = p->getWaveValue(iWaveIndex);
        if (p->getChromosome() <= iLastAutosomeChromosome)
        {
            dSum += (fValue * p->getLog2Ratio());
            mxTemp.element(iTempIndex, 0) = fValue; iTempIndex++;
        }
    }
    double dTempNorm = norm(mxTemp);
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
        p->setLog2Ratio(p->getLog2Ratio() - dSum * p->getWaveValue(iWaveIndex) / dTempNorm);
    }
    if (m_bMedianAutosomeMedianNormalization)
    {
        calculateMedianAutosomeMedian(getExperiment());
        normalizeLog2Ratios(getExperiment());
    }
    if (m_pEngine->isOptDefined("wave-correction-log2ratio-adjustment-method")     &&
        m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != "none" &&
        m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != ""     &&
        getExperiment())
    {
        getExperiment()->addWavinessAmplitude(iWaveIndex, dSum/dTempNorm);
    }
}

void CNLog2RatioAdjustmentMethodWaveCorrection::applyRunningMean(int iWaveIndex)
{
    std::vector<float> vWave(getProbeSets()->getCount());
    vector<int> vChromosomes = getChromosomes(getProbeSets());
    for (int i = 0; (i < vChromosomes.size()); i++)
    {
        int iChromosome = vChromosomes[i];
        if (iChromosome == 255)
        {
                continue;
        }
        int iProbeSetCount = getProbeSetCount(iChromosome, getProbeSets());
        if (iProbeSetCount == 0) continue; // probably 23 or Y
        int iStartChromosome = getChrBounds(iChromosome, getProbeSets()).first;
        double* p = new double[iProbeSetCount];
        double* pRunMean = new double[iProbeSetCount];
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            p[iProbeSetIndex] = getProbeSets()->getAt(iStartChromosome + iProbeSetIndex)->getWaveValue(iWaveIndex);
        }
        runmean(p, pRunMean, &iProbeSetCount, &m_iBandwidth);
        delete[] p;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            vWave[iStartChromosome + iProbeSetIndex] = (float)pRunMean[iProbeSetIndex];
        }
        delete[] pRunMean;
    }
    int iProbeSetCount = getProbeSets()->getCount();

    // Setup bins
    std::vector<int> vBinIndexes(iProbeSetCount);
    bin(vWave, vBinIndexes, m_iBinCount);

    // Apply adjustment.
    for (int iProbeSetType = 1; (iProbeSetType < 3); iProbeSetType++)
    {
        for (int iBinIndex = 0; (iBinIndex < m_iBinCount); iBinIndex++)
        {
            // Get autosome count per bin.
            int iProbeSetCount = 0;
            int iAutosomeProbeSetCount = 0;
            for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
                if ((iProbeSetType == 1) && (pobjProbeSet->getProcessFlag() != 1)) {continue;}
                if ((iProbeSetType == 2) && (pobjProbeSet->getProcessFlag() < 2)) {continue;}
                if (vBinIndexes[iRowIndex] == iBinIndex)
                {
                    iProbeSetCount++;
                    if (pobjProbeSet->getChromosome() < m_iXChromosome) // is autosome.
                    {
                        iAutosomeProbeSetCount++;
                    }
                }
            }
            std::vector<int> vIndexes(iProbeSetCount);
            std::vector<int> vAutosomeIndexes(iAutosomeProbeSetCount);
            // Get autosome count per bin.
            int i1 = 0;
            int i2 = 0;
            for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
                if ((iProbeSetType == 1) && (pobjProbeSet->getProcessFlag() != 1)) {continue;}
                if ((iProbeSetType == 2) && (pobjProbeSet->getProcessFlag() < 2)) {continue;}
                if (vBinIndexes[iRowIndex] == iBinIndex)
                {
                    vIndexes[i1] = iRowIndex; i1++;
                    if (pobjProbeSet->getChromosome() < m_iXChromosome) // is autosome.
                    {
                        vAutosomeIndexes[i2] = iRowIndex; i2++;
                    }
                }
            }
            if (iAutosomeProbeSetCount > 0)
            {
                // Calculate autosome median for each bin.
                AffxMultiDimensionalArray<float> vLog2Ratios(iAutosomeProbeSetCount);
                for (int iRowIndex = 0; (iRowIndex < vAutosomeIndexes.size()); iRowIndex++)
                {
                    CNProbeSet* pobjProbeSet = getProbeSets()->getAt(vAutosomeIndexes[iRowIndex]);
                    vLog2Ratios.set(iRowIndex, pobjProbeSet->getLog2Ratio());
                }
                // Perform actual adjustment.
                double dMedianLog2Ratio = vLog2Ratios.median();
                if (dMedianLog2Ratio == dMedianLog2Ratio) // Not NaN
                {
                    for (int iRowIndex = 0; (iRowIndex < vIndexes.size()); iRowIndex++)
                    {
                        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(vIndexes[iRowIndex]);
                        float fLog2Ratio = pobjProbeSet->getLog2Ratio();
                        fLog2Ratio -= dMedianLog2Ratio;
                        pobjProbeSet->setLog2Ratio(fLog2Ratio);
                    }
                }
            }
        }
    }
}


