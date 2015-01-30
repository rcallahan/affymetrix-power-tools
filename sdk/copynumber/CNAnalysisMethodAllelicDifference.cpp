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
 * @file CNAnalysisMethodAllelicDifference.cpp
 *
 * @brief This file contains the CNAnalysisMethodAllelicDifference class members.
 */
#include "copynumber/CNAnalysisMethodAllelicDifference.h"
//


/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodAllelicDifference::explainSelf()
{
    CNAnalysisMethodAllelicDifference obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodAllelicDifference::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt opt1 = {"outlier-trim", SelfDoc::Opt::Double, "3.0", "3.0", "NA", "NA", "AllelicDifference Outlier Trim"};
  opts.push_back(opt1);

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
SelfCreate* CNAnalysisMethodAllelicDifference::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodAllelicDifference* pMethod = new CNAnalysisMethodAllelicDifference();
    std::string strPrefix = getPrefix();

    pMethod->m_dAllelicDifferenceOutlierTrim = setupDoubleParameter("outlier-trim", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodAllelicDifference::CNAnalysisMethodAllelicDifference()
{
    m_iXChromosome = 0;
    m_iYChromosome = 0;
    m_dAllelicDifferenceOutlierTrim = 0;
    m_bLocalProbeSetsDetermined = false;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodAllelicDifference::run()
{
    Verbose::out(1, "CNAnalysisMethodAllelicDifference::run(...) start");
    isSetup();
	if (!m_pEngine->getOptBool("cyto2")) {determineLocalProbeSets();}
//    Verbose::progressBegin(1, "CNAnalysisMethodAllelicDifference::run(...) ", 1, 1, 1);
//    Verbose::progressStep(1);
    calculateAllelicDifferences();
//    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodAllelicDifference::run(...) end");
}

/**
 * Calcuate the Allelic Differences.
 * The data is processed by experiment meaning the all the probe set data is processed for a given experiment.
 * @param int - The col index.
 * @param int - The experiment index.
 */
void CNAnalysisMethodAllelicDifference::calculateAllelicDifferences()
{
//    Verbose::out(3, "CNAllelicDifferenceData::calculateAllelicDifferences");
    AffxMultiDimensionalArray<float>& vAAMedianSignal = getV1();
    AffxMultiDimensionalArray<float>& vABMedianSignal = getV2();
    AffxMultiDimensionalArray<float>& vBBMedianSignal = getV3();
    AffxMultiDimensionalArray<float>& vAAMinusABMedianSignal = getV4();
    AffxMultiDimensionalArray<float>& vBBMinusABMedianSignal = getV5();
    int iAAandABCount = 0;
    int iBBandABCount = 0;
    double NaN = numeric_limits<double>::quiet_NaN();
    for (int iRowIndex = 0; (iRowIndex < CNAnalysisMethod::getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = CNAnalysisMethod::getProbeSets()->getAt(iRowIndex);
        vAAMedianSignal.set(iRowIndex, pobjProbeSet->getAAMedianSignal());
        vABMedianSignal.set(iRowIndex, pobjProbeSet->getABMedianSignal());
        vBBMedianSignal.set(iRowIndex, pobjProbeSet->getBBMedianSignal());
        // Not AA, but AB and BB, so we need to set the AA value.
        if ((isNaN(vAAMedianSignal.get(iRowIndex))) && (!isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vBBMedianSignal.get(iRowIndex))))
        {
            vAAMedianSignal.set(iRowIndex, (vABMedianSignal.get(iRowIndex) - vBBMedianSignal.get(iRowIndex)));
        }
        // Not BB, but AA and AB, so we need to set the BB value.
        if ((isNaN(vBBMedianSignal.get(iRowIndex))) && (!isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vAAMedianSignal.get(iRowIndex))))
        {
            vBBMedianSignal.set(iRowIndex, (vABMedianSignal.get(iRowIndex) - vAAMedianSignal.get(iRowIndex)));
        }
        // Not AB, but AA and BB, so we need to set the AB value.
        if ((isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vAAMedianSignal.get(iRowIndex))) && (!isNaN(vBBMedianSignal.get(iRowIndex))))
        {
            vABMedianSignal.set(iRowIndex, (float)(0.5 * (vAAMedianSignal.get(iRowIndex) + vBBMedianSignal.get(iRowIndex))));
        }
        // Not AA and not AB, but BB, so we need to set the AA value and the AB value.
        if ((isNaN(vAAMedianSignal.get(iRowIndex))) && (isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vBBMedianSignal.get(iRowIndex))))
        {
            vABMedianSignal.set(iRowIndex, 0);
            vAAMedianSignal.set(iRowIndex, -vBBMedianSignal.get(iRowIndex));
        }
        // Not BB and not AB, but AA, so we need to set the AB and BB values.
        if ((isNaN(vBBMedianSignal.get(iRowIndex))) && (isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vAAMedianSignal.get(iRowIndex))))
        {
            vABMedianSignal.set(iRowIndex, 0);
            vBBMedianSignal.set(iRowIndex, -vAAMedianSignal.get(iRowIndex));
        }
        //AA and AB. Get the AA - AB median signals so we can calculate the median of them.
        if ((!isNaN(vAAMedianSignal.get(iRowIndex))) && (!isNaN(vABMedianSignal.get(iRowIndex))))
        {
            vAAMinusABMedianSignal.set(iAAandABCount, (vAAMedianSignal.get(iRowIndex) - vABMedianSignal.get(iRowIndex)));
            iAAandABCount++;
        }
        // BB and AB. Get the BB - AB median signals so we can calculate the median of them.
        if ((!isNaN(vBBMedianSignal.get(iRowIndex))) && (!isNaN(vABMedianSignal.get(iRowIndex))))
        {
            vBBMinusABMedianSignal.set(iBBandABCount, (vBBMedianSignal.get(iRowIndex) - vABMedianSignal.get(iRowIndex)));
            iBBandABCount++;
        }
    }
    // Calculate the median of the median singal, and use this to fill in where we have AB but not AA and BB.
    double dAAMinusABMedian = vAAMinusABMedianSignal.median(iAAandABCount);
    double dBBMinusABMedian = vBBMinusABMedianSignal.median(iBBandABCount);
    for (int iRowIndex = 0; (iRowIndex < CNAnalysisMethod::getProbeSets()->getCount()); iRowIndex++)
    {
        float thisABMedianSignal = vABMedianSignal.get(iRowIndex);
        if ((isNaN(vAAMedianSignal.get(iRowIndex))) && (isNaN(vBBMedianSignal.get(iRowIndex))) && (!isNaN(thisABMedianSignal)))
        {
            vAAMedianSignal.set(iRowIndex, (float)(dAAMinusABMedian + thisABMedianSignal));
            vBBMedianSignal.set(iRowIndex, (float)(dBBMinusABMedian + thisABMedianSignal));
        }
    }
    // For each experiment calculate the allelic differences using the AA, AB, and BB median signal values.
    for (int iRowIndex = 0; (iRowIndex < CNAnalysisMethod::getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = CNAnalysisMethod::getProbeSets()->getAt(iRowIndex);
        pobjProbeSet->setAllelicDifference((float)NaN);
        if ((!isNaN(vAAMedianSignal.get(iRowIndex))) && (!isNaN(vABMedianSignal.get(iRowIndex))) && (!isNaN(vBBMedianSignal.get(iRowIndex))))
        {
            double dDifference = (pobjProbeSet->getAAlleleSignal() - pobjProbeSet->getBAlleleSignal());
            if (!Util::isFinite(dDifference)) {dDifference = 0.0;}
            float fGcAdjustment = pobjProbeSet->getGcAdjustment();
            if (fGcAdjustment != 0)
            {
                dDifference *= pow(2.0, -((double)fGcAdjustment));
                float fMedianAutosomeMedian = getExperiment()->getMedianAutosomeMedian();
                if (fMedianAutosomeMedian != 0)
                {
                    dDifference *= pow(2.0, -((double)fMedianAutosomeMedian));
                }
            }
            double dAMinusB = 0;
            if (dDifference < vABMedianSignal.get(iRowIndex))
            {
                dAMinusB = -(dDifference - vABMedianSignal.get(iRowIndex)) / (vBBMedianSignal.get(iRowIndex) - vABMedianSignal.get(iRowIndex));
            }
            else
            {
                dAMinusB = (dDifference - vABMedianSignal.get(iRowIndex)) / (vAAMedianSignal.get(iRowIndex) - vABMedianSignal.get(iRowIndex));
            }
            ADOutlierTrim(dAMinusB);
            pobjProbeSet->setAllelicDifference((float)dAMinusB);
        }
    }
}

void CNAnalysisMethodAllelicDifference::ADOutlierTrim(double& ad)
{
    if (fabs(ad) > m_dAllelicDifferenceOutlierTrim) {
        ad = numeric_limits<double>::quiet_NaN();
    }
}

void CNAnalysisMethodAllelicDifference::determineLocalProbeSets()
{
    if(m_bLocalProbeSetsDetermined)
    {
        return;
    }
    const bool isCyto2 = m_pEngine->getOptBool("cyto2");
    int iNumberOfProbeSets = CNAnalysisMethod::getProbeSets()->getCount();
    for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        CNProbeSet* p = CNAnalysisMethod::getProbeSets()->getAt(iIndex);
        if (isCyto2 && !p->getValidSCARExists())
        {
            continue;
        }
		if (isCyto2)
		{
			if(p->processAsSNP() && !isNaN(p->getAllelicDifference()))
			{
				getProbeSets()->add(p);
			}
		}
		else
		{
			if(p->processAsVisualization() && !isNaN(p->getAllelicDifference()))
			{
				getProbeSets()->add(p);
			}
		}
    }
    m_bLocalProbeSetsDetermined = true;
}

void CNAnalysisMethodAllelicDifference::filterCutOff(
                                                CNProbeSetArray* probeSets,
                                                int chrStart,
                                                int chrEnd,
                                                float cutoff3,
                                                float cutoff4,
                                                vector<bool>& filteredProbeSets
                                                )
{
    // No filtering
    for (int i = chrStart; i <= chrEnd; i++) {
        filteredProbeSets[i] = true;
    }
}


void CNAnalysisMethodAllelicDifference::filterNoMansLand(
                                            CNProbeSetArray* probeSets,
                                            vector<bool>& filteredProbeSets,
                                            const vector<pair<int, int> >& stepWindowBds,
                                            const vector<vector<double> >& peaks
                                            )
{
    filterNoMansLandImpl<CNAnalysisMethodAllelicDifference>(
                                                probeSets,
                                                filteredProbeSets,
                                                stepWindowBds,
                                                peaks
                                                );
}

void CNAnalysisMethodAllelicDifference::filterShrinkTowardPeaks(
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
    filterShrinkTowardPeaksImpl<CNAnalysisMethodAllelicDifference>(
                                                probeSets,
                                                filteredProbeSets,
                                                stepWindowBds,
                                                peaks,
                                                shrinkFactor3,
                                                shrinkFactor4,
                                                closestPeak,
                                                processedValues,
                                                saveAllelePeaksFlag
                                                );
}
