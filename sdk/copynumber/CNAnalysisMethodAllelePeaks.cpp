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
 * @file CNAnalysisMethodAllelePeaks.cpp
 *
 * @brief This file contains the CNAnalysisMethodAllelePeaks class members.
 */

//
#include "copynumber/CNAnalysisMethodAllelePeaks.h"
//
#include "chipstream/HomHiLoCelListener.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "util/Fs.h"
//
using namespace std;


/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodAllelePeaks::explainSelf()
{
    CNAnalysisMethodAllelePeaks obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodAllelePeaks::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt Step = {"step", SelfDoc::Opt::Integer, "50", "50", "1", "NA", "Allele Peaks step size."};
  opts.push_back(Step);

  SelfDoc::Opt Window = {"window", SelfDoc::Opt::Integer, "250", "250", "30", "NA", "Allele Peaks window."};
  opts.push_back(Window);

  SelfDoc::Opt PointCount = {"point-count", SelfDoc::Opt::Integer, "128", "128", "30", "1024", "Allele Peaks number of points."};
  opts.push_back(PointCount);

  SelfDoc::Opt Bandwidth = {"bandwidth", SelfDoc::Opt::Double, "0.45", "0.45", "0", "1", "Allele Peaks bandwidth."};
  opts.push_back(Bandwidth);

  SelfDoc::Opt Cutoff = {"cutoff", SelfDoc::Opt::Double, "0.05", "0.05", "0", "0.5", "Allele Peaks cutoff."};
  opts.push_back(Cutoff);

  SelfDoc::Opt CleanThreshold = {"clean-threshold", SelfDoc::Opt::Double, "0.25", "0.25", "0", "1.0", "Allele Peaks clean threshold."};
  opts.push_back(CleanThreshold);

  SelfDoc::Opt Symmetry = {"symmetry", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Allele Peaks SCAR mirror flag."};
  opts.push_back(Symmetry);

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
SelfCreate* CNAnalysisMethodAllelePeaks::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodAllelePeaks* pMethod = new CNAnalysisMethodAllelePeaks();
    std::string strPrefix = getPrefix();

    pMethod->m_iStep = setupIntParameter("step", strPrefix, params, doc, opts);
    pMethod->m_iWindow = setupIntParameter("window", strPrefix, params, doc, opts);
    pMethod->m_iPointCount = setupIntParameter("point-count", strPrefix, params, doc, opts);
    pMethod->m_dBandwidth = setupDoubleParameter("bandwidth", strPrefix, params, doc, opts);
    pMethod->m_dCutoff = setupDoubleParameter("cutoff", strPrefix, params, doc, opts);
    pMethod->m_dCleanThreshold = setupDoubleParameter("clean-threshold", strPrefix, params, doc, opts);
    pMethod->m_bSymmetry = setupBoolParameter("symmetry", strPrefix, params, doc, opts);

    return pMethod;
}


/**
 * @brief Constructor
 */
CNAnalysisMethodAllelePeaks::CNAnalysisMethodAllelePeaks()
{
    m_bSNPAlreadyLoaded = false;
    m_bLocalProbeSetsDetermined = false;
}


bool CNAnalysisMethodAllelePeaks::loadSnpReferenceFile(const AffxString& strFileName)
{
    Verbose::out(3, "CNAnalysisMethodAllelePeaks::loadSnpReferenceFile");
    bool bSuccessful = false;

    if (affx::File5_File::isHdf5file(strFileName)) {
        try {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);

            if(!group5->name_exists("FLDCutoff3peaks") || !group5->name_exists("FLDCutoff4peaks")) {
                return false;
            }

            affx::File5_Tsv* tsv5;

            tsv5 = group5->openTsv("FLDCutoff3peaks", affx::FILE5_OPEN);
            double dThreePeakFLD_X = 0.0;
            double dThreePeakFLD_Y = 0.0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &dThreePeakFLD_X);  m_vThreePeakFLD_X.push_back(dThreePeakFLD_X);
                tsv5->get(0, 1, &dThreePeakFLD_Y);  m_vThreePeakFLD_Y.push_back(dThreePeakFLD_Y);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("FLDCutoff4peaks", affx::FILE5_OPEN);
            double dFourPeakFLD_X = 0.0;
            double dFourPeakFLD_Y = 0.0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &dFourPeakFLD_X);  m_vFourPeakFLD_X.push_back(dFourPeakFLD_X);
                tsv5->get(0, 1, &dFourPeakFLD_Y);  m_vFourPeakFLD_Y.push_back(dFourPeakFLD_Y);
            }
            tsv5->close();
            delete tsv5;

            if(!group5->name_exists("ShrinkFactor3peaks") || !group5->name_exists("ShrinkFactor4peaks")) {
                return false;
            }

            tsv5 = group5->openTsv("ShrinkFactor3peaks", affx::FILE5_OPEN);
            double dThreePeakShrink_X = 0.0;
            double dThreePeakShrink_Y = 0.0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &dThreePeakShrink_X);  m_vThreePeakShrink_X.push_back(dThreePeakShrink_X);
                tsv5->get(0, 1, &dThreePeakShrink_Y);  m_vThreePeakShrink_Y.push_back(dThreePeakShrink_Y);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("ShrinkFactor4peaks", affx::FILE5_OPEN);
            double dFourPeakShrink_X = 0.0;
            double dFourPeakShrink_Y = 0.0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &dFourPeakShrink_X);  m_vFourPeakShrink_X.push_back(dFourPeakShrink_X);
                tsv5->get(0, 1, &dFourPeakShrink_Y);  m_vFourPeakShrink_Y.push_back(dFourPeakShrink_Y);
            }
            tsv5->close();
            delete tsv5;

            group5->close();
            delete group5;
            file5.close();
            bSuccessful = true;

        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }
    return bSuccessful;
}


void CNAnalysisMethodAllelePeaks::filterCutOff( CNProbeSetArray* probeSets,
                                                int chrStart,
                                                int chrEnd,
                                                float cutoff3,
                                                float cutoff4,
                                                vector<bool>& filteredProbeSets
                                                )
{
    // Filter based on the FLD cutoff functions
    for (int i = chrStart; i <= chrEnd; i++) {
      if (((probeSets->getAt(i)->getMaxPeaks() <= 3)  &&
           (probeSets->getAt(i)->getFLD() > cutoff3))  ||
          ((probeSets->getAt(i)->getMaxPeaks() >= 4)  &&
           (probeSets->getAt(i)->getFLD() > cutoff4)))
        {
            filteredProbeSets[i] = true;
        }
    }
}

void CNAnalysisMethodAllelePeaks::computeSNPQC(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets)
{
    // Calculate SNPQC
    if (!m_pEngine->getOptBool("cyto2")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-use-contrast")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-k")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-em-threshold")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-bin-size")) {return;}
    AffxArray<AffxString> arSNPIds;
    if ((m_pEngine->isOptDefined("snp-qc-snp-list")) && (m_pEngine->getOpt("snp-qc-snp-list") != ""))
    {
      affx::TsvFile tsv;
      std::string snp;
      tsv.m_optAutoTrim = true;
      if (tsv.openTable(m_pEngine->getOpt("snp-qc-snp-list")) == affx::TSV_OK) {
        while( tsv.nextLevel(0) == affx::TSV_OK ) {
          tsv.get(0,0, snp);
          arSNPIds.add(new AffxString(snp) );
        }
        tsv.clear();
      }
    }
    arSNPIds.quickSort(0);
    std::vector<double> vData;

    if (m_pEngine->getOptBool("snp-qc-use-contrast")) {
        if (arSNPIds.empty()) {
            for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                if (pobjProbeSet->processAsSNP()) {
                    vData.push_back(pobjProbeSet->getIntensityContrast(m_pEngine->getOptDouble("snp-qc-k")));
                }
            }
        } else {
            for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                if (pobjProbeSet->processAsSNP()) {
                    AffxString str = pobjProbeSet->getProbeSetName();
                    if (arSNPIds.binarySearch(str, 0) != -1) {
                        vData.push_back(pobjProbeSet->getIntensityContrast(m_pEngine->getOptDouble("snp-qc-k")));
                    }
                }
            }
        }
    } else {
        if (arSNPIds.empty()) {
            for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                if (pobjProbeSet->getInformation() != 0.0 && pobjProbeSet->processAsSNP()) {
                    vData.push_back(pobjProbeSet->getSCAR());
                }
            }
        } else {
            for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                if (pobjProbeSet->getInformation() != 0.0 && pobjProbeSet->processAsSNP()) {
                    AffxString str = pobjProbeSet->getProbeSetName();
                    if (arSNPIds.binarySearch(str, 0) != -1) {
                        vData.push_back(pobjProbeSet->getSCAR());
                    }
                }
            }
        }
    }
    arSNPIds.deleteAll();
    for (int i = 0; (i < vData.size()); i++)
    {
        if (vData[i] != vData[i])
        {
            vData[i] = 0;
        }
    }
    objExperiment.setSNPQC(HomHiLoCelListener::computeStatistic(vData, m_pEngine->getOptDouble("snp-qc-em-threshold"), m_pEngine->getOptDouble("snp-qc-bin-size")));
    objExperiment.setIsSNPQCset(true);

    Verbose::out(1, "*");
    Verbose::out(1, "CNQC  (MAPD) = " + ToStr(objExperiment.getMadDiffCN()));
    if (m_pEngine->getOptBool("snp-qc-use-contrast"))
    {
        Verbose::out(1, "SNPQC (Contrast) = " + ToStr(objExperiment.getSNPQC()));
    }
    else
    {
        Verbose::out(1, "SNPQC (PVQC) = " + ToStr(objExperiment.getSNPQC()));
    }
    Verbose::out(1, "*");
}


/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodAllelePeaks::run()
{
    Verbose::out(1, "CNAnalysisMethodAllelePeaks::run(...) start");
    isSetup();

    resetAllelePeakInitialValues();
    determineLocalProbeSets();

    if(!m_bSNPAlreadyLoaded) {
        ///@todo  Add error messaging if no snp-reference file is input. Also exit run gracefully.
        AffxString strFileName = m_pEngine->getOpt("snp-reference-input-file");
        loadSnpReferenceFile(strFileName);
        m_bSNPAlreadyLoaded = true;
    }
    // Compute SNPQC
    if (!getExperiment()->isSNPQCset()) {
        computeSNPQC(*getExperiment(), *getProbeSets());
    }
    int iStepSize = m_iStep;
    if (iStepSize > m_iWindow) {
        Verbose::warn(0, "Step size (" + ToStr(m_iStep) + ") exceeds window size (" + ToStr(m_iWindow) +
                      "). Using window size instead.");
        iStepSize = m_iWindow;
    }
    fillSnpChrBounds();

    shrinkToPeaks(getProbeSets());

    Verbose::out(1, "CNAnalysisMethodAllelePeaks::run(...) end");
}


/**
  * Find chromosome bounds for SNP probe sets
  */
void CNAnalysisMethodAllelePeaks::fillSnpChrBounds()
{
    fillChrBoundsImpl(getProbeSets(), m_chrBounds);
}


void CNAnalysisMethodAllelePeaks::determineLocalProbeSets()
{
    if(m_bLocalProbeSetsDetermined)
    {
            return;
    }
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    const bool isCyto2 = m_pEngine->getOptBool("cyto2");

    int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
    for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        CNProbeSet* p = CNAnalysisMethod::getProbeSets()->getAt(iIndex);
        if (isCyto2 && !p->getValidSCARExists())
        {
            continue;
        }
		if (isCytoScanHD)
		{
			if(p->processAsVisualization())
			{
				getProbeSets()->add(p);
			}
		}
		else
		{
			if(p->processAsSNP())
			{
				getProbeSets()->add(p);
			}
		}
    }
    m_bLocalProbeSetsDetermined=true;
}


CNProbeSetArray* CNAnalysisMethodAllelePeaks::getProbeSets()
{
        return &m_vLocalProbeSets;
}


void CNAnalysisMethodAllelePeaks::filterNoMansLand(
                                            CNProbeSetArray* probeSets,
                                            vector<bool>& filteredProbeSets,
                                            const vector<pair<int, int> >& stepWindowBds,
                                            const vector<vector<double> >& peaks
                                            )
{
    filterNoMansLandImpl<CNAnalysisMethodAllelePeaks>(
                                                probeSets,
                                                filteredProbeSets,
                                                stepWindowBds,
                                                peaks
                                                );
}


void CNAnalysisMethodAllelePeaks::filterShrinkTowardPeaks(
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
    filterShrinkTowardPeaksImpl<CNAnalysisMethodAllelePeaks>(
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
