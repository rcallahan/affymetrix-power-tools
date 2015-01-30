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
 * @file CNAnalysisMethodCNCyto2Gender.cpp
 *
 * @brief This file contains the CNAnalysisMethodCNCyto2Gender class members.
 */
#include "copynumber/CNAnalysisMethodCNCyto2Gender.h"
//
#include "file/TsvFile/TsvFile.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCNCyto2Gender::explainSelf()
{
    CNAnalysisMethodCNCyto2Gender obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCNCyto2Gender::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;
  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt Cutoff = {"cutoff", SelfDoc::Opt::Double, "0.5", "0.5", "0", "0.5", "Allele Peaks cutoff."};
  opts.push_back(Cutoff);

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
SelfCreate* CNAnalysisMethodCNCyto2Gender::
newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCNCyto2Gender * pMethod = new CNAnalysisMethodCNCyto2Gender();
    std::string strPrefix = getPrefix();

    pMethod->m_dCutoff = setupDoubleParameter("cutoff", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodCNCyto2Gender::CNAnalysisMethodCNCyto2Gender()
{
    m_dCutoff = 0;
        m_bLocalProbeSetsDetermined=false;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodCNCyto2Gender::run()
{
  Verbose::out(1, "MAJOR PROGRESS UPDATE: Calculating gender.");
  Verbose::out(1, "CNAnalysisMethodCNCyto2Gender::run(...) start");
  isSetup();

  determineLocalProbeSets();

  AffxArray<AffxString> vYProbeSetsToProcess;

    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    AffxString strProbeSetName;

    if ( tsv.open(m_pEngine->getOpt("chrY-probes")) == affx::TSV_OK ) {
      int colIx = tsv.getColumnCount(0) -1;
      while (tsv.nextLevel(0) == affx::TSV_OK){
        tsv.get(0,colIx, strProbeSetName);
        if (!strProbeSetName.empty()) {
          vYProbeSetsToProcess.add(new AffxString(strProbeSetName));
        }
      }
      tsv.clear();
    }

    vYProbeSetsToProcess.quickSort(0);

    CNProbeSetArray * allPsets = getProbeSets();

    // The gender call is all done from chromosome Y.  If it is there,
    // it is a male.  If it is not, it is a female.
    int not_zero = 0;
    int ycount = 0;
    pair<int,int>chrY = getChrBounds(m_iYChromosome, getProbeSets());
    for (int i=chrY.first; i<chrY.second; i++)
        {
        // Stay out of autosomal regions.
        if (allPsets->at(i)->isPseudoAutosomalRegion()) continue;
        AffxString strProbeSetName = allPsets->at(i)->getProbeSetName();
        if ((vYProbeSetsToProcess.getCount() > 0) && (vYProbeSetsToProcess.binarySearch(strProbeSetName, 0) == -1)) {continue;}

        // If the probe has copy number > 0 it hints toward a male.
        if (allPsets->at(i)->getCNState()) not_zero++;
        ycount++;
        }
    vYProbeSetsToProcess.deleteAll();

    // This may need to be fixed.  For now the algorithm is as simple
    // as possible.
    double confY = double(not_zero)/ycount;
	getExperiment()->setCNCallGenderNotZero(not_zero);
	getExperiment()->setCNCallGenderYCount(ycount);
    if (confY < m_dCutoff)
        {
        getExperiment()->setCNCallGender(affx::Female);
        getExperiment()->setCNCallGenderConfidence(1.0 - confY);
        getExperiment()->setCNCallGenderRatio(confY);
        getExperiment()->setCNCallGenderComputed(true);
        }

    else
        {
        getExperiment()->setCNCallGender(affx::Male);
        getExperiment()->setCNCallGenderConfidence(confY);
        getExperiment()->setCNCallGenderRatio(confY);
        getExperiment()->setCNCallGenderComputed(true);
        }

	Verbose::out(1, "CNAnalysisMethodCNCyto2Gender::run(...)\t" + getExperiment()->getExperimentName() + "\tCNGenderRatio\t" + ::getDouble(confY, 6) + "\t" + ((getExperiment()->getCNCallGender() == "male") ? "XY" : "XX") + "\tConfidence\t" + ::getDouble(getExperiment()->getCNCallGenderConfidence(), 6));
    Verbose::out(1, "CNAnalysisMethodCNCyto2Gender::run(...) end");
    }


void CNAnalysisMethodCNCyto2Gender::determineLocalProbeSets()
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

