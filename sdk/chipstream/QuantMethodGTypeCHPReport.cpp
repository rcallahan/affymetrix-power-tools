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
 * @file   QuantMethodGTypeCHPReport.cpp
 * @author Chuck Sugnet
 * @date   Wed Mar 15 23:01:47 2006
 * 
 * @brief Reporter for genotyping probe sets that outputs chp files. There is
 * one CHP file produced for each cel file and the order of the CHP file is
 * guarenteed to be the same as the cdf file. It is very important that the
 * number and order of the probeset results is the exact same as the CDF
 * file. There have been a number of high impact bugs on this front and it is
 * important to be very careful.
 */

//
#include "chipstream/QuantMethodGTypeCHPReport.h"
//
#include "stats/stats.h"
#include "util/Fs.h"

using namespace affx;

/** Swap out the .cel for .chp */
void QuantMethodGTypeCHPReport::setupFileNames(const IntensityMart &iMart) 
{
  /* If the chp files names are empty chip off cel suffix and add "algname.chp" */
  m_CELFileNames = iMart.getCelFileNames();
  if (m_CHPFileNames.size() == 0) {
    m_CHPFileNames = Fs::changeDirAndExt(m_Prefix,m_CELFileNames,"."+m_Info.m_AlgName+".chp");
  }
  if(m_CHPFileNames.size() != iMart.getCelFileNames().size()) 
    Err::errAbort("Must be same number of output CHP files (" + 
      ToStr(m_CHPFileNames.size()) + ") as input CEL files (" + ToStr(iMart.getCelFileNames().size()));
}

/** Fill in a chp header file info. */
void QuantMethodGTypeCHPReport::setupChpFile(affxchpwriter::CCHPFileWriter &chp,
                                             AnalysisInfo &info) 
{
	assert(info.m_ParamNames.size() == info.m_ParamValues.size());
	chp.SetAlgorithmName(info.m_AlgName.c_str());
	chp.SetAlgorithmVersion(info.m_AlgVersion.c_str());
	chp.SetProgID(info.m_ProgID.c_str());
	chp.AddAlgorithmParameter("program-name", info.m_ProgramName.c_str());
	chp.AddAlgorithmParameter("program-version", info.m_ProgramVersion.c_str());
	chp.AddAlgorithmParameter("program-company", info.m_ProgramCompany.c_str());
	for(unsigned int i = 0; i < info.m_ParamNames.size(); i++) 
	{
		chp.AddAlgorithmParameter(info.m_ParamNames[i].c_str(), info.m_ParamValues[i].c_str());
	}
	// for(unsigned int i = 0; i < info.m_ClientInfoNames.size(); i++) 
	// {
	// 	chp.AddAlgorithmAppMetaInfo(info.m_ClientInfoNames[i].c_str(), info.m_ClientInfoValues[i].c_str());
	// }

}

/** 
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 * 
 * @param qMethod - Quantification method to be used.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeCHPReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {

	setupFileNames(iMart);
  
	// Make sure our output directory exists.
	if (!Fs::isWriteableDir(m_Prefix.c_str()) &&
            (Fs::mkdirPath(m_Prefix, false) != APT_OK)) {
          APT_ERR_ABORT("Can't make or write to directory: " + m_Prefix);
	}

	// Prepare headers for all CHP files.
    Verbose::progressBegin(1, "Initializing " + ToStr(m_CHPFileNames.size()) + " XDA CHP files", m_CHPFileNames.size(), 0, m_CHPFileNames.size());
	for(int i=0; i<m_CHPFileNames.size(); i++) 
	{
		// Remove any existing CHP file
		std::string chpName = m_CHPFileNames[i];
                if(Fs::isReadable(chpName.c_str())) {
                  if(Fs::rm(chpName, false)){ 
                    Err::errAbort("Unable to remove old CHP file, " + chpName);
                  }
                }

		// Create tmp chp file
        Verbose::progressStep(1);
		affxchpwriter::CCHPFileWriter chp;
    std::string tmp_name = m_CHPFileNames[i] + ".tmp";
    std::string tmp_unc_name=Fs::convertToUncPath(tmp_name);
		chp.SetFileName(tmp_unc_name.c_str());
		if(!chp.CreateNewFile()) 
		{
			Err::errAbort("QuantMethodGTypeCHPReport::prepare() - Can't open CHP file: "+
                    FS_QUOTE_PATH(tmp_unc_name)+" to write.");
		}
		chp.InitializeForWriting(
			m_Info.m_NumRows, 
			m_Info.m_NumCols, 
			m_Info.m_NumProbeSets,
			m_Info.m_ChipType.c_str(), 
			m_Info.m_ProbeSetType, 
			false);
		std::string fileRoot = Fs::basename(m_CELFileNames[i]);
		chp.SetParentCelFileName(fileRoot.c_str());
		setupChpFile(chp, m_Info);
		Err::check(chp.SaveHeader(), "QuantMethodGTypeCHPReport::prepare() - "
               "Couldn't save header for file '" + tmp_name);
		chp.Close();
    m_FilesForWriter.push_back(tmp_name);
	}
  Verbose::progressEnd(1, "Done.");

	// Initialize genotype entry buffer writer.
	m_GenotypeEntryBufferWriter.Initialize(&m_FilesForWriter, true);	// true = Genotype

	return true;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeCHPReport::report(ProbeSetGroup &psGroup, 
                                       QuantMethod &qMethod, 
                                       const IntensityMart &iMart, 
                                       std::vector<ChipStream *> &iTrans, 
                                       PmAdjuster &pmAdjust) {
	QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
	if(gMethod == NULL) { Err::errAbort("Can only use a QuantMethodGTypeReport with a QuantGTypeMethod."); }
	assert(psGroup.probeSets.size() > 0);

	// Check to make sure that we are getting probe sets in the right order.
	if(!Util::sameString(psGroup.probeSets[0]->name, m_Info.m_ProbesetNames[m_CurrentProbeSetCount])) 
	{
		Err::errAbort("QuantMethodGTypeCHPReport::report() - Expecting probeset: " + 
			ToStr(m_Info.m_ProbesetNames[m_CurrentProbeSetCount]) + ToStr(" got: ") +
			ToStr(psGroup.probeSets[0]->name) + " at index: " + ToStr(m_CurrentProbeSetCount));
	}

	// Create a genotype CHP entry from analysis results.
	affxchp::CGenotypeProbeSetResults entry;
	for(int target=0; target<gMethod->getNumCalls(); target++) 
	{
        ///@todo should check that the call encoding is the expected GType encoding
		GType call = GType_from_int(gMethod->getCall(target));
		double confidence = gMethod->getConfidence(target);
		// Translate our representation to the CHP call representation.
		switch (call) 
		{
			case NN: entry.AlleleCall = ALLELE_NO_CALL; break;
			case AA: entry.AlleleCall = ALLELE_A_CALL; break;
			case AB: entry.AlleleCall = ALLELE_AB_CALL; break;
			case BB: entry.AlleleCall = ALLELE_B_CALL; break;
			default: Err::errAbort("Don't recognize call of type: " + ToStr(call));
		}		
		entry.Confidence = confidence;

		/*  
		We do not modify the following to take integer values since the accessor function
		getConfGenotype accesses m_Distances which maintains data indexed by GType not int. 
		*/	
		entry.pvalue_AA = gMethod->getConfGenotype(AA, target);
		entry.pvalue_AB = gMethod->getConfGenotype(AB, target);
		entry.pvalue_BB = gMethod->getConfGenotype(BB, target);
		entry.pvalue_NoCall = gMethod->getConfGenotype(NN, target);
		m_GenotypeEntryBufferWriter.WriteGenotypeEntry(target, entry);
  }

  m_CurrentProbeSetCount++;

  return true;
}

/** 
 * For genotyping CHP files GTYPE likes to have the median of the raw intensity
 * probes stuffed into the confidence field of the CHP files. I'm not sure why
 * but Richard C asked for it so here it is.
 * 
 * @param psGroup - Groupt of probes.
 * @param iMart - Raw intensity data for chips.
 * @param chipIx - Which chip we are taking the median for.
 * 
 * @return median of all raw intensities for the PM probes in psGroup.
 */
float QuantMethodGTypeCHPReport::medianOfPmProbes(ProbeSetGroup &psGroup, const IntensityMart &iMart, int chipIx) {
  vector<float> pm;
  float med = 0;
  for(unsigned int psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    for(unsigned int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom *a = ps->atoms[atomIx];
      for(unsigned int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
        Probe *p = a->probes[probeIx];
        if(Probe::isPm(*p)) {
          pm.push_back(iMart.getProbeIntensity(p->id, chipIx));
        }
      }
    }
  }
  if(pm.empty()) {
    Err::errAbort("QuantMethodGTypeCHPReport::medianOfPmProbes() - No Pm Probes.");
  }
  med = median(pm.begin(), pm.end());
  return med;
}

/** 
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeCHPReport::reportFailure(ProbeSetGroup &psGroup, 
                                              QuantMethod &qMethod,
                                              const IntensityMart &iMart, 
                                              std::vector<ChipStream *> &iTrans, 
                                              PmAdjuster &pmAdjust) {
	QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
	if(gMethod == NULL) { Err::errAbort("Can only use a QuantMethodGTypeReport with a QuantGTypeMethod."); }
	assert(psGroup.probeSets.size() > 0);
	// Check to make sure that we are getting probe sets in the right order.
	if(!Util::sameString(psGroup.probeSets[0]->name,m_Info.m_ProbesetNames[m_CurrentProbeSetCount]))
	{
		Err::errAbort("QuantMethodGTypeCHPReport::reportFailure() - Expecting probeset: " + 
			ToStr(m_Info.m_ProbesetNames[m_CurrentProbeSetCount]) + ToStr(" got: ") +
			ToStr(psGroup.probeSets[0]->name) + " at index: " + ToStr(m_CurrentProbeSetCount));
	}

	// Create a blank entry in order to be consistent with psGroups list.
	affxchp::CGenotypeProbeSetResults entry;
	for (int target=0; target<m_CHPFileNames.size(); target++)
	{
		entry.AlleleCall = ALLELE_NO_CALL;
		entry.Confidence = medianOfPmProbes(psGroup, iMart, target);
		entry.RAS1 = 0.0f;
		entry.RAS2 = 0.0f;
		entry.pvalue_AA = 0.0f;
		entry.pvalue_AB = 0.0f;
		entry.pvalue_BB = 0.0f;
		entry.pvalue_NoCall = 0.0f;
		m_GenotypeEntryBufferWriter.WriteGenotypeEntry(target, entry);
	}

	m_CurrentProbeSetCount++;

	return true;
}

/** 
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * @param qMethod - Quantification method that was used.
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeCHPReport::finish(QuantMethod &qMethod) 
{
	// Sanity to check we saw all the probe sets we were expecting.
	if(m_CurrentProbeSetCount != m_Info.m_NumProbeSets) 
	{
		Err::errAbort("QuantMethodGTypeCHPReport::finish() - Expecting: " + ToStr(m_Info.m_NumProbeSets) +
			" but got: " + ToStr(m_CurrentProbeSetCount) + ". GCOS XDA CHP file wil be corrupt.");
	}

	// Flush remaining signal entries in the buffer.
	m_GenotypeEntryBufferWriter.FlushBuffer();

    // Remove .tmp extension
	for(unsigned int i = 0; i < m_CHPFileNames.size(); i++) 
	{
		std::string orig = m_CHPFileNames[i] + ".tmp";
        if(!Fs::fileRename(orig.c_str(),m_CHPFileNames[i].c_str())) {
            // Try and remove all the CHP files
            Fs::rm(orig, false);
	        for(unsigned int j = 0; j < m_CHPFileNames.size(); j++)
                    Fs::rm(m_CHPFileNames[j], false);

            Err::errAbort("Unable to rename " + orig + " to " + 
                        m_CHPFileNames[i]);
        }
	}
	return true;
}
