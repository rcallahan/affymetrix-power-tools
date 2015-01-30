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
 * @file   QuantMethodExprCHPReport.cpp
 * @author David Le
 * @date   Mon May 15 12:09:42 2006
 * 
 * @brief  Class for reporting results of quantification methods.
 */

//
#include "chipstream/QuantMethodExprCHPReport.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "file/CHPFileData.h"
#include "util/Fs.h"
//
#include <cstdlib>
#include <iostream>
//

using namespace std;

/** Constructor. */
QuantMethodExprCHPReport::QuantMethodExprCHPReport(AnalysisInfo& chpInfo, 
                                                   const std::string& prefix, 
                                                   const std::string& algName)
{
  m_Info = chpInfo;
  m_AlgName = algName;
  m_Prefix = prefix;
  m_CurrentProbeSetCount = 0;
  Err::check(chpInfo.m_ProbesetNames.size() == chpInfo.m_NumProbeSets, "Error: QuantMethodExprCHPReport::QuantMethodExprCHPReport() - m_NumProbeSets(" + ToStr(chpInfo.m_NumProbeSets) + ") != m_ProbesetNames.size(" + ToStr(chpInfo.m_ProbesetNames.size()) + ")");
}

/**
 * Set the CHP filenames to use for output. This overrides any prefix and allows
 * the caller to place the chp files anyplace they are desired.
 * @param fileNames - Vector of fileNames with full path to output
 * file, one for each cel file.
 */
void QuantMethodExprCHPReport::setChpFileNames(std::vector<std::string> &fileNames) {
  m_CHPFileNames = fileNames;
}

/** Destructor. */
QuantMethodExprCHPReport::~QuantMethodExprCHPReport() {
}

/** Swap out the .cel for .chp */
void QuantMethodExprCHPReport::setupFileNames(const IntensityMart &iMart) 
{
	m_CELFileNames = iMart.getCelFileNames();
  if (m_CHPFileNames.size() == 0) {
    m_CHPFileNames = Fs::changeDirAndExt(m_Prefix,m_CELFileNames,"."+m_Info.m_AlgName+".chp");
  }
  if (m_CHPFileNames.size() != iMart.getCelFileNames().size()) 
    Err::errAbort("Must be same number of output CHP files (" + 
                  ToStr(m_CHPFileNames.size()) + ") as input CEL files (" + ToStr(iMart.getCelFileNames().size()));
}

/** Fill in a chp header file info. */
void QuantMethodExprCHPReport::setupChpFile(affxchpwriter::CCHPFileWriter &chp, AnalysisInfo &info) 
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
bool QuantMethodExprCHPReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) 
{
	setupFileNames(iMart);

	// Make sure our output directory exists.
  Fs::ensureWriteableDirPath(m_Prefix, false);

	// Prepare headers for all CHP files.
    Verbose::progressBegin(1, "Initializing " + ToStr(m_CHPFileNames.size()) + " XDA CHP files", m_CHPFileNames.size(), 0, m_CHPFileNames.size());
	for(int target=0; target<m_CHPFileNames.size(); target++) 
	{
    Verbose::progressStep(1);

		// Remove any existing CHP file
		std::string chpName = m_CHPFileNames[target];
    Fs::rmIfExists(chpName, false);
    if (Fs::exists(chpName)) {
      Verbose::out(1,"Fs::exists returned true?");
      Err::errAbort("Unable to remove old CHP file, " + chpName);
    }

    // gen the name.
    std::string tmp_name=m_CHPFileNames[target] + ".tmp";

		// Create tmp chp file
		affxchpwriter::CCHPFileWriter chp;

    std::string tmp_unc_name=Fs::convertToUncPath(tmp_name);
		chp.SetFileName(tmp_unc_name.c_str());
		if(!chp.CreateNewFile()) {
      Err::errAbort("QuantMethodExprCHPReport::prepare() - Can't open CHP file '" +
                    FS_QUOTE_PATH(tmp_unc_name) + " to write."); 
    }
		chp.InitializeForWriting(m_Info.m_NumRows, 
                             m_Info.m_NumCols, 
                             m_Info.m_NumProbeSets, 
                             m_Info.m_ChipType.c_str(), 
                             m_Info.m_ProbeSetType, 
                             false);
		std::string fileRoot = Fs::basename(chpName);
		chp.SetParentCelFileName(fileRoot.c_str());
		setupChpFile(chp, m_Info);
		Err::check(chp.SaveHeader(), "QuantMethodExprCHPReport::prepare() - "
               "Couldn't save header for file: "+FS_QUOTE_PATH(tmp_unc_name));
		chp.Close();
    m_FilesForWriter.push_back(tmp_name);
	}
    Verbose::progressEnd(1, "Done.");

	// Initialize expression entry buffer writer.
	m_ExpressionEntryBufferWriter.Initialize(&m_FilesForWriter, false);	// false = Expression

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
bool QuantMethodExprCHPReport::report(ProbeSetGroup &psGroup, 
                                      QuantMethod &qMethod, 
                                      const IntensityMart &iMart, 
                                      std::vector<ChipStream *> &iTrans, 
                                      PmAdjuster &pmAdjust) 
{
	QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
	if(eMethod == NULL) { Err::errAbort("Can only use a QuantMethodExprReport with a QuantExprMethod."); }
	assert(!psGroup.probeSets.empty());

	// We only report expression probe sets.
	if(psGroup.probeSets[0]->psType != ProbeSet::Expression) { return false; }

	// Create an expression CHP result entry.
	affxchp::CExpressionProbeSetResults entry;
	for(int target=0; target<eMethod->getNumTargets(); target++) 
	{
		entry.Detection = ABS_NO_CALL;
		entry.DetectionPValue = 0.0f;
		entry.Signal = eMethod->getSignalEstimate(target);
		entry.NumPairs = m_Info.m_NumProbeSets;
		entry.NumUsedPairs = m_Info.m_NumProbeSets;
		m_ExpressionEntryBufferWriter.WriteExpressionEntry(target, entry);
	}

	m_CurrentProbeSetCount++;

	return true;
}

/** 
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodExprCHPReport::reportFailure(ProbeSetGroup &psGroup, 
                                             QuantMethod &qMethod, 
                                             const IntensityMart &iMart, 
                                             std::vector<ChipStream *> &iTrans, 
                                             PmAdjuster &pmAdjust) 
{
	QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
	if(eMethod == NULL) { Err::errAbort("Can only use a QuantMethodExprReport with QuantExprMethod."); }
	assert(!psGroup.probeSets.empty());

	// We only report expression probe sets.
	if(psGroup.probeSets[0]->psType != ProbeSet::Expression) { return false; }

	// Create a blank entry in order to be consistent with psGroups list.
	affxchp::CExpressionProbeSetResults entry;
	for(int target=0; target<m_CHPFileNames.size(); target++) 
	{
		entry.Detection = ABS_NO_CALL;
		entry.DetectionPValue = 0.0f;
		entry.Signal = 0.0f;
		entry.NumPairs = 0;
		entry.NumUsedPairs = 0;
		m_ExpressionEntryBufferWriter.WriteExpressionEntry(target, entry);
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
bool QuantMethodExprCHPReport::finish(QuantMethod &qMethod) 
{
	// Sanity to check we saw all the probe sets we were expecting.
	if(m_CurrentProbeSetCount != m_Info.m_NumProbeSets) 
	{
		Err::errAbort("QuantMethodExprCHPReport::finish() - Expecting: " + ToStr(m_Info.m_NumProbeSets) +
			" but got: " + ToStr(m_CurrentProbeSetCount) + ". GCOS XDA CHP file will be corrupt.");
	}

	// Flush remaining signal entries in the buffer.
	m_ExpressionEntryBufferWriter.FlushBuffer();

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
