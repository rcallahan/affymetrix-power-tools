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
   @file   MAS5Engine.cpp
   @author Alan Williams
   @date   Wed Apr 15 15:29:23 PDT 2009

   @brief Core routines for MAS5.
*/

//
//#include "ExpressionAlgorithmImplementation.h"
//#include "ExpressionControlsParameterExtraction.h"
#include "mas5-stat/apt-mas5/MAS5Engine.h"
//
#include "mas5-stat/workflow/MAS5ParameterExtraction.h"
//
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/EngineUtil.h"
#include "exp_report/src/ExpressionControlsParameterExtraction.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;

MAS5Engine::Reg MAS5Engine::reg;

MAS5Engine * MAS5Engine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == MAS5Engine::EngineName())
		return (MAS5Engine *)engine;
	return NULL;
}

MAS5Engine::MAS5Engine() {
    defineOptions();
}

MAS5Engine::~MAS5Engine() {
}

void MAS5Engine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
  defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets.",
                     "");
  defineOption("", "algo-params", PgOpt::STRING_OPT,
                     "File with MAS5 algorithm parameters.",
                     "");
  defineOption("", "control-params", PgOpt::STRING_OPT,
                     "File with MAS5 control parameters.",
                     "");
  defOptMult("", "chip-type", PgOpt::STRING_OPT,
                     "Chip types to check library and CEL files against. "
                     "Can be specified multiple times. "
                     "The first one is propigated as the chip type in the output files. "
                     "Warning, use of this option will override the usual check between chip types "
                     "found in the library files and cel files. You should use this option instead "
                     "of --force when possible. ",
                     "");

  defineOptionSection("Output Options");
  defineOption("", "cc-chp-output", PgOpt::BOOL_OPT,
                    "Output results in AGCC format.",
                    "true");
  defineOption("", "xda-chp-output", PgOpt::BOOL_OPT,
                    "Output results in XDA/GCOS format.",
                    "false");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "cels", PgOpt::STRING_OPT,
                     "Cel files to process.",
                     "");
  defOptMult("", "result-files", PgOpt::STRING_OPT,
                     "CHP file names to output. Must be paired with cels.",
                     "");
}

void MAS5Engine::defineStates() {
  defineOption("", "num-rows", PgOpt::INT_OPT,
                     "The number of rows on the chip.",
                     "-1");
  defineOption("", "num-cols", PgOpt::INT_OPT,
                     "The number of cols on the chip.",
                     "-1");
  defineOption("", "probe-count", PgOpt::INT_OPT,
                     "The number of probes on the chip.",
                     "-1");
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void MAS5Engine::checkOptionsImp() {

  defineStates();

  setLibFileOpt("cdf-file");
  setLibFileOpt("algo-params");
  setLibFileOpt("control-params");

  std::string cdfFile = getOpt("cdf-file");

  if(cdfFile == "")
      Err::errAbort("Must specify either a cdf file");

  bool ccOut = getOptBool("cc-chp-output");
  bool xdaOut = getOptBool("xda-chp-output");
  if(ccOut && xdaOut) 
      Err::errAbort("You can only set cc-chp-output or xda-chp-output. Not both.");

  if(!ccOut && !xdaOut)
      Err::errAbort("You must set either cc-chp-output or xda-chp-output.");

  string outDir = getOpt("out-dir");

  /* Read in cel file list from other file if specified. */
  vector<string> celFiles;
  EngineUtil::getCelFiles(celFiles, this);
  if(celFiles.size() == 0)
      Err::errAbort("No cel files specified.");
  setOpt("cels",celFiles);
    
  vector<string> resultFiles = getOptVector("result-files");
  if(resultFiles.size() > 0) {
      if(celFiles.size() != resultFiles.size())
        Err::errAbort("result-files option used but is not the same size as cel file listing");
  }

  // setup output folder
  if(!Fs::isWriteableDir(outDir)) {
    if(Fs::mkdirPath(outDir, false) != APT_OK) {
      Err::errAbort("Can't make or write to directory: " + ToStr(outDir));
    }
  }

  // Check chip types
  vector<string> chipTypesInLayout;

  /* Get the intial info about the chip and check cel files to make sure
     they match. */
  colrow_t numRows = 0, numCols = 0;
  int probeCount=0;
  int probeSetCount=0;
  if (cdfFile != "") {
    EngineUtil::getCdfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, cdfFile);
  }

  setOpt("num-rows", ToStr(numRows));
  setOpt("num-cols", ToStr(numCols));
  setOpt("probe-count", ToStr(probeCount));

  if (chipTypesInLayout.empty() || chipTypesInLayout[0] == "" || probeCount == 0) {
    Err::errAbort("Problem determining ChipType in file: " + cdfFile);
  }

  /* Did the user "force" a set of chip types via options? */
  vector<string> chipTypesSupplied = getOptVector("chip-type");

  /* Figure out what chip type to report */
  if(chipTypesSupplied.size() > 0) {
      setOpt("chip-type", chipTypesSupplied[0]);
  } else if(chipTypesInLayout.size() > 0) {
      setOpt("chip-type", chipTypesInLayout[0]);
  } else {
      Err::errAbort("Unable to figure out a chip type.");
  }

  /* Do Chip Type Check */
  if(!getOptBool("force")) {
    vector<string> chipTypeJustPrimary;
    vector<string> chipTypesToCheck;

    if(chipTypesSupplied.size() > 0) {
        chipTypesToCheck = chipTypesSupplied;
        EngineUtil::checkChipTypeVectors(chipTypesSupplied, chipTypesInLayout);
    } else {
        chipTypesToCheck = chipTypesInLayout;
    }

    chipTypeJustPrimary.push_back(chipTypesToCheck[0]);
    EngineUtil::checkCelChipTypes(chipTypesToCheck, probeCount, celFiles, numRows, numCols);

  }
}


void MAS5Engine::runImp() {
	m_Mas5.SaveToLegacyFile() = getOptBool("xda-chp-output");

    string algParams = getOpt("algo-params");
    if(algParams != "")
        MAS5ParameterExtraction::ExtractParameters(algParams.c_str(), m_Mas5.AlgParameters());

    string ctrlParams = getOpt("control-params");
    if(ctrlParams != "")
        ExpressionControlsParameterExtraction::ExtractParameters(
                ctrlParams.c_str(), m_Mas5.ProbePairThreshold(), m_Mas5.ReportControls());

    m_Mas5.ProgramName() = StringUtils::ConvertMBSToWCS(getOpt("program-name"));
    m_Mas5.ProgramCompany() = StringUtils::ConvertMBSToWCS(getOpt("program-version"));
    ///@todo what to report for ProgramId?
    //m_Mas5.ProgramId() = newVal;

    string baselineFile;
    string cdfFile = getOpt("cdf-file");
    vector<string> celFiles = getOptVector("cels");
    vector<string> resultFiles = getOptVector("result-files");
    Verbose::progressBegin(1, ToStr("Processing cel files"), celFiles.size(), 0, celFiles.size());
    for(size_t i=0; i<celFiles.size(); i++) {
        Verbose::progressStep(1);
        string celFile = celFiles[i];
        string resultFile;
        if(resultFiles.size() > 0) {
            resultFile = resultFiles[i];
        }
        else {
          resultFile = Fs::join(getOpt("out-dir"),Fs::basename(celFiles[i]));
          resultFile = Fs::noextname1(resultFile) + ".chp";
        }
	    bool success = m_Mas5.Execute(celFile, baselineFile, cdfFile, resultFile);
        if(!success)
            Err::errAbort("Failure in MAS5 execution");
    }
    Verbose::progressEnd(1, ToStr("Done."));
}


