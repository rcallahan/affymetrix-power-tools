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

//
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-snp-summary/SnpSummaryEngine.h"
#include "chipstream/apt-snp-summary/SnpSummaryReporter.h"
#include "util/Fs.h"
#include "util/FsPath.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
//
SnpSummaryEngine::Reg SnpSummaryEngine::reg;

SnpSummaryEngine * SnpSummaryEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == SnpSummaryEngine::EngineName())
		return (SnpSummaryEngine *)engine;
	return NULL;
}

// constructor
SnpSummaryEngine::SnpSummaryEngine() {
  defineOptions();
}

SnpSummaryEngine::~SnpSummaryEngine() {
  clear();
}

void SnpSummaryEngine::clear()
{
	chpFiles.clear();  
}

void SnpSummaryEngine::defineOptions() 
{

  defineOptionSection("Input Options");
  defineOption("", "chp-files", PgOpt::STRING_OPT,
                     "Text file specifying chp files to process, one per line with the first line being 'chp_files'.",
                     "");

  defineOption("", "call-file", PgOpt::STRING_OPT,
                     "TSV calls file (table input) from apt-probeset-genotype.",
                     "");

   defineOption("b", "buffer-size", PgOpt::INT_OPT,
                    "Size of buffer used to process files",
                    "10000");
 
  defineOptionSection("Output Options");

  defineOption("", "output-format", PgOpt::STRING_OPT,
	  "Set the format of output file:  Valid vales: txt for text output, db for SQLite database format, bin for binary.",
                     "txt");

  defineOption("s", "summary-out-file", PgOpt::STRING_OPT,
                     "Name to use for the output file.",
                     "apt-snp-summary.report.txt");

  defineOption("", "files-col", PgOpt::STRING_OPT,
		"The column name for the input files.",
        "chp_files");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "chps", PgOpt::STRING_OPT,
                     "CHP files to process.",
                     "text output");
}

void SnpSummaryEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void SnpSummaryEngine::checkOptionsImp()
{
    defineStates();
	string file;
	file = getOpt("summary-out-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output summary file.");
	if (getOpt("output-format") == "")
		setOpt("output-format", "txt");
	std::string outputFormat = getOpt("output-format");
	if ((outputFormat != "txt") && (outputFormat != "bin") && (outputFormat != "db"))
		Err::errAbort("Must specify the the format of the output file.  Set the --output option to txt for text output, db for SQLite database format, bin for binary.");
    if((outputFormat == "bin") && (getOpt("call-file") != "")) {
        Err::errAbort("Cannot use table input with bin output.");
    }

	// Read in list files list from other file if specified.
    vector<string> chpFiles;
    EngineUtil::getChpFiles(chpFiles, this);
    for(int i=0; i < chpFiles.size(); i++) {
      FsPath chpPath(chpFiles[i]);
      if(!chpPath.isReadable()) {
        Err::errAbort("Unable to read " + chpFiles[i] + " " + Fs::getErrMsg());
      }
    }

    if((chpFiles.size() == 0) && getOpt("call-file")=="")
        Err::errAbort("No chp files specified.");

	setOpt("chps",chpFiles);
}

// Compute the various metrics over all chips.
void SnpSummaryEngine::runImp() {

	SnpSummaryReporter reporter;
	int size = getOptInt("buffer-size");

	bool rc = reporter.CreateSnpSummaryReport(getOpt("summary-out-file"), getOptVector("chps"), getOpt("output-format"), size, getOpt("call-file"));	
    if(!rc)
        Err::errAbort("Unable to extract SNP summaries.");
}




