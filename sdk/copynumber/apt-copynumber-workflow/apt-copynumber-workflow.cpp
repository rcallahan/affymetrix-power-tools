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
 * @file apt-copynumber-workflow.cpp
 *
 * @brief This file contains the housekeeping code for the apt-copynumber-workflow module.
 */

#ifdef _DEBUG
#ifndef _WIN64
#ifdef _USE_VLD
#include "../external/vld/vld.h"
#endif
#endif
#endif

//
#include "copynumber/CNAnalysisMethodChipstream.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNWorkflowEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"


using namespace std;

int main(int argc,const char* argv[]) 
{
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {
    // make sure the file5 library has been started.
    affx::File5_open();
    
    CNWorkflowEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage(
	"apt-copynumber-workflow - A program to compute copy number \n"
    "results from DNA analysis arrays.\n\n"
	"usage:\n"
	"	./apt-copynumber-workflow  \\\n"
    "       --adapter-type-normalization true \\\n"
    "       --text-output false \\\n"
    "       --reference-output CNReference.a5 \\\n"
    "       --set-analysis-name TestReference \\\n"
    "       --cdf-file GenomeWideSNP_6.cdf \\\n"
    "       --chrX-probes GenomeWideSNP_6.chrXprobes \\\n"
    "       --chrY-probes GenomeWideSNP_6.chrYprobes \\\n"
    "       --special-snps GenomeWideSNP_6.specialsnps \\\n"
    "       --netaffx-snp-annotation-file snp_annot_2.csv \\\n"
    "       --netaffx-cn-annotation-file cn_annot_2.csv \\\n"
    "       --o results --cel-files celfiles.txt");


    engine.parseArgv(argv);
    Verbose::setLevel(engine.getOptInt("verbose"));
    const string progName = Fs::basename(engine.getProgName());
    engine.setOpt("command-line",engine.commandLine());
    engine.setOpt("program-name",progName);
    engine.setOpt("program-company","Affymetrix");
    engine.setOpt("program-version",version);
    engine.setOpt("program-cvs-id",cvsId);
    engine.setOpt("version-to-report",versionToReport);
    engine.setOpt("exec-guid",execGuid);
    if(argc == 1) { engine.setOpt("help","true"); }
    
    // Check Options. Will print out version/help if requested then exit.
    engine.checkOptions();

    engine.openStandardLog("apt-copynumber-workflow.log",logOut,log);
    
    engine.run();
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      // Close log files
      logOut.close();
      return 1;
  }
  // Close log files
  logOut.close();

  return 0;
}
