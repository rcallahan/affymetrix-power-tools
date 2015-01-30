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
 * @file apt-copynumber-wave.cpp
 *
 * @brief This file contains the housekeeping code for the apt-copynumber-wave module.
 */

#ifdef _DEBUG
#ifndef _WIN64
#ifdef _USE_VLD
#include "../../../external/vld/vld.h"
#endif
#endif
#endif

//
#include "chipstream/EngineUtil.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNWaveEngine.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/TmpFileFactory.h"
#include "util/Util.h"


using namespace std;

/* main */
int main(int argc,const char* argv[]) 
{
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {
    // make sure the file5 library has been started.
    affx::File5_open();
    
    CNWaveEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */

    engine.setUsage(
      "apt-copynumber-wave - A program to compute additional copy number waves.\n\n"
      "usage:\n"
      " ./apt-copynumber-wave \\\n"
        "     --cn-reference-input /CytoFullV2.1-lib/na30.1Gold/Cytogenetics_Array.na30.1.v1.REF_MODEL\\\n"
        "     --cn-reference-output /CytoFullV2.1-lib/na30.1Gold/newref/Cytogenetics_Array.na30.1.v1.modified.REF_MODEL\\\n"
        "     --analysis additional-waves-reference-method.additional-wave-count=1.trim=2.0\\\n"
        "     .percentile=0.75.demean=false.cn-qc-cutoff=0.27.snp-qc-cutoff=1.1.force=false.keep-temp-data=false\\\n"
        "     .waviness-seg-count-cutoff=100.use-high-waviness-seg-count=true \\\n"
        "     --cychp-files /cychp_files/cychpList.txt\\\n"
        "     --temp-dir =/localData/WaveTemp");    
    
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
    try {engine.checkOptions();} catch(...) {}
    
    /* Set up the logging and message handlers. */
	engine.setOpt("out-dir", Fs::dirname(engine.getOpt("cn-reference-output")));
	AffxString str = engine.getOpt("out-dir");
	if ((str.endsWith("/")) || (str.endsWith("\\")))
	{
		str = str.substring(0, str.length() - 1);
		engine.setOpt("out-dir", str);
	}

  engine.openStandardLog("apt-copynumber-wave.log",logOut,log);
    
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
