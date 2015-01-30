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
 * @file apt-copynumber-cyto.cpp
 *
 * @brief This file contains the housekeeping code for the apt-copynumber-cyto module.
 */

#ifdef _DEBUG
#ifndef _WIN64
#ifdef _USE_VLD
#include "../../../external/vld/vld.h"
#endif
#endif
#endif

//
#include "copynumber/CNAnalysisMethodChipstream.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNCytoEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/TmpFileFactory.h"
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
    
    CNCytoEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage(
	    "apt-copynumber-cyto - A program to compute copy number \n"
        "results from DNA analysis arrays.\n\n"
	    "sample usage for CytoScanHD:\n"
	    "	./apt-copynumber-cyto "
            " --cyto2 false"
	    " --run-geno-qc true"
	    " <OPTIONS>"
        " \n\n"
        "sample usage for Cytogenetics_Array:\n"
	    "	./apt-copynumber-cyto"
            " --cyto2 true"
	    " --run-geno-qc false"
	    " <OPTIONS>"
        " \n\n"
        "where <OPTIONS> are as described bellow.\n"
        "(See apt-copynumber-cyto.html manual and vignettes for recommended parameters.)");
    
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

    engine.openStandardLog("apt-copynumber-cyto.log",logOut,log);
    
    // Check Options. Will print out version/help if requested then exit.
    try {engine.checkOptions();}
    catch(const exception& e) {
        Err::errAbort("Exception caught in main checkOptions(): " + ToStr(e.what()));
    }
    catch(...) {
        throw;
    }

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
