////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
#include "chipstream/apt-impute/ImputeEngine.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/LogStream.h"
//
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
//

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#ifdef _DEBUG
#ifndef _WIN64
#include "../external/vld/vld.h"
#endif
#endif



using namespace std;

int main(int argc,const char* argv[]) 
{
  ofstream logOut;
  LogStream log;
  string logName;

  try {
    // make sure the file5 library has been started.
    
    ImputeEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    
    /* Parse options. *///$$
    engine.setUsage(
    "apt-snp-summary - a tool that generates a snp summary report in,\n"
    "                  text or binary format \n"
    "\n"
    "usage:\n"
    "  apt-snp-summary --summary-out-file my_output.txt --chp-files my_chp_files.txt\n"
    "\n"
    "  apt-snp-summary --summary-out-file my_output.bin --text-format false - \\\n"
    "    -chp-files my_chp_files.bin\n");
    
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
    
    // Set up the logging and message handlers.
    engine.openStandardLog("apt-impute.log",logOut,log);
    
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
