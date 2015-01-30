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
#include "canary/CanaryEngine.h"
#include "canary/CanaryOptions.h"
#include "chipstream/EngineUtil.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/Util.h"
//
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>


///@todo need reference for Broad paper on canary

using namespace std;
using namespace affx;

int main(int argc, char *argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {
    CanaryEngine engine;

    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();

    /* Parse options. */
	engine.setUsage("apt-canary - Call copy number states for defined "
			"regions using the canary algorithm");

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

    engine.openStandardLog("apt-canary.log",logOut,log);
    
    // Run the engine. Will exit with non-zero exit statis if fails.
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
