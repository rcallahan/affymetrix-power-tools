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
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-export-impute/ExportDataFilesEngine.h"
#include "util/AptVersionInfo.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"

#ifndef WIN32
#include <unistd.h>
#endif



/** Everybody's favorite function... */
int main(int argc,const char* argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {

    ExportDataFilesEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();

    
    /* Parse options. */
    engine.setUsage("apt-chp-to-table - A program for converting genotype CHP files to the\n"
                    "various table output files for calls and confidences.\n"
                    "usage:\n"
                    "   apt-chp-to-table ... \\\n"
                    "        ...");
    
    // Setup Options
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
    engine.openStandardLog("apt-chp-to-table.log",logOut,log);
    
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
