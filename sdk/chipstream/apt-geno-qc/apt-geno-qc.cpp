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


#include "chipstream/apt-geno-qc/GenoQC.h"
#include "chipstream/EngineUtil.h"
//

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
#ifdef _USE_VLD
#include "../external/vld/vld.h"
#endif
#endif
#endif


using namespace affxcdf;
using namespace affxcel;

using namespace std;

int main(int argc,const char* argv[]) 
{
  ofstream logOut;
  LogStream log;
  string logName;

  try {
    // make sure the file5 library has been started.
    affx::File5_open();
    
    GenoQC engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage(
    "apt-geno-qc - a single-chip-based genotyping chip quality control tool,\n"
                    "it reports one or more metrics for each chip analyzed.  The metrics\n"
                    "computed are determined by one of the supplied input files (the qca file).\n"
                    "\n"
                    "usage:\n"
                    "  apt-geno-qc --cdf-file my_chip.cdf --qcc-file my_chip.qcc \\\n"
                    "    --qca-file my_chip.qca --out-file my_output.txt \\\n"
                    "    --cel-files my_cel_files.txt\n"
                    "\n"
                    "  apt-geno-qc -c my_chip.cdf -q my_chip.qcc -a my_chip.qca \\\n"
                    "    --out-file my_output.txt --cel-files my_cel_files.txt\n");
    
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
    std::string outDir = Fs::dirname(engine.getOpt("out-file"));
    if(outDir == "") {
        outDir = ".";
    } else {
        if(!Fs::isWriteableDir(outDir)) {
            if(Fs::mkdirPath(outDir, false) != APT_OK) {
                Err::errAbort("Can't make or write to directory: " + outDir);
            }
        }
    }    
    // Set up the logging and message handlers.
    engine.openStandardLog("apt-geno-qc.log",logOut,log);
    
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
