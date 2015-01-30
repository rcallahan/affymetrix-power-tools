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
 * @file apt-copynumber-gender.cpp
 *
 * @brief This file contains the housekeeping code for the apt-copynumber-gender module.
 */

#ifdef _DEBUG
#ifndef _WIN64
#ifdef _USE_VLD
#include "../../../external/vld/vld.h"
#endif
#endif
#endif

//
#include "copynumber/CNGenderEngine.h"
//
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
    
    CNGenderEngine engine;
    
    /* Parse options. */
    engine.setUsage(
	    "apt-copynumber-gender - A program to compute gender. \n"
	    "sample usage for CytoScanHD:\n"
	    "	./apt-copynumber-gender \\\n"
	    "		 --verbose 3 \\\n"
	    "		 --chrX-probes ./lib/CytoScanHD_Array.chrXprobes \\\n"
	    "		 --chrY-probes ./lib/CytoScanHD_Array.chrYprobes \\\n"
	    "		 --male-gender-ratio 1.5 \\\n"
	    "		 --female-gender-ratio 0.9 \\\n"
	    "		 --cel-files ./celFile"
	    "		 --out-dir outputDir"
        " \n\n");
    engine.parseArgv(argv); 
    Verbose::setLevel(engine.getOptInt("verbose"));
    const string progName = Fs::basename(engine.getProgName());
/*
    engine.setOpt("command-line",engine.commandLine());
    engine.setOpt("program-name",progName);
    engine.setOpt("program-company","Affymetrix");
    engine.setOpt("program-version",version);
    engine.setOpt("program-cvs-id",cvsId);
    engine.setOpt("version-to-report",versionToReport);
    engine.setOpt("exec-guid",execGuid);
*/
    if(argc == 1) { engine.setOpt("help","true"); }

    engine.openStandardLog("apt-copynumber-gender.log",logOut,log);
    
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
