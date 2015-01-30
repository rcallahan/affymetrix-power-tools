////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   apt-summary-genotype.cpp
 * @author Alan Williams
 * @date   Fri Dec  7 09:10:30 PST 2007
 * 
 * @brief Program to generically transform cel file data using
 * chipstreams and write the data to a new cel file. Also serves as relatively simple
 * code sample to illustrate usage of chipstream code.
 */

//
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-summary-genotype/SummaryGenotypeEngine.h"
//
#include "util/AptVersionInfo.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"

//
#ifndef WIN32
#include <unistd.h>
#endif /* WIN32 */


/** Everybody's favorite function... */
int main(int argc,const char* argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {
    SummaryGenotypeEngine engine;
  
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage("apt-summary-genotype - Utility program to feed an allele summary file directly into\n"
		            "a genotype calling Quant method. Currenlty on brlmm-p is supported.\n"
		            "\n"
		            "usage:\n"
		            "   ./apt-summary-genotype -p CM=1.bins=100.K=2.SB=0.003.MS=0.05 \\\n"
		            "       -o results \\\n"
		            "       -s brlmm-p.plier.summary.txt \\\n"
		            "       --read-models-brlmmp GenomeWideSNP_5.models \\\n"
		            "       --chrX-snps GenomeWideSNP_5.chrx \\\n"
		            "       --read-genders genders.txt\n");
    
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
    engine.openStandardLog("apt-summary-genotype.log",logOut,log);
    
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
