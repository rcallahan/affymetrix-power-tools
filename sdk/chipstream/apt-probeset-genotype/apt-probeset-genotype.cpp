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
 * @file   apt-probeset-genotype.cpp
 * @author Chuck Sugnet
 * @date   Thu Jan 19 14:18:57 2006
 *
 * @brief  Initial version of program that uses chipstream architecture to
 *         analyze SNP chips.
 * @todo What happens to chp files when they are partially written? What should happen?
 * @todo More documentation for code, including high level summary. - Some done.
 * @todo Output "sampling" reports on summarization on qc reports.
 * @todo Read parameters from single file.
 * @todo Use robust estimators for cluster centers. - Not going to be implemented.
 * @todo Enable outputs easy to plot in R.
 * @todo Make gender estimates 'unknown' when no chrX SNPs are supplied.
 */

//
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-probeset-genotype/ProbesetGenotypeEngine.h"
#include "util/AptVersionInfo.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/MsgSocketHandler.h"
#include "util/Util.h"


/** Everybody's favorite function... */
int main(int argc,const char* argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;
//#ifndef WIN32

//#endif
  try {
    ProbesetGenotypeEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage(
                    "apt-probeset-genotype - program for determining genotype calls\n"
                    "from Affymetrix SNP microarrays. The model based algorithms for\n"
                    "making calls (brlmm/brlmm-p/birdseed) require multiple cel files\n"
                    "to be analyzed at once to learn the parameters for each SNP. \n"
                    "\n"
                    "usage:\n"
                    "  BRLMM (500K arrays):\n"
                    "    apt-probeset-genotype -c chip.cdf --chrX-snps chip.chrx \\ \n"
                    "         -o out-dir/ *.cel\n"
                    "\n"
                    "  BRLMM-P (GenomeWide SNP 5.0 arrays):\n"
                    "    apt-probeset-genotype -c chip.cdf --chrX-snps chip.chrx \\ \n"
                    "         -o out-dir/ -a brlmm-p --read-models-brlmmp chip.models \\\n"
                    "          *.cel\n"
                    "\n"
                    "  Birdseed (GenomeWide SNP 6.0 arrays):\n"
                    "    apt-probeset-genotype -c chip.cdf --special-snps chip.specialSNPs \\ \n"
                    "         -o out-dir/ -a birdseed --read-models-birdseed chip.birdseed.models \\\n"
                    "          *.cel\n"
                    "\n"
                    "  See the apt-probeset-genotype manual for more information about\n"
                    "  birdseed including the latest improvements from The Broad.\n"
                    );
    
    
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
    engine.checkDiskSpace();

    // Set up the logging and message handlers.
    engine.openStandardLog("apt-probeset-genotype.log",logOut,log);
    
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
