////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#include "chipstream/apt-snp-model-converter/SnpModelConverterEngine.h"
//
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/Guid.h"
//
#include <fstream>
#include <string>
//


int main(int argc, const char* argv[])
{
  std::ofstream logOut;
  LogStream log(3, NULL, false);
  std::string logName;

  try {

    SnpModelConverterEngine engine;

    const std::string version = AptVersionInfo::version();
    const std::string cvsId = AptVersionInfo::cvsId();
    const std::string versionToReport = AptVersionInfo::versionToReport();
    const std::string execGuid = affxutil::Guid::GenerateNewGuid();


    /* Parse options. *///$$
    engine.setUsage(
      "apt-snp-model-converter - a tool for accessing chipstream/SnpModelConverter\n"
      "                  class directly. \n"
      "\n"
      "usage:\n"
      "  apt-snp-model-converter --dump-headers [raw|gtc4.1] my_models_file.[txt,a5,tsv,db]\n"
      "  NOTE: this utility currently only supports converting to sql database files."
      "\n"
    );

    engine.parseArgv(argv);
    Verbose::setLevel(engine.getOptInt("verbose"));
    const std::string progName = Fs::basename(engine.getProgName());
    engine.setOpt("command-line", engine.commandLine());
    engine.setOpt("program-name", progName);
    engine.setOpt("program-company", "Affymetrix");
    engine.setOpt("program-version", version);
    engine.setOpt("program-cvs-id", cvsId);
    engine.setOpt("version-to-report", versionToReport);
    engine.setOpt("exec-guid", execGuid);
    if(argc == 1) {
      engine.setOpt("help", "true");
    }

    // Check Options. Will print out version/help if requested then exit.
    engine.checkOptions();

    // Set up the logging and message handlers.
    engine.openStandardLog("apt-snp-model-converter.log", logOut, log);

    engine.run();
  } catch(...) {
    Verbose::out(1, "Unexpected Error: uncaught exception.");
    // Close log files
    logOut.close();
    return 1;
  }
  // Close log files
  logOut.close();

  return 0;
}
