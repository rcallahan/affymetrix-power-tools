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

#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/LogStream.h"

#include "AnnotationConverterEngine.h"

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
#define _CRT_SECURE_NO_WARNINGS // To disable localtime function warning
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

  try 
  {

    AnnotationConverterEngine engine;
    
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    
    /* Parse options. *///$$
    engine.setUsage(
    "apt-annotation-converter - A program for importing an annotation text files,\n"
    "                           into the annotation database. \n"
    "\n"
    "usage:\n"
    "Import a single file:\n"
    "apt-annotation-converter.exe -db-template Homo_sapiens.dbtemplate -db-file GenomeWideSNP_6.na30.annot.db  -tsv-file GenomeWideSNP_6.na30.annot.csv -array-set GenomeWideSNP_6 -array-config GenomeWideSNP_6.arrayconfig\n\n"
    "Import a list of files:\n"
    "apt-annotation-converter.exe -db-template Homo_sapiens.dbtemplate -db-file GenomeWideSNP_6.na30.annot.db  -tsv-files ListOfTsvToLoad.txt -array-set GenomeWideSNP_6 -array-config GenomeWideSNP_6.arrayconfig\n\n"
    "Import using the analysis_job job file:\n"
    "apt-annotation-converter.exe -xml-file  Axiom_GW_Hu.analysis_job\n");
    
    
    engine.programPath = argv[0]; 
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
    
    /*std::string outdir = Fs::dirname(engine.getOpt("summary-out-file"));
    Fs::ensureWriteableDirPath(outdir);*/

    // Set up the logging and message handlers.
    std::stringstream logFileName;
    
    time_t now = time(NULL);
    tm *timeinfo = localtime (&now);
    char buf[256];	
	strftime(buf, sizeof(buf), "%Y-%m-%dT%H_%M_%S", timeinfo);

    logFileName << "apt-annotation-converter." << buf << ".log";        
    engine.openStandardLog(logFileName.str(),logOut,log);
    
    engine.run();
  } 
  catch(...) {
      Verbose::out(0,"Unexpected Error: uncaught exception.");
      // Close log files
      logOut.close();
      return 1;
  }
  // Close log files
  logOut.close();
    
  return 0;
}
