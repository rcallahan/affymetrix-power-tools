////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


//
#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/converters/cel/src/CELFileConverter.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AptVersionInfo.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifndef WIN32
#include <unistd.h>
#endif


using namespace std;             // Who isn't?
using namespace affymetrix_cel_converter;
using namespace affx;
using namespace affxutil;

class AOptions {
    public:
      bool help;
      int verbose;
      string outDir;
      bool inPlace;
      CELFileVersionType format;
      string progName;
      string version;
      string execGuid;
      string timeStr;
      string commandLine;
      string chipType;
      bool setDatName;
      int verbosity;
      vector<string> celFiles;
};

void define_aptcelconvert_options(PgOptions* opts)
{
  opts->setUsage("apt-cel-convert - program to convert cel files to different types.\n"
                 "usage:\n"
                 "   apt-cel-convert -f text -o text-cels *.CEL\n"
                 "   \n"
                 "   apt-cel-convert --format xda --out-dir text-cels --cel-files cel-files.txt\n");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Display program options and extra documentation about possible analyses.",
                     "false");
 opts->defineOption("i", "in-place", PgOpt::BOOL_OPT,
                    "Convert the file in place. Over-write existing file.",
                    "false");
 opts->defineOption("f", "format", PgOpt::STRING_OPT,
                    "Set the output cel file format type. Valid values: xda, text, agcc.",
                    "");
 opts->defineOption("", "log-file", PgOpt::STRING_OPT,
                    "The output log file. Defaults to location of output with name apt-cel-convert.log.",
                    "");
 opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                    "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                    "1");
 opts->defineOption("", "version", PgOpt::BOOL_OPT,
                    "Display version information.",
                    "false");
 opts->defineOption("", "set-dat-name", PgOpt::BOOL_OPT,
                    "Set the DAT file name to match that of the cel file name.",
                    "false");
 opts->defineOption("", "cel-files", PgOpt::STRING_OPT,
                    "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                    "");
 opts->defineOption("o", "out-dir", PgOpt::STRING_OPT,
                    "Directory to write result files into.",
                    "");
 opts->defineOption("", "chip-type", PgOpt::STRING_OPT,
                    "Force the new cel file to be this chip type.",
                    "");
}

/** Fill in options given command line arguments. */
void fillInOptions(PgOptions *opts, AOptions &o, int argc) {
  if(opts->getBool("help") || argc == 1) 
    o.help = true;
  else
    o.help = false;

  o.progName = opts->getProgName();
  o.progName = Fs::basename(o.progName);
  o.version = AptVersionInfo::versionToReport();
  o.execGuid = affxutil::Guid::GenerateNewGuid();

  //
  o.format = Unknown_Version;
  if (opts->get("format")=="text")
    o.format = GCOS_Version3;
  else if(opts->get("format")=="xda")
    o.format = GCOS_Version4;
  else if (opts->get("format")=="agcc")
    o.format = Calvin_Version1;
  else if (opts->get("format")!="")
    Err::errAbort("Invalid cel format: " + opts->get("format"));

  o.verbosity = opts->getInt("verbose");
  o.outDir = opts->get("out-dir");
  o.chipType = opts->get("chip-type");
  o.setDatName = opts->getBool("set-dat-name");
  o.inPlace = opts->getBool("in-place");
  
  /* Read in cel file list from other file if specified. */
  if (opts->get("cel-files")!="") {
    affx::TsvFile tsv;
#ifdef WIN32
    tsv.m_optEscapeOk = false;
#endif
    std::string celFiles = opts->get("cel-files");
    if(tsv.open(celFiles) != TSV_OK) {
      Err::errAbort("Couldn't open cell-files file: " + celFiles);
    }
    std::string file;
    tsv.bind(0, "cel_files", &file, TSV_BIND_REQUIRED);
    tsv.rewind();
    while(tsv.nextLevel(0) == TSV_OK) {
      o.celFiles.push_back(Util::cloneString(file.c_str()));
    }
    tsv.close();
    Verbose::out(1, "Read " + ToStr(o.celFiles.size()) + " cel files from: " + Fs::basename(celFiles));
  }
  else {
    for(int i = 0; i < opts->getArgCount(); i++) 
      o.celFiles.push_back(opts->getArg(i));
  }

  /* post processing. */
  o.timeStr = Util::getTimeStamp();

  if(!o.help) {
    /* Some sanity checks. */
    if(o.format == Unknown_Version) 
      Err::errAbort("Must provide a valid cel file format with --format option.");
    if(o.celFiles.empty()) 
      Err::errAbort("Must specify at least one cel file to convert.");
    if(Util::sameString(o.outDir.c_str(), "") && !o.inPlace)
      Err::errAbort("Must specify an output directory (--out-dir) or in place conversion (--in-place)");
    if(!Util::sameString(o.outDir.c_str(), "") && o.inPlace)
      Err::errAbort("Must specify an output directory (--out-dir) OR in place conversion (--in-place). Not both.");

    if(!o.inPlace) {
        if(!Fs::isWriteableDir(o.outDir.c_str())) {
            if(Fs::mkdirPath(o.outDir) != APT_OK) {
                Err::errAbort("Can't make or write to directory: " + ToStr(o.outDir));
            }
        }
    }
  }
}

/** 
 * Report out some basics about what we think the run is:
 * 
 * @param o - Program options.
 */
void reportBasics(AOptions &o, const string &version) {
  Verbose::out(3, "version=" + version);
#ifndef WIN32
  char name[8192];
  if(gethostname(name, ArraySize(name)) == 0) 
    Verbose::out(3, "host=" + ToStr(name));
  else 
    Verbose::out(3, "host=unknown");
  if(getcwd(name, ArraySize(name)) != NULL) 
    Verbose::out(3, "cwd=" + ToStr(name));
  else
    Verbose::out(3, "cwd=unknown");
#endif /* WIN32 */
  Verbose::out(3, "command-line=" + o.commandLine);
}

/** Everybody's favorite function... */
int main(int argc, char *argv[]) {
  try {
    const string version ("NON-OFFICIAL-RELEASE");
    ofstream logOut;
    string logName;
    PgOptions *opts = NULL;
    unsigned int i = 0;
    AOptions ourOpts;
    
    /* Parse options. */
    opts = new PgOptions();
    define_aptcelconvert_options(opts);
    opts->parseArgv(argv);
    /* Need to check for the version option before fillInOptions
        to avoid conflicts with those consistency checks */
    if(opts->getBool("version")) {
        cout << "version: " << version << endl;
        exit(0);
    }
    
    fillInOptions(opts, ourOpts, argc);
    ourOpts.commandLine = opts->commandLine();
    
    Verbose::setLevel(ourOpts.verbosity);
    // Do we need help? (I know I do...)
    if(ourOpts.help) {
        opts->usage();
        cout << "version: " << version << endl;
        exit(0);
    }
    else {
        time_t startTime = time(NULL);
    
        /* Set up the logging and message handlers. */
        if (opts->get("log-file") != "") {
            logName = opts->get("log-file");
        }
        else if (ourOpts.inPlace) {
          logName = Fs::join(".","apt-cel-convert.log");
        }
        else {
          logName = Fs::join(ourOpts.outDir,"apt-cel-convert.log");
        }
        Fs::mustOpenToWrite(logOut, logName.c_str());
        LogStream log(3, &logOut);
        Verbose::pushMsgHandler(&log);
        Verbose::pushProgressHandler(&log);
        Verbose::pushWarnHandler(&log);
    
        reportBasics(ourOpts, version);
    
        /* Do the heavy lifting */
	    CELFileConverter converter;
        for(i = 0; i < ourOpts.celFiles.size(); i++) {
    
            /* Figure out extra conversion options */
            CELFileConversionOptions convertOpts;
            // Copy (ie no format change) will fail if options are supplied -- so only provide
            // them if we need to
            CELFileConversionOptions *optsPtr = NULL;
            if(ourOpts.chipType != "") {
                convertOpts.m_ChipType = Util::cloneString(ourOpts.chipType.c_str());
                optsPtr = &convertOpts;
            }
            if(ourOpts.setDatName) {
                string newName = Fs::basename(ourOpts.celFiles[i]);
                newName = newName.substr(0,newName.rfind("."));
                convertOpts.m_DATFileName = Util::cloneString(newName.c_str());
                optsPtr = &convertOpts;
            }
    
            if(ourOpts.inPlace) {
                Verbose::out(1, "Converting " + ToStr(ourOpts.celFiles[i]));
	            if (converter.ConvertFile( ourOpts.celFiles[i].c_str(), ourOpts.format, optsPtr ) == false)
                    Err::errAbort("Could not convert cel file " + ToStr(ourOpts.celFiles[i]) + ": " + 
		                CELFileConverterErrorMessage(converter.ErrorCode()));
            } 
            else {
              string outfile = Fs::join(ourOpts.outDir,Fs::basename(ourOpts.celFiles[i]));
    
              Verbose::out(1, "Converting " + ToStr(ourOpts.celFiles[i]) + " to " + outfile);
    
	            if (converter.ConvertFile( ourOpts.celFiles[i].c_str(), outfile.c_str(), ourOpts.format, optsPtr ) == false)
                    Err::errAbort("Could not convert cel file " + ToStr(ourOpts.celFiles[i]) + ": " + 
		                CELFileConverterErrorMessage(converter.ErrorCode()));
            }
    
            if(convertOpts.m_DATFileName != NULL)
                delete[] convertOpts.m_DATFileName;
            if(convertOpts.m_ChipType != NULL)
                delete[] convertOpts.m_ChipType;
        }
    
        /* Close out log(s) */
        time_t endTime = time(NULL);
        int t = int(  (float)(endTime - startTime) / 60.0 * 100); // convert to minutes
        Verbose::out(1, ToStr("Run took approximately: ") + ToStr((float)t/100) + ToStr(((float)t/100) > 1 ? " minutes." : " minute."));
        if(!ourOpts.inPlace) 
            logOut.close();
    }
    
    delete opts;
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

