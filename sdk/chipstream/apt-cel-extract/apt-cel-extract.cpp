////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
/// @file   apt-cel-extract.cpp
/// @brief  Main for extracting probe level intensities.
*/

//
#include "chipstream/apt-cel-extract/CelExtract.h"
//
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"

using namespace std;
using namespace affx;

CelExtractOptions ourOpts;


void define_aptcelextract_options(PgOptions* opts)
{
  opts->setUsage("apt-cel-extract - Extract probe level intensities from cel files.\n"
                 "Usage:\n"
                 "   apt-cel-extract -o intensities.txt \\\n"
                 "       [-c chip.clf -p chip.pgf] || [-d chip.cdf] \\\n"
                 "       [--probeset-ids norm-exon.txt] \\\n"
                 "       [--probe-ids probelist.txt] \\\n"
                 "       [--analysis analysis-string [--target-sketch-file sketchfile.txt]] \\\n"
                 "       [--pm-only] | [--pm-with-mm-only] \\\n"
                 "       [--report-background] \\\n"
                 "       [--cel-files celfiles.txt] | [ file1.CEL ... ] \n"
                 "\n"
                 "Synopsis:\n" 
                 "\n"
                 "   Simple extraction with no library files: \n"
                 "       apt-cel-extract -o intensities.txt *.CEL\n"
				 "  If the output file name includes the text \"SaturationReport\"\n"
				 "  then instead of dumping all the cel intensities, a cel saturation report is written.\n"
                 "\n"
                 "   Simple extraction of myprobes.txt from WT-based expression array: \n"
                 "       apt-cel-extract -o intensities.txt -c chip.clf \\\n"
                 "           -p chip.pgf --probeset-ids myprobes.txt *.CEL\n"
                 "\n"
                 "   Simple extraction from CDF-based expression array: \n"
                 "       apt-cel-extract -o intensities.txt -d chip.cdf *.CEL\\\n"
                 "\n"
                 "   Extract quantile(sketch) normalized PM with GCBG values: \n"
                 "       apt-cel-extract -o intensities.txt -c chip.clf \\\n"
                 "           -p chip.pgf -b chip.bgp -a quant-norm,pm-gcbg \\\n"
                 "           --report-background --cel-files celfiles.txt \n"
                 "\n");

  opts->defOptMult("", "probeset-ids", PgOpt::STRING_OPT,
                   "File containing probeset ids to extract probe level data for. "
                   "If no probeset-ids and no probe-ids file is provided, "
                   "information will be extracted for all probes. "
                   "May be specified multiple times.",
                   "");
  opts->defOptMult("", "probe-ids", PgOpt::STRING_OPT,
                   "File containing probe ids to extract probe level data for. "
                   "May be specified multiple times.",
                   "");
  opts->defineOption("c", "clf-file", PgOpt::STRING_OPT,
                     "The cel layout file, describing where a probe is within the cel file.",
                     "");
  opts->defineOption("p", "pgf-file", PgOpt::STRING_OPT,
                     "The probe group file, describing what probes are included in what probe sets.",
                     "");
  opts->defineOption("d", "cdf-file", PgOpt::STRING_OPT,
                     "Alternate method for describing probe sets. "
                     "Either -d or both -c and -p is required.",
                     "");
  opts->defineOption("s", "spf-file", PgOpt::STRING_OPT,
                     "Use simple probe file (SPF) for chip layout."
                     "[Experimental]",
                     "");
  opts->defineOption("b", "bgp-file", PgOpt::STRING_OPT,
                     "The background probes file, describing what probes are to be used for "
                     "gc background adjustment. Required if --pm-gcbg extraction is requested.",
                     "");
  opts->defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Optional file containing the names of cel files to extract data from.",
                     "");
  opts->defineOption("o", "out-file", PgOpt::STRING_OPT,
                     "Output file to contain the extraction output. If not provided, "
                     "the output will go to stdout (the console).",
                     "");
  opts->defineOption("","output-precision", PgOpt::INT_OPT,
                     "Number of decimal places in output.",
                     "2");
  opts->defineOption("", "target-sketch-file", PgOpt::STRING_OPT,
                     "A target sketch to normalize to when using quant-norm.",
                     "");
  opts->defineOption("a", "analysis", PgOpt::STRING_OPT,
                     "An analysis string (no quant method) to use to transform intensities. "
                     "[EXPERIMENTAL]",
                     "");
  opts->defineOption("", "log-file", PgOpt::STRING_OPT,
                     "The name of the log file to generate.",
                     "");
  opts->defineOption("f", "force", PgOpt::BOOL_OPT,
                     "Override sanity checks, for instance, requiring the same "
                     "lib file set/version in pgf, clf, and bgp files.",
                     "false");
  opts->defineOption("", "pm-only", PgOpt::BOOL_OPT,
                     "Only report PM probes. Requires chip layout information.",
                     "false");
  opts->defineOption("", "pm-with-mm-only", PgOpt::BOOL_OPT,
                     "Only report PM probes with have an MM. "
                     "Note only the PM probes are reported. "
                     "Requires chip layout information.",
                     "false");
  opts->defineOption("", "report-background", PgOpt::BOOL_OPT,
                     "Report the background value associated with each probe. "
                     "Requires an analysis string and chip layout information. "
                     "[EXPERIMENTAL]",
                     "false");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Print help message.",
                     "false");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("","temp-dir", PgOpt::STRING_OPT,
                     "Directory for temporary files when working off disk. Using network mounted drives is not advised. When not set, the output folder will be used.",
                     "");
  opts->defineOption("", "use-disk", PgOpt::BOOL_OPT,
                     "Store CEL intensities to be analyzed on disk.", "true");
  opts->defineOption("", "disk-cache", PgOpt::INT_OPT,
                     "Size of intensity memory cache in millions of intensities (when --use-disk=true).",
                     "50");

}

/** Fill in options given command line args **/
void fillInOptions(PgOptions *opts, CelExtractOptions &o, int argc) {
  //uint64_t free = 0, total = 0, swap = 0, avail = 0;
    o.version = AptVersionInfo::versionToReport();
    o.cvsId = AptVersionInfo::cvsId();
    o.progName = opts->getProgName();
    o.progName = Fs::basename(o.progName);
    o.execGuid = affxutil::Guid::GenerateNewGuid();
    o.timeStr = Util::getTimeStamp();
    if(opts->getBool("help") || argc == 1) {
        o.help = "true";
        return;
    }
    else
        o.help = "";

    if(opts->getBool("version")) {
        Verbose::out(0, "version: " + o.version);
        exit(0);
    }
    o.force = opts->getBool("force");
    o.verbosity = opts->getInt("verbose");
    o.diskDir = opts->get("temp-dir");
    o.useDisk = opts->getBool("use-disk");
    o.diskCache = opts->getInt("disk-cache");
    o.cdfFile = opts->get("cdf-file");
    o.pgfFile = opts->get("pgf-file");
    o.clfFile = opts->get("clf-file");
    o.spfFile = opts->get("spf-file");
    o.bgpFile = opts->get("bgp-file");
    //
    opts->mustFindOpt("probeset-ids")->push_user_values_into(o.probesetIdsFiles);
    opts->mustFindOpt("probe-ids")->push_user_values_into(o.probeIdsFiles);
    //
    o.outFile = opts->get("out-file");
    o.analysisString = opts->get("analysis");
    o.sketchInFile = opts->get("target-sketch-file");

    o.m_output_precision=opts->getInt("output-precision");

    o.pmOnly = opts->getBool("pm-only");
    o.ignoreProbesWithoutMM = opts->getBool("pm-with-mm-only");
    o.pairWithBackground = opts->getBool("report-background");

    /* Read in cel file list from other file if specified. */
    ///@todo this should go into Engine::Util as it is used in a number of apps
    if(opts->get("cel-files")!="") {
        affx::TsvFile tsv;
#ifdef WIN32
        tsv.m_optEscapeOk = false;
#endif
        std::string celFiles = opts->get("cel-files");
        string file;
        tsv.bind(0, "cel_files", &file, TSV_BIND_REQUIRED);
        if(tsv.open(celFiles) != TSV_OK) {
            Err::errAbort("Couldn't open cel-files file: " + celFiles);
        }
        tsv.rewind();
        while(tsv.nextLevel(0) == TSV_OK) {
            o.celFiles.push_back(Util::cloneString(file.c_str()));
        }
        tsv.close();
        Verbose::out(1, "Read " + ToStr(o.celFiles.size()) + " cel filenames from: " + Fs::basename(celFiles));
    }
    else {
        for(int i = 0; i < opts->getArgCount(); i++)
            o.celFiles.push_back(opts->getArg(i));
    }
}

/**
 * Report out some basics about what we think the run is:
 *
 * @param o - Program options.
 */
void reportBasics(CelExtractOptions &o) {
    Verbose::out(3, "version=" + o.version);
    Verbose::out(3, "command-line=" + o.commandLine);
    Verbose::out(3, "exec-guid=" + o.execGuid);
    Verbose::out(3, "time=" + o.timeStr);
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
}

/// @brief     The entry point for apt-cel-extract
/// @param     argc      arg count
/// @param     argv      arg strings
/// @return    
int main (int argc, char* argv[])
{
  try {
    PgOptions *opts = NULL;

    ofstream logOut;

    /* Parse options. */
    opts = new PgOptions();
    define_aptcelextract_options(opts);
    opts->parseArgv(argv);
    fillInOptions(opts, ourOpts, argc);
    Verbose::setLevel(ourOpts.verbosity);

    /* Print our help message if necessary. */
    if(ourOpts.help == "true" || argc == 1) {
        map<string,string>::iterator iter;
        opts->usage();
        cout << endl << "version: " << ourOpts.version << endl;
        exit(0);
    }
    /* Off we go to work. */
    else {
        time_t startTime = time(NULL);

        LogStream log;
        if(opts->get("log-file")!="") {
            string logName = opts->get("log-file");
            Fs::mustOpenToWrite(logOut, logName.c_str());
            log.setStream(&logOut);
            Verbose::pushMsgHandler(&log);
            Verbose::pushProgressHandler(&log);
            Verbose::pushWarnHandler(&log);
        }

  // Run the engine. Will exit with non-zero exit statis if fails.
        reportBasics(ourOpts);

        /* Do the work */
        CelExtract ce(ourOpts);
        ce.extract();

        /* Hmm, precision doesn't seem to work on cout, get 2 decimal places the hard way. */
        time_t endTime = time(NULL);
        int t = int(  (float)(endTime - startTime) / 60.0 * 100); // convert to minutes
        Verbose::out(1, ToStr("Run took approximately: ") + ToStr((float)t/100) + ToStr(" minutes."));

        if(opts->get("log-file")!="") {
            logOut.close();
        }
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
