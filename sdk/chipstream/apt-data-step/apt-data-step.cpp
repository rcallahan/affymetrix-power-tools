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
/**
 * @file   apt-data-step
 * @author Chuck Sugnet
 */

//
#include "calvin_files/converters/cel/src/CELFileVersion.h"
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipStreamDataTransform.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/EngineUtil.h"
#include "file/CELFileWriter.h"
#include "file/TsvFile/ClfFile.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;

class AptCelTransformerOptions {
public:
  bool help;
  int verbose;
  string targetSketch;
  string m_a5_targetSketch;
  string chipstream;
  string probeNormFile;
  int cacheSize;
  string cdfFile;
  string pgfFile;
  string spfFile;
  string clfFile;
  bool force;
  string killListFile;
  bool writeSketch;
  bool m_a5_writeSketch;
  string outDir;
  string tempDir;
  string f5input;
  string f5output;
};

void define_aptceltransformer_options(PgOptions* opts)
{
  opts->setUsage("apt-data-step - Trying out chipstream stages\n");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("", "target-sketch", PgOpt::STRING_OPT,
                     "File specifying a target distribution to use for quantile normalization.",
                     "");
  opts->defineOption("", "a5-target-sketch", PgOpt::STRING_OPT,
                     "File specifying an target distribution to use for quantile normalization in A5 format.",
                     "");
  opts->defineOption("c", "chipstream", PgOpt::STRING_OPT,
                     "String representing tranformation desired. "
                     "For example: 'quant-norm.sketch=50000' does a quantile normalization using 50000 points.",
                     "");
  opts->defineOption("", "probe-norm-file", PgOpt::STRING_OPT,
                     "File specifying index of probes (1 based) to use for normalization.",
                     "");
  opts->defineOption("d", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets",
                     "");
  opts->defineOption("", "clf-file", PgOpt::STRING_OPT,
                     "File defining probes",
                     "");
  opts->defineOption("", "kill-list", PgOpt::STRING_OPT,
                     "File defining probes to remove",
                     "");
  opts->defineOption("", "pgf-file", PgOpt::STRING_OPT,
                     "File defining probe sets",
                     "");
  opts->defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "File defining probe sets",
                     "");
  opts->defineOption("", "write-sketch", PgOpt::BOOL_OPT,
                     "Write the quantile normalization distribution (or sketch) "
                     "to a file in output-dir for reuse with target-sketch option.",
                     "false");
  opts->defineOption("", "a5-write-sketch", PgOpt::BOOL_OPT,
                     "Write the quantile normalization distribution (or sketch) "
                     "to a file in output-dir for reuse with target-sketch option.",
                     "false");
  opts->defineOption("", "force", PgOpt::BOOL_OPT,
                     "Disable chip type checks.",
                     "false");
  opts->defineOption("o", "out-dir", PgOpt::STRING_OPT,
                     "Directory to write result files into.",
                     ".");
  opts->defineOption("","temp-dir", PgOpt::STRING_OPT,
                     "Directory for temporary files when working off disk. Using network mounted drives is not advised. When not set, the output folder will be used.",
                     "");
  opts->defineOption("", "input-file", PgOpt::STRING_OPT,
                     "Input data store file.",
                     "");
  opts->defineOption("", "output-file", PgOpt::STRING_OPT,
                     "Output data store file.",
                     "");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
}

/** 
 * Figure out chip dimensions from the cel file specified.
 * 
 * @param fileName - Name of cel file to read header from.
 * @param numCols - Number of x columns.
 * @param numRows - Number of y rows.
 */
void getDimensions(string fileName, colrow_t& numCols, colrow_t& numRows) {
  affxcel::CCELFileData cel;
  cel.SetFileName(fileName.c_str());
  if(!cel.ReadHeader()) 
    Err::errAbort("Can't read cel file: " + cel.GetFileName());
  numCols = cel.GetCols();
  numRows = cel.GetRows();
  cel.Close();
}



/** 
 * Transform all of the cel files specified using the chipstream modification stream
 * specified by the chipStreamStr and save the new modified cel files in the output directory.
 */
void transformCels(AptCelTransformerOptions &o) {
  /* get basic layout info */
  colrow_t numCols, numRows;
  int probeCount, psCount;
  vector<string> chipTypes;

  if(o.cdfFile != "") {
      EngineUtil::getCdfChipType(chipTypes, numRows, numCols, probeCount, psCount, o.cdfFile);
  }else if(o.pgfFile != "") {
      EngineUtil::getPgfChipType(chipTypes, numRows, numCols, probeCount, o.pgfFile, o.clfFile);
  }else if(o.spfFile != "") {
      EngineUtil::getSpfChipType(chipTypes, numRows, numCols, probeCount, psCount, o.spfFile);
  }

  /* Get the layout for the chip. */
  // Specifies probesets, locations of features on chip.
  ChipLayout layout; 
  // empty list of probesets, ie we want to load all
  std::set<const char *, Util::ltstr> probeSetsToLoad; 
  // empty list, ie load everything
  std::vector<bool> probeSubset; 
  // we have already check chip type, so use empty string;
  string chipType; 
  //
  probeidmap_t killList;

  if(o.cdfFile != "") {
    Verbose::out(1, "Opening cdf file: " + Fs::basename(o.cdfFile));
    if(!layout.openCdf(o.cdfFile, probeSetsToLoad, NULL, probeSubset, chipType, killList, false)) 
        Err::errAbort("Couldn't open cdf file: " + o.cdfFile);
  }else if(o.pgfFile != "") {
    affx::ClfFile clf;

    /* open clf file. */
    if(o.clfFile=="") 
        Err::errAbort("Must specify a clf file for chip.");
    Verbose::out(1, "Opening clf file: " +  Fs::basename(o.clfFile));
    if(!clf.open(o.clfFile))
        Err::errAbort("Couldn't open clf file: " + o.clfFile);
    layout.setDimensions(clf.getXMax() + 1, clf.getYMax() + 1);
    
    /* Make sure our assumptions about CLF file are true */
    Err::check(clf.getSequential() == 1,
                "ProbesetSummarizeEngine::loadPgfLayout() - unable to handle clf file without sequential set to 1.");
    Err::check(clf.getOrder().compare("col_major") == 0 || clf.getOrder().compare("row_major") == 0,
               "Unable to handle clf file without order set to row_major (old mislabeled 'col_major' accepted due to earlier bug.)");
   
    /* open pgf */
    Verbose::out(1, "Opening pgf file: " + Fs::basename(o.pgfFile));
    if(!layout.openPgf(o.pgfFile, clf.getXMax() + 1, clf.getYMax() + 1, probeSetsToLoad, NULL, NULL, probeSubset, chipType, killList, false, false))
        Err::errAbort("Couldn't open PGF file: " + o.pgfFile);
  }else if(o.spfFile != "") {
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    Verbose::out(1, "Opening spf file: " + Fs::basename(o.cdfFile));
    layout.openSpf(o.spfFile, probeSetsToLoad, NULL, probeSubset, chipType, false, psTypesToLoad);
  } else {
    Verbose::out(1,"Setting layout dimensions from cel file. Assuming all probes are PM.");
    layout.setDimensions(numCols, numRows);
    vector<bool> mask(layout.getProbeCount(), true);
    layout.setPmProbeMask(mask);
  }

  /* get list of probes to kill */
  /* we do this after the chip layout load as we don't want these to drop on the floor */
  if(o.killListFile != "")
    ChipLayout::fillInKillList(o.killListFile, killList, numRows, numCols);

  /* Setup chipstream factory */
  ChipStreamFactory cFactory;

  if(o.targetSketch != "") 
    cFactory.readTargetSketchFromFile(o.targetSketch);
  if(o.m_a5_targetSketch != "") 
    cFactory.readTargetSketchFromFile_a5(o.m_a5_targetSketch);
  if(o.writeSketch) 
    cFactory.setWriteSketchDir(o.outDir);
  if(o.m_a5_writeSketch) 
    cFactory.setWriteSketchDir_a5(o.outDir);
  if(o.probeNormFile != "")
    cFactory.readNormProbesFromFile(o.probeNormFile);
  if(o.killListFile != "") {
    /// @todo allow provision of kill list to Chipstream Factor
    ///cFactory.setKillList(killList);
    Err::errAbort("Kill list not yet implemented. Try probe-norm-file instead.");
  }

  /* Create our chipstream objects. */
  std::string chipStreamStr = o.chipstream;
  ChipStream *cStream = NULL;      // Our chipstream for transforming data.
  string dummy;
  cStream = cFactory.chipStreamForString(chipStreamStr, layout, dummy );
  ChipStreamDataTransform transform(cStream);
  DataStore in;
  in.initFromFile(o.f5input);
  DataStore out(o.f5output);
  out.initFromDataStore(in);
  PsBoard board;
  transform.transformData(board, in, out);
  Verbose::out(1, "Done.");
}


void fillInOptions(PgOptions *opts, AptCelTransformerOptions &o) {
  o.help = opts->getBool("help");
  o.verbose = opts->getInt("verbose");
  o.targetSketch = opts->get("target-sketch");
  o.chipstream = opts->get("chipstream");
  o.probeNormFile = opts->get("probe-norm-file");
  o.writeSketch = opts->getBool("write-sketch");
  o.m_a5_writeSketch = opts->getBool("a5-write-sketch");
  o.force = opts->getBool("force");
  o.outDir = opts->get("out-dir");
  o.tempDir = opts->get("temp-dir");
  if (o.tempDir == "") {
    o.tempDir = Fs::join(o.outDir,"temp"); 
  }
  Util::chompLastIfSep(o.outDir);
  o.cdfFile = opts->get("cdf-file");
  o.clfFile = opts->get("clf-file");
  o.pgfFile = opts->get("pgf-file");
  o.spfFile = opts->get("spf-file");
  o.killListFile = opts->get("kill-list");
  o.f5input = opts->get("input-file");
  o.f5output = opts->get("output-file");

  /* Sanity checks... */
  if(o.chipstream == "") 
    Err::errAbort("Must specify a chipstream transformer.");
  if(o.probeNormFile != "" && o.killListFile != "") 
    Err::errAbort("Cannot provide a kill list and a probe norm list.");
}

/** Everbody's favorite function... */
int main(int argc, const char *argv[]) {
  try {
    const string version = AptVersionInfo::versionToReport();

    AptCelTransformerOptions ourOpts;
    PgOptions *opts = new PgOptions();
    define_aptceltransformer_options(opts);
    opts->parseArgv(argv);
    
    if(opts->getBool("help") == true || argc == 1)
    {
        opts->usage();
        cout << "version: " << version << endl;
    }
    else if(opts->getBool("version"))
        cout << "version: " << version << endl;
    else {
        fillInOptions(opts, ourOpts);
        transformCels(ourOpts);
    }
    
    //
    delete opts;
    
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
