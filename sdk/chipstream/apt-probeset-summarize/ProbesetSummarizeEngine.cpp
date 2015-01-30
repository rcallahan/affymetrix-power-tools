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
 * @file   ProbesetSummarizeEngine.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 12:21:35 2005
 * 
 * @brief  Program for summarizing probeset results from the command line. Requires
 * three main inputs: 1) layout file decribing individual probes and their grouping into 
 * probe sets. 2) clf file telling where probes are located on chip. 3) Data input files,
 * usually cel files.
 *
 * @todo runImp is getting long, time to refactor.
 * @todo Windows documentation and example quantification method.
 * @todo Figure out what debugging information necessary for remote troubleshooting and include it somewhere.
 * @todo Enable quantile normalization withing groups and median normalization across groups.
 * @todo DABG flag for -logP output
 * @todo Only read in required columns from map files.
 * @todo For QC Reports: consider making the P-value threshold for detection configurable
 *       from the command line, with the default set from the pgf file.
 */

//
#include "chipstream/apt-probeset-summarize/ProbesetSummarizeEngine.h"
//
#include "chipstream/CelReader.h"
#include "chipstream/CelStatListener.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/SparseMart.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/MetaProbeset.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/QuantMethodExprCCCHPReport.h"
#include "chipstream/QuantMethodExprCHPReport.h"
#include "chipstream/QuantMethodExprChipSummary.h"
#include "chipstream/QuantMethodMultiDataCCCHPReport.h"
#include "chipstream/QuantMethodRunReport.h"
#include "chipstream/QuantMethodSamplerReport.h"
//
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/TsvFile/ClfFile.h"
#include "file/TsvFile/TsvFile.h"

#include "util/Fs.h"
//
#include "newmat.h"


using namespace std; // Who isn't?
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affx;

ProbesetSummarizeEngine::Reg ProbesetSummarizeEngine::reg;

ProbesetSummarizeEngine * ProbesetSummarizeEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == ProbesetSummarizeEngine::EngineName())
        return (ProbesetSummarizeEngine *)engine;
    return NULL;
}

ProbesetSummarizeEngine::ProbesetSummarizeEngine() {
    m_a5_global_input_file = NULL;
    m_a5_global_input_group = NULL;
    m_a5_global_output_file = NULL;
    m_a5_global_output_group = NULL;
    defineStdMethods();
    defineOptions();
}

ProbesetSummarizeEngine::~ProbesetSummarizeEngine() {
}

void ProbesetSummarizeEngine::defineStdMethods() {
    m_stdMethods["plier-gcbg"] = "quant-norm.sketch=0.bioc=false,pm-gcbg,plier";
    m_stdMethods["plier-mm"] = "quant-norm.sketch=0.bioc=false,pm-mm,plier";
    m_stdMethods["rma"] = "rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish";
    m_stdMethods["plier-gcbg-sketch"] = "quant-norm.sketch=-1.bioc=false,pm-gcbg,plier";
    m_stdMethods["plier-mm-sketch"] = "quant-norm.sketch=-1.bioc=false,pm-mm,plier";
    m_stdMethods["rma-sketch"] = "rma-bg,quant-norm.sketch=-1.usepm=true.bioc=true,pm-only,med-polish";
    m_stdMethods["dabg"] = "pm-only,dabg";
}

void ProbesetSummarizeEngine::defineOptions() {
    defineOptionSection("Input Options");
    defineOption("", "cel-files", PgOpt::STRING_OPT,
                 "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                 "");
    defineOption("d", "cdf-file", PgOpt::STRING_OPT,
                 "File defining probe sets. Use either --cdf-file, --spf-file, or --pgf-file and --clf-file. "
                 "Automatically sets --names.",
                 "");
    defineOption("", "spf-file", PgOpt::STRING_OPT,
                 "File defining probe sets in spf (simple probe format) which is like a text cdf file.",
                 "");
    defineOption("p", "pgf-file", PgOpt::STRING_OPT,
                 "File defining probe sets.",
                 "");
    defineOption("c", "clf-file", PgOpt::STRING_OPT,
                 "File defining x,y <-> probe id conversion. Required when using PGF file.",
                 "");
    defineOption("b", "bgp-file", PgOpt::STRING_OPT,
                 "File defining probes to be used for GC background.",
                 "");
    defineOption("s", "probeset-ids", PgOpt::STRING_OPT,
                 "File specifying probe sets to summarize.",
                 "");
    defineOption("m", "meta-probesets", PgOpt::STRING_OPT,
                 "File containing meta probeset definitions. "
                 "File must contain a probeset_id column and a probeset_list column.",
                 "");
    defineOption("", "probe-class-file", PgOpt::STRING_OPT,
                 "File containing probe_id (1-based) of probes and a 'class' designation. "
                 "Used to compute mean probe intensity by class for report file.",
                 "");
    defineOption("", "qc-probesets", PgOpt::STRING_OPT,
                 "File with probeset_id(name) and group_name columns specifying subsets of probesets to compute qc stats for.",
                 "");
    defOptMult("", "chip-type", PgOpt::STRING_OPT,
               "Chip types to check library and CEL files against. "
               "Can be specified multiple times. "
               "The first one is propigated as the chip type in the output files. "
               "Warning, use of this option will override the usual check between chip types "
               "found in the library files and cel files. You should use this option instead "
               "of --force when possible. ",
               "");

    defineOptionSection("Output Options");
    defineOption("", "use-pgf-names", PgOpt::BOOL_OPT,
                 "Use the probeset_names instead of probeset_id column in the PGF file for output.",
                 "false");
    defineOption("", "cc-chp-output", PgOpt::BOOL_OPT,
                 "Output results in directory called 'cc-chp' under out-dir. "
                 "This makes one AGCC Expression CHP file per cel file analyzed.",
                 "false");
    defineOption("", "xda-chp-output", PgOpt::BOOL_OPT,
                 "Output resulting calls in directory called 'chp' under out-dir. "
                 "This makes one GCOS XDA CHP file per cel file analyzed.",
                 "false");
    defineOption("", "cc-md-chp-output", PgOpt::BOOL_OPT,
                 "Output resulting calls in directory called 'cc-md-chp' under out-dir. "
                 "This makes one AGCC Multi Data CHP file per cel file analyzed.",
                 "false");
    defineOption("", "cc-chp-out-dir", PgOpt::STRING_OPT,
                 "Over-ride the default location for chp output.",
                 "");
    defineOption("", "xda-chp-out-dir", PgOpt::STRING_OPT,
                 "Over-ride the default location for chp output.",
                 "");
    defineOption("", "cc-md-chp-out-dir", PgOpt::STRING_OPT,
                 "Over-ride the default location for chp output.",
                 "");
    defineOption("", "subsample-report", PgOpt::BOOL_OPT,
                 "Output subsamples of the data intensities, summaries and residuals for error checking downstream.",
                 "false");
    defineOption("", "report-file", PgOpt::STRING_OPT,
                 "Over-ride the default report file name.",
                 "");

    defineOptionSection("Analysis Options");
    defOptMult("a", "analysis", PgOpt::STRING_OPT,
               "String representing analysis pathway desired. For example: 'quant-norm,pm-gcbg,plier'. "
               "Prepackaged analysis such as 'plier-gcbg-sketch', 'plier-gcbg', 'plier-mm-sketch', "
               "'plier-mm', 'rma-sketch', and 'rma' can be specified. "
               "Multiple analysis allowed at same time. "
               "When using quantile normalization, you may need to use the sketch option "
               "to avoid running out of memory.",
               "");
    defineOption("", "summaries", PgOpt::BOOL_OPT,
                 "Output expression summaries in text table format.",
                 "true");
    defineOption("", "feat-effects", PgOpt::BOOL_OPT,
                 "Output feature effects when available. By convention med-polish feature effects have total probeset median added to them, see RMA module for details",
                 "false");
    defineOption("", "use-feat-eff", PgOpt::STRING_OPT,
                 "File defining a plier feature effect for each probe. "
                 "Note that precomputed effects should only be used for an appropriately similar analysis "
                 "(i.e. feature effects for pm-only may be different than for pm-mm). "
                 "Currently a probe is expected to be in only a single probeset. "
                 "This option does not work for IterPlier or SEA.",
                 "");
    defineOption("", "feat-details", PgOpt::BOOL_OPT,
                 "Output probe by chip specific details (often residuals) when available.",
                 "false");
    defineOption("", "target-sketch", PgOpt::STRING_OPT,
                 "File specifying a target distribution to use for quantile normalization.",
                 "");
    defineOption("", "write-sketch", PgOpt::BOOL_OPT,
                 "Write the quantile normalization distribution (or sketch) to a file for reuse with target-sketch option. "
                 "WARNING: If more than one -a option generates a target sketch file, it is not deterministic which file "
                 "will be retained by the OS if the target sketch files have the same name.",
                 "false");
    defineOption("","reference-profile",PgOpt::STRING_OPT, "Reference profile",
                 "");
    defineOption("","write-profile", PgOpt::BOOL_OPT,
                 "write reference profile.","false");
    defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
                 "Explicitly set the analysis name. "
                 "This affects output file names (ie prefix) and various meta info.",
                 "");
    defineOption("x", "precision", PgOpt::INT_OPT,
                 "How many digits of precision to use after decimal.",
                 "5");

    defineOptionSection("Misc Options");
    defineOption("", "explain", PgOpt::STRING_OPT,
                 "Explain a particular operation (i.e. --explain rma-bg).",
                 "");

    defineOptionSection("Advanced Options");
    defineOption("", "kill-list", PgOpt::STRING_OPT,
                 "Do not use the PM probes specified in file for computing results. "
                 "[experimental]",
                 "");

    defineOptionSection("Execution Control Options");
    defineOption("", "use-disk", PgOpt::BOOL_OPT,
                 "Store CEL intensities to be analyzed on disk.", "true");
    defineOption("", "disk-cache", PgOpt::INT_OPT,
                 "Size of intensity memory cache in millions of intensities (when --use-disk=true).",
                 "50");
    defineOption("", "store-duplicate-probes", PgOpt::BOOL_OPT, "Store intensities for probes appearing in multiple probesets in memory (Prevents page thrashing.  Is a bad idea for Axiom.  Turned on automatically when using meta-probesets)","false");
    defineOptionSection("A5 output options");

    defineOption("","a5-global-file",PgOpt::STRING_OPT,
                 "Filename for the A5 global output file. "
                 "[Experimental]",
                 "");
    defineOption("","a5-global-file-no-replace",PgOpt::BOOL_OPT,
                 "Append or create rather than replace. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-group",PgOpt::STRING_OPT,
                 "Group name where to put results in the A5 output files. "
                 "Defaults to '/'. "
                 "[Experimental]",
                 "");
    defineOption("", "a5-summaries", PgOpt::BOOL_OPT,
                 "Output the summary values from the quantifcation method for each allele in A5 format. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-summaries-use-global",PgOpt::BOOL_OPT,
                 "Use the global A5 file for summaries. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-feature-effects",PgOpt::BOOL_OPT,
                 "Output feature effects in A5 format. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-feature-effects-use-global",PgOpt::BOOL_OPT,
                 "Use the global A5 file for feature effects."
                 "[Experimental]",
                 "false");
    defineOption("","a5-feature-details",PgOpt::BOOL_OPT,
                 "Output feature level residuals in A5 format. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-feature-details-use-global",PgOpt::BOOL_OPT,
                 "Use the global A5 file for residuals. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-sketch",PgOpt::BOOL_OPT,
                 "Output normalization sketch in A5 format. "
                 "--write-sketch option will override this option. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-sketch-use-global",PgOpt::BOOL_OPT,
                 "Put the sketch in the global A5 output file. "
                 "[Experimental]",
                 "false");
  
    defineOptionSection("A5 input options");

    defineOption("","a5-global-input-file",PgOpt::STRING_OPT,
                 "Filename for the group in the global input file."
                 "[Experimental]",
                 "");
    defineOption("","a5-input-group",PgOpt::STRING_OPT,
                 "Group name for input. Defaults to --a5-group or if that "
                 "is not set, then '/'. "
                 "[Experimental]",
                 "");

    defineOption("","a5-sketch-input-global",PgOpt::BOOL_OPT,
                 "Read the sketch from the global A5 input file. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-sketch-input-file",PgOpt::STRING_OPT,
                 "Read the sketch from the an A5 input file. "
                 "[Experimental]",
                 "");
    defineOption("","a5-sketch-input-group",PgOpt::STRING_OPT,
                 "Group name to read the sketch from. Defaults to --a5-input-group. "
                 "[Experimental]",
                 "");
    defineOption("","a5-sketch-input-name",PgOpt::STRING_OPT,
                 "The name of the data section. Defaults to 'target-sketch'. "
                 "[Experimental]",
                 "");

    defineOption("","a5-feature-effects-input-global",PgOpt::BOOL_OPT,
                 "Read the feature effects global A5 input file. "
                 "[Experimental]",
                 "false");
    defineOption("","a5-feature-effects-input-file",PgOpt::STRING_OPT,
                 "Read the feature effects from the an A5 input file. "
                 "[Experimental]",
                 "");
    defineOption("","a5-feature-effects-input-group",PgOpt::STRING_OPT,
                 "Group name to read the feature effects from. Defaults to --a5-input-group. "
                 "[Experimental]",
                 "");
    defineOption("","a5-feature-effects-input-name",PgOpt::STRING_OPT,
                 "The name of the data section. Defaults to XXX.feature-response "
                 "where XXX is the analysis name and quant method. IE 'brlmm-p.plier'. "
                 "[Experimental]",
                 "");

    defineOptionSection("Engine Options (Not used on command line)");
    defOptMult("", "cels", PgOpt::STRING_OPT,
               "Cel files to process.",
               "");
    defOptMult("", "result-files", PgOpt::STRING_OPT,
               "CHP file names to output. Must be paired with cels.",
               "");
    defineOption("","annotation-file", PgOpt::STRING_OPT,
                 "NetAffx Annotation database file.",
                 "");
}

void ProbesetSummarizeEngine::defineStates() {
    defineOption("", "num-rows", PgOpt::INT_OPT,
                 "The number of rows on the chip.",
                 "-1");
    defineOption("", "num-cols", PgOpt::INT_OPT,
                 "The number of cols on the chip.",
                 "-1");
    defineOption("", "probe-count", PgOpt::INT_OPT,
                 "The number of probes on the chip.",
                 "-1");
    defineOption("", "probeset-count", PgOpt::INT_OPT,
                 "The number of probesets on the chip.",
                 "-1");
    defineOption("", "channel-count", PgOpt::INT_OPT,
                 "The number of data channels in the CEL file.",
                 "-1");
}

/** 
 * Make sure that our options are sane.
 * 
 * @param o - Options to be checked.
 */
void ProbesetSummarizeEngine::checkOptionsImp() {

    defineStates();

    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
    setLibFileOpt("pgf-file");
    setLibFileOpt("clf-file");
    setLibFileOpt("bgp-file");
    setLibFileOpt("probeset-ids");
    setLibFileOpt("meta-probesets");
    setLibFileOpt("probe-class-file");
    setLibFileOpt("qc-probesets");
    setLibFileOpt("use-feat-eff");
    setLibFileOpt("target-sketch");
    setLibFileOpt("reference-profile");
    setLibFileOpt("kill-list");
    setLibFileOpt("a5-global-input-file");
    setLibFileOpt("a5-sketch-input-file");
    setLibFileOpt("a5-feature-effects-input-file");

    if (getOpt("explain") != "") { explain(); exit(0); }

    if (getOpt("out-dir") == "") {Err::errAbort("Must specify an output directory.");}
    if (getOpt("temp-dir") == "") { setOpt("temp-dir", Fs::join(getOpt("out-dir"),"temp")); }

    bool cdfMode = false;
    bool pgfMode = false;

    if (getOpt("cdf-file")!="" || getOpt("spf-file")!="") 
        cdfMode = true;
    if (getOpt("pgf-file")!="")
        pgfMode = true;

    if (cdfMode && pgfMode) {
        Err::errAbort("Must specify only one layout spec: cdf, pgf/clf, or spf file.");
    }

    string cdfFile = getOpt("cdf-file");
    string spfFile = getOpt("spf-file");
    string pgfFile = getOpt("pgf-file");
    string clfFile = getOpt("clf-file");
    string bgpFile = getOpt("bgp-file");
    string killListFile = getOpt("kill-list");
    vector<string> analysis = getOptVector("analysis");

    // Check that a global file was given if the user wants global a5 output
    if (getOpt("a5-global-file") == "") {
        vector<string> optionNames;
        getOptionNames(optionNames);
        for (size_t i=0; i<optionNames.size(); i++) {
            string name = optionNames[i];
            if ( name.find("use-global") != name.npos ) {
                if ( getOptBool(name) ) {
                    Err::errAbort("--" + name + " option given, but no global file. Must specify --a5-global-file.");
                }
            }
        }
    }
    // Check that a global input file was given if the user wants global a5 input
    if (getOpt("a5-global-input-file") == "") {
        vector<string> optionNames;
        getOptionNames(optionNames);
        for (size_t i=0; i<optionNames.size(); i++) {
            string name = optionNames[i];
            if ( name.find("input-global") != name.npos ) {
                if ( getOptBool(name) ){
                    Err::errAbort("--" + name + " option given, but no global input file. Must specify --a5-global-input-file.");
                }
            }
        }
    }

    // make sure global input and output files are different
    if (getOpt("a5-global-input-file") == getOpt("a5-global-file"))
        if (getOpt("a5-global-input-file") != "")
            Err::errAbort("Unable to use the same A5 file for global input and output.");

    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if (celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);

    vector<string> resultFiles = getOptVector("result-files");
    if (resultFiles.size() > 0) {
        if (celFiles.size() != resultFiles.size())
            Err::errAbort("result-files option used but is not the same size as cel file listing");
        bool ccOut = getOptBool("cc-chp-output");
        bool mdOut = getOptBool("md-chp-output");
        bool xdaOut = getOptBool("xda-chp-output");
        if ((ccOut && mdOut) || (ccOut && xdaOut) || (mdOut && xdaOut)) 
            Err::errAbort("cannot use multiple CHP output formats with the result-files option");
    }

    if ((cdfMode) && (cdfFile=="" && spfFile=="")) {
        Err::errAbort("Must specify a layout file like cdf, spf, or pgf/clf files.");
    }
    if ((!cdfMode || pgfMode) && getOptBool("xda-chp-output")) {
        Err::errAbort("Must use a CDF file to generate an XDA CHP file.");
    }
    if (getOpt("probeset-ids") != "" && getOptBool("xda-chp-output")) {
        Err::errAbort("Can't specify a subset and XDA CHP output at same time. Must analyze all probesets for XDA CHP file output.");
    }
    if (killListFile != "" && !pgfMode && cdfFile == "") {
        Err::errAbort("Can't use a kill list with spf file");
    }
    if (killListFile != "" && getOptBool("xda-chp-output")) {
        Err::errAbort("Can't specify a kill list and XDA CHP output at same time.");
    }
    if (killListFile != "" && getOptBool("cc-md-chp-output")) {
        Err::errAbort("Can't specify a kill list and Expression CC CHP output at same time.");
    }
    if (killListFile != "" && getOptBool("cc-chp-output")) {
        Err::errAbort("Can't specify a kill list and Multi Data CC CHP output at same time.");
    }
    if (bgpFile!="" && !pgfMode)
        Err::errAbort("Can only use BGP file with PGF/CLF file. Not CDF or SPF.");

    /* Some sanity checks. */
    if (analysis.empty())
        Err::errAbort("Must specify at least one analysis to preform.");

    if (getOpt("probeset-ids") != "" && getOpt("meta-probesets") != "")
        Err::errAbort("Can't specify both probeset-ids and meta-probesets options.");

    // Check chip types
    vector<string> chipTypesInLayout;

    /* Get the intial info about the chip and check cel files to make sure
       they match. */
    colrow_t numRows = 0, numCols = 0;
    int probeCount = 0;
    int probeSetCount=0;
    int channelCount = 1;

    if (cdfFile!="")
        EngineUtil::getCdfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, cdfFile);
    else if (spfFile!="")
        EngineUtil::getSpfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, spfFile);
    else if (pgfMode) 
        EngineUtil::getPgfChipType(chipTypesInLayout, numRows, numCols, probeCount, pgfFile, clfFile);
    else
        Err::errAbort("Must specify a cdf file, spf file, or PGF and CLF files.");

    setOpt("num-rows", ToStr(numRows));
    setOpt("num-cols", ToStr(numCols));
    setOpt("probe-count", ToStr(probeCount));
    setOpt("probeset-count", ToStr(probeSetCount));
    setOpt("channel-count", ToStr(channelCount));

    if (chipTypesInLayout.empty() || chipTypesInLayout[0] == "" || probeCount == 0) 
        Err::errAbort("Problem determining ChipType in file: " + 
                      ( cdfFile != "" ? cdfFile : 
                        ( spfFile != "" ? spfFile : pgfFile + ToStr(", ") + clfFile) 
                          ) );

    /* Did the user "force" a set of chip types via options? */
    vector<string> chipTypesSupplied = getOptVector("chip-type");


    /* Figure out what chip type to report */
    if (chipTypesSupplied.size() > 0) {
        setOpt("chip-type", chipTypesSupplied[0]);
    } 
    else if (chipTypesInLayout.size() > 0) {
        setOpt("chip-type", chipTypesInLayout[0]);
    } 
    else {
        Err::errAbort("Unable to figure out a chip type.");
    }

    /* Do Chip Type Check */
    if (!getOptBool("force")) {
        if (chipTypesSupplied.size() > 0) {
            EngineUtil::checkCelChipTypes(chipTypesSupplied, probeCount, celFiles, numRows, numCols);
            EngineUtil::checkChipTypeVectors(chipTypesSupplied, chipTypesInLayout);
        }
        else {
            EngineUtil::checkCelChipTypes(chipTypesInLayout, probeCount, celFiles, numRows, numCols);
        }
    }
}


/** 
 * @brief High level function to set up analysis streams and then
 * start the gears turning.
 * 
 * @param asFactory - Object for making analysis streams.
 * @param errMsg - If there was a problem, error message gets filled in here.
 * @return true if successful, false otherwise.
 */
void ProbesetSummarizeEngine::runImp() {

    // probeset names act as a common pool for a number of different
    // things so we don't have to copy the relatively large (in RAM)
    // probe names over and over. They are known to be used for the chp
    // file writers and meta probesets, but may also be used for other
    // things as time goes on. Change or delete them prematurely at your
    // own risk.
    std::vector<const char *> probesetNames;
    std::vector<const char *> probesetDisplayNames;
    vector<bool> probeSubset;
    ChipLayout *layout = NULL;
    vector<MetaProbeset *> metaSets;
    vector<AnalysisStreamExpression *> analysisStreams;
    vector<int> probesetIds;
    ProbeListFactory plFactory;
    //vector<ProbeListPacked> probeListVec;
    vector<bool> psToLoad;
    /* Define the currently hard-coded p-value threshold for detection, used
       by QC reporters.  Comment in sample-probe-summarize.pl: this should be
       configurable from the commandline and the default set from pgf file. */
    unsigned int numProbes = 0;
    IntensityMart *iMart = NULL;
    CelReader reader;

    try { // outer try for handling the exception.
        try { // inner try for memory clean up.
            if (!Fs::isWriteableDir(getOpt("out-dir")))
                if (Fs::mkdirPath(getOpt("out-dir"), false) != APT_OK)
                    Err::errAbort("Can't make or write to directory: " + getOpt("out-dir"));
            makeTempDir(getOpt("temp-dir"));
  
            // Do we have a global output file?
            if (getOpt("a5-global-file") != "") {
                m_a5_global_output_file=affx::TsvReport::openA5File(getOpt("a5-global-file"), !getOptBool("a5-global-file-no-replace"));
                // default the group name
                string groupName;
                if (getOpt("a5-group")!="") {
                    groupName=getOpt("a5-group");
                } 
                else if (getOpt("set-analysis-name") != "") {
                    groupName=getOpt("set-analysis-name");
                } 
                else {
                    ///@todo we cannot call get name on analysis stream yet -- perhaps we need a static method to get the name given an analysis spec?
                    groupName="/";
                }
                m_a5_global_output_group=
                    m_a5_global_output_file->openGroup(groupName,affx::FILE5_CREATE|affx::FILE5_OPEN);
    
            }

            // ...a global input file?
            if (getOpt("a5-global-input-file") != "") {
                // is the input and output files the same?
                if (getOpt("a5-global-input-file") == getOpt("a5-global-file")) {
                    ///@todo allow same file for global input and output
                    Err::errAbort("Unable to use the same A5 file for global input and output.");
                } 
                else {
                    m_a5_global_input_file=new affx::File5_File();
                    m_a5_global_input_file->open(getOpt("a5-global-input-file"),affx::FILE5_OPEN_RO);
                }
                // default the group name
                string groupName;
                if (getOpt("a5-input-group") != "") {
                    groupName = getOpt("a5-input-group");
                } 
                else if (getOpt("set-analysis-name") != "") {
                    groupName=getOpt("set-analysis-name");
                } 
                else {
                    ///@todo we cannot call get name on analysis stream yet -- perhaps we need a static method to get the name given an analysis spec?
                    groupName="/";
                }
                m_a5_global_input_group=
                    m_a5_global_input_file->openGroup(groupName,affx::FILE5_OPEN_RO);
            }

            layout = new ChipLayout();
            ///@todo we only need this if pm-sum is used
            layout->setNeedPmAlleleMatch(true);

            // Load up some option values into local variables to make coding easier
            vector<string> celFiles = getOptVector("cels");
            vector<string> analysisStrings = getOptVector("analysis");
            unsigned int numRows = Convert::toUnsignedInt(getOpt("num-rows"));
            unsigned int numCols = Convert::toUnsignedInt(getOpt("num-cols"));

            // Load up probes to ignore
            probeidmap_t killList; // map of probe IDs to exclude
            string killListFile = getOpt("kill-list");
            if (killListFile != "") {
                ChipLayout::fillInKillList(killListFile, killList, numRows, numCols);
            }

            // If we have a subset or meta probeset list then read it.
            bool haveMetaProbeset = (getOpt("probeset-ids") != "" || getOpt("meta-probesets") != "");
            if (haveMetaProbeset) {
                bool metaFile = false;
                string fileName;
                if (getOpt("probeset-ids") != "" && getOpt("meta-probesets") != "")
                    Err::errAbort("Can't specify both probeset-ids and meta-probesets options.");
                if (getOpt("meta-probesets") != "") {
                    fileName = getOpt("meta-probesets");
                    metaFile = true;
                }
                else {
                    fileName = getOpt("probeset-ids");
                    metaFile = false;
                }
                MetaProbeset::readMetaProbesets(fileName, metaFile, metaSets);
            }

            // To satisfy the function signature we create an empty set of Probes.
            std::vector<Probe *> emptySetOfProbes;      

            // will fill in metaSets if empty
            loadChipLayout( *layout, 
                            metaSets, 
                            //layout->m_PlFactory.m_probelist_vec,
                            killList, 
                            &probesetNames,
                            (getOptBool("use-pgf-names") == true ? &probesetDisplayNames : NULL),
                            emptySetOfProbes, 
                            probeSubset, 
                            false, 
                            false);

            AnalysisStreamFactory asFactory(QuantMethodFactory::Expression);
            AnalysisStreamFactory::setupAnalysisStreamFactory(
                asFactory,
                m_a5_global_input_file,
                m_a5_global_input_group,
                m_a5_global_output_file,
                m_a5_global_output_group,
                getOptVector("cels").size(),
                getOptInt("probe-count"),
                getOpt("out-dir"),

                getOptBool("write-sketch"),
                getOpt("target-sketch"),

                getOptBool("write-profile"),
                getOpt("reference-profile"),

                getOptBool("a5-sketch"),
                getOptBool("a5-sketch-use-global"),

                getOptBool("a5-sketch-input-global"),
                getOpt("a5-sketch-input-file"),
                getOpt("a5-sketch-input-name"),
                getOpt("a5-sketch-input-group"),

                getOpt("use-feat-eff"),
                getOptBool("a5-feature-effects-input-global"),
                getOpt("a5-feature-effects-input-file"),
                getOpt("a5-feature-effects-input-name"),
                getOpt("a5-feature-effects-input-group"),

                getOpt("a5-input-group"),
                getOpt("a5-group"),

                getOpt("set-analysis-name"),
                ///@todo bit of a hack -- assumes a simple string search for med-polish or plier
                getOptVector("analysis")[0],

                getOpt("bgp-file"),
                getOpt("annotation-file"),
                *layout
                );


            int toRun = metaSets.size();

            mergeProbeSubset(layout->m_PlFactory.m_probelist_vec, probeSubset);

            vector<int> desiredOrder;
            determineDesiredOrder(*layout, desiredOrder, asFactory, haveMetaProbeset, metaSets);

            numProbes = layout->getProbeCount();  

            string diskDir = getOpt("temp-dir");
            if (getOptBool("use-disk")) {
                DiskIntensityMart* diskMart = 
                    new DiskIntensityMart(desiredOrder, 
                                          celFiles, 
                                          getOptInt("disk-cache") * 1048576,
                                          diskDir, 
                                          "apt-summarize.tmp",
                                          true);
                if (haveMetaProbeset || getOptBool("store-duplicate-probes")) {
                    diskMart->setUseAuxMemCache(true);
                }
                iMart = diskMart;
            }
            else {
                SparseMart* sparseMart = new SparseMart(desiredOrder, celFiles, true);
                iMart = sparseMart;
            }
        
            assert(celFiles.size() > 0);

            // set up cel stat listener
            map<string, vector<bool> > masks;

            vector<bool> pmMask = layout->getPmProbeMask();
            masks["pm"] = pmMask;

            vector<bool> mmMask = layout->getMmProbeMask();
            int mmCount = 0;
            for (size_t i=0; i<mmMask.size(); i++)
                if (mmMask[i])
                    mmCount++;
            if (mmCount > 5000) {
                masks["mm"] = mmMask;
            }

            // if there are gc control probes, save a background mask and
            // also add probes to probeSubset.
            const vector<Probe*> &probesToLoad = asFactory.getGcControlProbes(); 
            if (probesToLoad.size() > 0) {
                vector<bool> bgMask(pmMask.size(), false);
                //
                int probeSubset_size=probeSubset.size();
                int bgMask_size=bgMask.size();
                for (unsigned int i = 0; i < probesToLoad.size(); i++) {
                    int p_id = probesToLoad[i]->id;
                    if ((p_id>=probeSubset_size)||(p_id>=bgMask_size)) {
                        Err::errAbort("The probeid "+ToStr(p_id)+" from the BgpFile is out of range."
                                      " (0-"+ToStr(probeSubset_size)+")");
                    }
                    probeSubset.at(p_id) = true;
                    bgMask.at(p_id) = true;
                }
                masks["bgrd"] = bgMask;
            }

            Verbose::out(2, "Loading: " + ToStr(probeSubset.size()) + " probes");
            int subsetSum = 0;
            for (size_t blah = 0; blah < probeSubset.size(); blah++) {
                if (probeSubset[blah]) {
                    subsetSum++;
                }
            }

            if (getOpt("probe-class-file") != "")
                EngineUtil::readProbeClassFile(getOpt("probe-class-file"), layout->getProbeCount(), masks);

            map<string, vector<bool> >::iterator maskIter;
            for (maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
                int count = 0;
                for (vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++)
                    if (*mIx) 
                        count++;
                Verbose::out(2,"Number of " + maskIter->first + " probes is: " + ToStr(count));
            }

            CelStatListener celStats(masks);

            /* Make the analysis pathways. */
            for (size_t i = 0; i < analysisStrings.size(); i++) {
                AnalysisStreamExpression *as = NULL;

                if (analysisStrings[i] == "")
                    Err::errAbort("Invalid (blank) analysis specification");
                as = asFactory.constructExpressionAnalysisStream(analysisStrings[i], *layout, m_stdMethods, getOpt("set-analysis-name"));
                as->setOutPrefix(getOpt("out-dir"));
                QuantExprMethod *eMethod = NULL;
                eMethod = dynamic_cast<QuantExprMethod *>(as->getQuantMethod());
      
                AnalysisInfo analInfo = makeAnalysisInfo(as, *layout, metaSets);
                as->setInfo(analInfo);

                setupReporters(as, eMethod, metaSets, &celStats);
                analysisStreams.push_back(as);

                // %HACK% This is the least painful way I could think of to pass
                // the bgp filename to the QuantDabg QuantExprMethod.
                if (InstanceOf(eMethod, QuantDabg)) {
                    QuantDabg* dabgMethod = dynamic_cast<QuantDabg*>(eMethod);
                    dabgMethod->setBgpFileName(getOpt("bgp-file"));
                }
            }


            /* Set up the cel reader */
            reader.setFiles(celFiles);
            for (unsigned int analysisIx = 0; analysisIx < analysisStreams.size(); analysisIx++)
                analysisStreams[analysisIx]->registerChipStreamObjs(reader);
            reader.registerCelListener(&celStats);
            reader.registerIntensityMart(iMart);

            /* Read CEL files */
            reader.readFiles();

            writeHeaders(*iMart, analysisStreams);
            /* Write headers after iMart is populated as reports key off cels loaded into iMart */
            for (unsigned int analysisIx = 0; analysisIx < analysisStreams.size(); analysisIx++)
                analysisStreams[analysisIx]->prepare(*iMart);


            /* Do the analysis */
            doSummaries(*layout, *iMart, layout->m_PlFactory.m_probelist_vec, metaSets, analysisStreams, toRun);

            /* Let the streams do their post summaries cleanup. */
            Verbose::out(1,"Flushing output reporters. Finalizing output.");
            for (unsigned int analysisIx = 0; analysisIx < analysisStreams.size(); analysisIx++)
                analysisStreams[analysisIx]->finish();
            closeGlobalA5();

            Verbose::out(1,"Done.");
        } // inner try end
        catch(...) {
            for (uint32_t setIx = 0; setIx < metaSets.size(); setIx++) {
                delete metaSets[setIx];
            }
            metaSets.clear();
            for (uint32_t psIx = 0; psIx < probesetNames.size(); psIx++) {
                FreezArray(probesetNames[psIx]);
            }
            for (uint32_t psIx = 0; psIx < probesetDisplayNames.size(); psIx++) {
                FreezArray(probesetDisplayNames[psIx]);
            }
            for (uint32_t i = 0; i < analysisStreams.size(); i++) 
                delete analysisStreams[i];
            analysisStreams.clear();
            Freez(layout);
            Freez(iMart);
            closeGlobalA5();
            removeTempDir(getOpt("temp-dir"));

            // re-throw the unchanged exception.
            throw;
        }
    } // end outer try
    /* When things go wrong see if we can die gracefully here. */
    catch(Except &e) {
        Verbose::out(0,"");
        Err::errAbort(e.what());
    }
    catch(const std::bad_alloc &e) {
        Verbose::out(0,"");
        Err::errAbort("Ran out of memory. "
                      "Try using quitting other applications.");
    }
    catch(CalvinException &ce) {
        Verbose::out(0,"");
        Err::errAbort("Affymetrix GeneChip Command Console library has thrown an exception. "
                      "Calvin Error: '" +StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
    }
    catch(const std::exception &e) {
        Verbose::out(0,"");
        Err::errAbort("Exception caught. "
                      "Most problems with this program are memory related. "
                      "Message is: " + ToStr(e.what()));
    }
    catch(const BaseException &e) {
        Verbose::out(0,"");
        Err::errAbort("newmat issue: " + ToStr(e.what()));
    }
    catch(...) {
        Verbose::out(0,"");
        Err::errAbort("Unknown exception caught "
                      "(most problems with this program are memory related).");
    }
    // cleanup -- KEEP IN SYNC WITH MEMORY FREEZ ABOVE IN CATCH
    for (uint32_t setIx = 0; setIx < metaSets.size(); setIx++) {
        delete metaSets[setIx];
    }
    metaSets.clear();
    for (uint32_t psIx = 0; psIx < probesetNames.size(); psIx++) {
        FreezArray(probesetNames[psIx]);
    }
    for (uint32_t psIx = 0; psIx < probesetDisplayNames.size(); psIx++) {
        FreezArray(probesetDisplayNames[psIx]);
    }
    for (uint32_t i = 0; i < analysisStreams.size(); i++) 
        delete analysisStreams[i];
    analysisStreams.clear();
    Freez(layout);
    Freez(iMart);
    removeTempDir(getOpt("temp-dir"));
}

/** 
 * @brief Write the headers for the various reports.
 * 
 * @param layout - Probe and probe set info.
 * @param iMart - Memory reprentation of data.
 * @param analysis - Collection of all the analysis that we are going to run.
 * @param startTime - Time of the beginning of processing.
 */
void ProbesetSummarizeEngine::writeHeaders(IntensityMart &iMart,
                                           std::vector<AnalysisStreamExpression *> &analysis) {
  
    for (unsigned int analysisIx = 0; analysisIx < analysis.size(); analysisIx++) {
        string qMethodGuid = affxutil::Guid::GenerateNewGuid();
        AnalysisStreamExpression *as = analysis[analysisIx];
        assert(as);
        AnalysisInfo info = as->getInfo();
        as->addStdHeaders(getOpt("exec-guid"),
                          getOpt("time-start"),
                          getOpt("command-line"),
                          getOpt("version-to-report"),
                          info);
    }
}

ProbeSetGroup *ProbesetSummarizeEngine::makeProbesetGroupFromMeta(MetaProbeset &meta, ChipLayout &layout) {
    ProbeSetGroup *group = new ProbeSetGroup();
    group->ownMem = true;
    group->name = Util::cloneString(meta.name);
	group->displayName = meta.displayName;
    group->probeSets.resize(meta.probesets.size());
    int type = 0;
    int cnt = 0;
    for (int metaIx = 0; metaIx < meta.probesets.size(); metaIx++) {
        ProbeListPacked sub = layout.getProbeListByName(meta.probesets[metaIx]);
        if (sub.isNull()) {
            Verbose::out(2,"Skipping probeset " + ToStr(meta.probesets[metaIx]) + " in meta probeset " + ToStr(meta.name) + ". (Due to kill list?)");
        }
        else {
            if (cnt == 0) {
                type = sub.get_type();
            }
            if (type != sub.get_type()) {
                Err::errAbort("For meta probeset id: '" + ToStr(meta.name) + "'. Can't combine different types of probesets.");
            }
            ProbeSet *ps = ProbeListFactory::asProbeSet(sub);
            group->probeSets[cnt] = ps;
            cnt++;
        }
    }
    if (cnt > 0) {
        group->probeSets.resize(cnt);
        return group;
    } 
    else {
        Verbose::out(2,"Skipping meta probeset probeset " + ToStr(meta.name) + ". (Due to kill list?)");
        delete group;
        return NULL;
    }
}

/** 
 * @brief Loop through the groups and do the requested analysis for
 * each one.
 * 
 * @param layout - Probe and probe set info.
 * @param iMart - Memory reprentation of data.
 * @param psGroups - Groups of probe sets to be analyzed.
 * @param analysis - Analysis pathways to be computed.
 */
void ProbesetSummarizeEngine::doSummaries(ChipLayout &layout,
                                          IntensityMart &iMart,
                                          const vector<ProbeListPacked>& plVec,
                                          vector<MetaProbeset *> &metaToRun,
                                          vector<AnalysisStreamExpression *> &analysis,
                                          int numPsSets) {
    vector<QuantMethodReport *> qReports;
    if (metaToRun.size() == 0) {
        unsigned int dotMod = max((int)plVec.size()/20, 1);
        Verbose::progressBegin(1, ToStr("Processing Probesets"), 20, (int)dotMod, (int)plVec.size());
        /* Heavy lifting. */
        for (unsigned int groupIx = 0; groupIx < plVec.size(); groupIx++) {
            /* Little UI */
            Verbose::progressStep(1);
            /* Run each summary. */
            ProbeSet* ps = ProbeListFactory::asProbeSet(plVec[groupIx]);
            ProbeSetGroup psGroup(ps); 
            for (unsigned int analysisIx = 0; analysisIx < analysis.size(); analysisIx++) {
                analysis[analysisIx]->doAnalysis(psGroup, iMart, true);
            }
            // psGroup should delete the memory for ps...
        }
        Verbose::progressEnd(1, ToStr("Done."));
    }
    else {
        unsigned int dotMod = max((int)metaToRun.size()/20, 1);
        Verbose::progressBegin(1, ToStr("Processing Probesets"), 20, (int)dotMod, (int)metaToRun.size());
        /* Heavy lifting. */
        for (unsigned int groupIx = 0; groupIx < metaToRun.size(); groupIx++) {
            /* Little UI */
            Verbose::progressStep(1);
            ProbeSetGroup *psGroup = makeProbesetGroupFromMeta(*(metaToRun[groupIx]), layout);
            if (psGroup != NULL) {
                /* Run each summary. */
                for (unsigned int analysisIx = 0; analysisIx < analysis.size(); analysisIx++) {
                    analysis[analysisIx]->doAnalysis(*psGroup, iMart, true);
                }
                delete psGroup;
            }
        }
        Verbose::progressEnd(1, ToStr("Done."));
    }
}


/** 
 * Set all of the probes in the psGroups file to be loaded.
 * 
 * @param psGroups - Vector of probests containing probes to be set as true.
 * @param probeSubset - Bit vector in which probe ids should be set to true.
 */
void ProbesetSummarizeEngine::mergeProbeSubset(const std::vector<ProbeListPacked> &plVec, std::vector<bool> &probeSubset) {
    for (unsigned int plIx = 0; plIx < plVec.size(); plIx++) {
        for (int pIx = 0; pIx < plVec[plIx].probe_cnt(); pIx++) {
            int id = plVec[plIx].get_probeId(pIx);
            if (id != ProbeList::null_probe) {
                probeSubset[id] = true;
            }
        }
    }
}

/** 
 * Construct an object for output that contains the program an
 * algorithm information.
 * @param as - AnalysisStream to get information from.
 * @param layout - Layout to get information from.
 * @param metaSets - List of probeset names to output
 *
 * @return - Information for making CHP file.
 */
AnalysisInfo ProbesetSummarizeEngine::makeAnalysisInfo(AnalysisStream *as,
                                                       ChipLayout &layout, 
                                                       vector<MetaProbeset *> &metaSets) {

    AnalysisInfo info;
    info.m_NumExpression = layout.getNumExpressionPSets();
    info.m_NumGenotyping = layout.getNumGenotypeingPSets();
    info.m_NumRows = layout.getXCount();
    info.m_NumCols = layout.getYCount();
    info.m_NumProbeSets = metaSets.size();
    info.m_ProbeSetType = affxcdf::ExpressionProbeSetType;
    fillInAnalysisInfo(info, as);

    // If CHP Output, then we need list of probeset names to output
    // Otherwise do not waste the memory
    if (getOptBool("cc-md-chp-output") || getOptBool("cc-chp-output") || getOptBool("xda-chp-output")) {
        vector<const char *> probesetNames;
	vector<const char *> probesetDisplayNames;
        probesetNames.reserve(metaSets.size());
        vector<MetaProbeset *>::iterator mIx;
        for (mIx = metaSets.begin(); mIx != metaSets.end(); ++mIx) {
            probesetNames.push_back((*mIx)->name);
            probesetDisplayNames.push_back((*mIx)->displayName);
        }
        info.m_ProbesetNames = probesetNames;
	info.m_ProbesetDisplayNames = probesetDisplayNames;
    }

    // Sanity check
    Err::check(info.m_ParamValues.size() == info.m_ParamNames.size(), 
   
               "AnalysisInfo - Names and values out of sync.");
    return info;
}

void ProbesetSummarizeEngine::addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                            AnalysisInfo& info,
                            AnalysisStream *as, int precision) {

    QuantMethodExprReport *report = new QuantMethodExprReport(info.getNumCols());
    report->m_qMethod=as->getQuantMethod();
    report->setIsHeaderBuffer(1);
    report->setDirPath(getOpt("out-dir"));
    report->setPrecision(precision);
    report->m_summary_tsv_precision=precision;
    report->m_feffects_tsv_precision=precision-1;
    report->m_residuals_tsv_precision=precision;
    report->setFileprefix(as->getName());
    report->setFilename("QuantMethodExprReport-DEBUG"); // filename not used -- for debugging
    report->setWriteOldStyleFeatureEffectsFile(false);

    if (report_format == affx::TsvReport::FMT_TSV) {
        report->m_DoSummary=getOptBool("summaries");
        report->m_DoResiduals=getOptBool("feat-details");
        report->m_DoFeatureEffects=getOptBool("feat-effects");
    } 
    else {
        report->m_DoSummary=getOptBool("a5-summaries");
        report->m_DoResiduals=getOptBool("a5-feature-details");
        report->m_DoFeatureEffects=getOptBool("a5-feature-effects");

        if (getOptBool("a5-summaries-use-global")) {
            if (m_a5_global_output_group == NULL)
                Err::errAbort("--a5-summaries-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_summary_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for allele summary output.");
        }
    
        if (getOptBool("a5-feature-details-use-global")) {
            if (m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-details-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_residuals_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for residuals output.");
        }
    
        if (getOptBool("a5-feature-effects-use-global")) {
            if (m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-effects-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_feffects_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for feature effects output.");
        }
    }

    report->setFormat(report_format);

    report->addStdHeaders(report,
                          getOpt("exec-guid"),
                          as->getGuid(),
                          getOpt("time-start"),
                          getOpt("command-line"),
                          getOpt("version-to-report"),
                          info);
    as->addReporter(report);
}

/** 
 * Add the requested reporters for a given AnalysisStream.
 * 
 * @param as - AnalysisStream that needs the reporters.
 * @param eMethod - Expression quantification method that is being used to summarize probesets in this AnalysisStream.
 * @param metaSets - list of probesets to output.
 * @param stats -
 */
void ProbesetSummarizeEngine::setupReporters(AnalysisStreamExpression *as, 
                                             QuantExprMethod *eMethod, 
                                             vector<MetaProbeset *> &metaSets,
                                             CelStatListener *stats) {

    int precision = getOptInt("precision");
    if (precision > 16) {
        precision = 16;
        Verbose::out(2, "Changing precision for reports from "+ ToStr(precision) +" to 16");
    }
    string outDir = getOpt("out-dir");
    vector<string> celFiles = getOptVector("cels");

    AnalysisInfo info = as->getInfo();
  
    /* Add our primary text reporter. */
    if (precision < 0) {
        precision = 0;
        if (as->getQuantMethod()->getType() == "dabg") {
            precision = 5;
        }
    }
    // Add txt expression reporter
    if (getOptBool("summaries") || getOptBool("feat-effects") || getOptBool("feat-details")) {
        addReporters_Expression(affx::TsvReport::FMT_TSV,info,as,precision);
    }
    if (getOptBool("a5-summaries") || getOptBool("a5-feature-effects") || getOptBool("a5-feature-details")) {
        addReporters_Expression(affx::TsvReport::FMT_A5,info,as,precision);
    }

    bool isDetection = static_cast<QuantExprMethod *>(as->getQuantMethod())->getQuantType() == QuantExprMethod::Detection;

    QuantMethodExprChipSummary *runSummary = new QuantMethodExprChipSummary(celFiles, isDetection, 0, 0.01, getOpt("qc-probesets"), metaSets);
    runSummary->setFilename("runSummary");
  
    QuantMethodRunReport *runReport = new QuantMethodRunReport(celFiles);
    if (getOpt("report-file") != "") {
        string dir = Fs::dirname(getOpt("report-file"));
        if (dir == "")
            dir = outDir;
        string name = Fs::basename(getOpt("report-file"));
        runReport->setDirPath(dir);
        runSummary->setDirPath(dir);
        runReport->setFilename(name);
    } 
    else {
        runSummary->setDirPath(outDir);
        runReport->setDirPath(outDir);
        runReport->setFilename(as->getName()+".report");
    }
    runReport->setFormat(affx::TsvReport::FMT_TSV);
    runReport->setPrecision(6);

    runReport->registerChipSummary(dynamic_cast<ChipSummary *>(stats));
    runReport->registerChipSummary(dynamic_cast<ChipSummary *>(runSummary));

    set<string> group;
    set<uint32_t> groupId;
    as->addReporter(runSummary);
    as->addReporter(runReport);

    // Setup Multi Data CHP Output
    if (getOptBool("cc-md-chp-output")){
        string ccchpDir;
        /* If chp directory not specified use  output/cc-chp/ by default. */
        if (getOpt("cc-md-chp-out-dir") == "") { ccchpDir = Fs::join(outDir,"cc-md-chp"); }
        else { ccchpDir = getOpt("cc-md-chp-out-dir"); }
        Util::chompLastIfSep(ccchpDir);
        QuantMethodMultiDataCCCHPReport *ccchpReporter = new QuantMethodMultiDataCCCHPReport(ccchpDir);
        vector<string> resultFiles = getOptVector("result-files");
        if (resultFiles.size() > 0)
            ccchpReporter->setChpFileNames(resultFiles);
        ccchpReporter->SetGenoTypeOnly(false);
        ccchpReporter->registerChipSummary(dynamic_cast<ChipSummary *>(stats));
        ccchpReporter->registerChipSummary(dynamic_cast<ChipSummary *>(runSummary));
        as->addReporter(ccchpReporter);
    }

    // Setup Old Expression CHP Output
    if (getOptBool("cc-chp-output")) {
        string ccchpDir;
        /* If chp directory not specified use  output/cc-chp/ by default. */
        if (getOpt("cc-chp-out-dir") == "") { ccchpDir = Fs::join(outDir,"cc-chp"); }
        else { ccchpDir = getOpt("cc-chp-out-dir"); }
        Util::chompLastIfSep(ccchpDir);
        QuantMethodExprCCCHPReport *ccchpReporter = new QuantMethodExprCCCHPReport(info, ccchpDir, eMethod->getType());
        vector<string> resultFiles = getOptVector("result-files");
        if (resultFiles.size() > 0)
            ccchpReporter->setChpFileNames(resultFiles);
        ccchpReporter->registerChipSummary(dynamic_cast<ChipSummary *>(stats));
        ccchpReporter->registerChipSummary(dynamic_cast<ChipSummary *>(runSummary));
        as->addReporter(ccchpReporter);
    }

    // setup XDA CHP output if neccesary.
    if (getOptBool("xda-chp-output")) {
        string chpDir;
        /* If chp directory not specified use  output/chp/ by default. */
        if (getOpt("xda-chp-out-dir") == "") { chpDir =  Fs::join(outDir,"chp"); } 
        else { chpDir = getOpt("xda-chp-out-dir"); }
        Util::chompLastIfSep(chpDir);
        QuantMethodExprCHPReport *chpReporter = new QuantMethodExprCHPReport(info, chpDir, eMethod->getType());
        vector<string> resultFiles = getOptVector("result-files");
        if (resultFiles.size() > 0)
            chpReporter->setChpFileNames(resultFiles);
        as->addReporter(chpReporter);
    }
  
    if (getOptBool("subsample-report")) {
        QuantMethodSamplerReport *sReport = NULL;
        /* Add our sampling reporter, just reports on 1 out of every 10 probesets. */
        string prefix = Fs::join(outDir,as->getName());
        sReport = new QuantMethodSamplerReport(9, prefix, 5, true, true, true);
        as->addReporter(sReport);
    }

}

/** 
 * Read the design of the chip from the cdf file. Also can read the
 * design from an spf file is that is supplied.
 * 
 * @param layout - Object to be filled in.
 * @param psGroups - Probe sets to be filled in.
 * @param probesToLoad - Additional probes to be loaded (i.e. background probes.)
 * @param probeSubset - Bitmap to remember the probes that were actually loaded.
 */
void ProbesetSummarizeEngine::loadCdfLayout(ChipLayout &layout, 
                                            vector<MetaProbeset *> &metaToLoad,
                                            const vector<Probe *> &probesToLoad, 
                                            std::vector<const char *> *probesetNames,
                                            vector<bool> &probeSubset, 
                                            bool justStats, 
                                            bool doAll,
                                            probeidmap_t &killList)
{
    set<const char *, Util::ltstr> psNameLoadMap;
    string cdfFile = getOpt("cdf-file");
    string spfFile = getOpt("spf-file");

    // Are we doing a subset of all the probe sets? 
    if (!metaToLoad.empty() && !doAll) {
        for (unsigned int i = 0; i < metaToLoad.size(); i++) {
            for (unsigned int psIx = 0; psIx < metaToLoad[i]->probesets.size(); psIx++) { 
                psNameLoadMap.insert(metaToLoad[i]->probesets[psIx]);
            }
        }
    }
    // Load up the data from the file.
    string fileName = cdfFile!="" ? cdfFile : spfFile;
    try {
        if (cdfFile!="") {
            Verbose::out(1, ToStr("Opening cdf file: ") + Fs::basename(cdfFile));
            if (!layout.openCdf(cdfFile, 
                                psNameLoadMap, 
                                probesetNames,
                                probeSubset, 
                                "", 
                                killList, 
                                justStats)) 
                Err::errAbort("Couldn't open cdf file: " + cdfFile);
        }
        else if (spfFile!="") {
            std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
            if (getOpt("kill-list") != "")
                Err::errAbort("Cannot specify probe kill list with spf file");
            Verbose::out(1, ToStr("Opening layout file: ") + Fs::basename(spfFile));
            layout.openSpf(spfFile, 
                           psNameLoadMap, 
                           probesetNames, 
                           probeSubset,
                           "", 
                           justStats, 
                           psTypesToLoad);
        }
        else {
            Err::errAbort("Must have either a cdf file or spf file in cdf mode.");
        }
    }
    catch(const Except &e) {
        Err::errAbort(e.what());
    }
    catch(...) {
        Err::errAbort("Uncaught exception opening file: '" + fileName + "'");
    }
}

/** 
 * Read the design of the chip from the pgf & clf files.
 * 
 * @param layout - Object to be filled in.
 * @param metaToLoad
 * @param probesToLoad - Additional probes to be loaded (i.e. background probes.)
 * @param probeSubset - Bitmap to remember the probes that were actually loaded.
 * @param killList  
 */
void ProbesetSummarizeEngine::loadPgfLayout(ChipLayout &layout, 
                                            std::vector<MetaProbeset *> &metaToLoad,
                                            const std::vector<Probe *> &probesToLoad, 
                                            std::vector<const char *> *probesetNames,
                                            std::vector<const char *> *probesetDisplayNames,
                                            std::vector<bool> &probeSubset, 
                                            probeidmap_t &killList,
                                            bool justStats, bool doAll) {
    affx::ClfFile clf;
    vector<string> requiredCols;
    unsigned int i = 0;
    unsigned int numProbes = 0;
    std::set<const char *, Util::ltstr> psIdLoadMap;
    vector<unsigned int> psIdNotInGroup;
    std::string pgfFile = getOpt("pgf-file"); 
    std::string clfFile = getOpt("clf-file");

    /* open clf file. */
    if (clfFile=="")
        Err::errAbort("Must specify a clf-file for chip.");
    Verbose::out(1, ToStr("Opening clf file: ") +  Fs::basename(clfFile));
    if (!clf.open(clfFile))
        Err::errAbort(ToStr("Couldn't open clf file: ") + clfFile);
    layout.setDimensions(clf.getXMax() + 1, clf.getYMax() + 1);

    /* Make sure our assumptions about CLF file are true */
    Err::check(clf.getSequential() == 1,
               "ProbesetSummarizeEngine::loadPgfLayout() - unable to handle clf file without sequential set to 1.");
    Err::check(clf.getOrder().compare("col_major") == 0 || clf.getOrder().compare("row_major") == 0,
               "Unable to handle clf file without order set to row_major (old mislabeled 'col_major' accepted due to earlier bug.)");

    /* What probes are we loading? */
    numProbes = (clf.getXMax() + 1) * (clf.getYMax() + 1);
    probeSubset.resize(numProbes);
    for (i = 0; i < probeSubset.size(); i++) 
        probeSubset[i] = false;

    // Are we doing a subset of all the probe sets? 
    if (!metaToLoad.empty() && !doAll) {
        for (i = 0; i < metaToLoad.size(); i++) {
            for (unsigned int psIx = 0; psIx < metaToLoad[i]->probesets.size(); psIx++) { 
                psIdLoadMap.insert(metaToLoad[i]->probesets[psIx]);
            }
        }
    }
  
    /* Open actual pgf file. */
    if (Util::sameString(pgfFile,"")) 
        Err::errAbort("Must specify PGF for chip.");
    Verbose::out(1, ToStr("Opening pgf file: ") + Fs::basename(pgfFile));
    time_t startTime = time(NULL);
    if (!layout.openPgf(pgfFile, clf.getXMax() + 1, clf.getYMax() + 1, psIdLoadMap, probesetNames, probesetDisplayNames,
                        probeSubset, "", killList, justStats, false)) {
        Err::errAbort("Couldn't open PGF file: " + ToStr(pgfFile));
    }
    time_t endTime = time(NULL);
    int t = endTime - startTime;
    Verbose::out(2, ToStr("Pgf load took approximately: ") + ToStr(t) + ToStr(" seconds."));
}


/** 
 * Load up the probe set definitions.
 * 
 * @param layout - Object to be filled in.
 * @param psGroups - Probe sets to be filled in.
 * @param metaToLoad
 * @param killList  
 * @param probesToLoad - Additional probes to be loaded (i.e. background probes.)
 * @param probeSubset - Bitmap to remember the probes that were actually loaded.
 * @param probetIds 
 * @param probesetNames
 */
void ProbesetSummarizeEngine::loadChipLayout(ChipLayout &layout, 
                                             vector<MetaProbeset *> &metaToLoad,  
                                             probeidmap_t &killList,
                                             std::vector<const char *> *probesetNames,
                                             std::vector<const char *> *probesetDisplayNames,
                                             const vector<Probe *> &probesToLoad, 
                                             vector<bool> &probeSubset, 
                                             bool justStats, 
                                             bool doAll)
{
    if (getOpt("cdf-file")!="" || getOpt("spf-file")!="") 
        loadCdfLayout(layout,
                      metaToLoad, 
                      probesToLoad, 
                      probesetNames, 
                      probeSubset, 
                      justStats, 
                      doAll, 
                      killList);
    else if (getOpt("pgf-file")!="")
        loadPgfLayout(layout, 
                      metaToLoad, 
                      probesToLoad, 
                      probesetNames,
                      probesetDisplayNames,
                      probeSubset, 
                      killList, 
                      justStats, 
                      doAll);
    else
        Err::errAbort("Must be in either pgf or cdf mode.");
    if (metaToLoad.empty()) {
        metaToLoad.reserve(probesetNames->size());
        for (unsigned int psIx = 0; psIx < probesetNames->size(); psIx++) {
            MetaProbeset *set = new MetaProbeset();
            set->freeSubs = false;
            set->name = (*probesetNames)[psIx];
            set->displayName = (probesetDisplayNames != NULL ? (*probesetDisplayNames)[psIx] : NULL);
            set->probesets.resize(1);
            set->probesets[0] = set->name;
            metaToLoad.push_back(set);
        }
    }
  else if (probesetDisplayNames != NULL)
  {
	  map<string, MetaProbeset*> mpsmap;
	  for(unsigned int mIx = 0; mIx < metaToLoad.size(); mIx++) {
		  mpsmap[metaToLoad[mIx]->name] = metaToLoad[mIx];
	  }
	  for(unsigned int psIx = 0; psIx < probesetNames->size(); psIx++) {
		  const char *name = (*probesetNames)[psIx];
		  if (mpsmap.find(name) != mpsmap.end()) {
			  mpsmap[name]->displayName = (*probesetDisplayNames)[psIx];
		  }
	  }
  }
//  unsigned int psCount = layout.getProbeSetCount();
//  plVec.resize(psCount);
//  for (unsigned int psIx = 0; psIx < psCount; psIx++) {
//    plVec[psIx] = layout.getProbeListAtIndex(psIx);
//  }
}

void
ProbesetSummarizeEngine::determineDesiredOrder(ChipLayout& layout,
                                               vector<int> &desiredOrder, 
                                               AnalysisStreamFactory &asFactory,
                                               bool haveSubset,
                                               vector<MetaProbeset *> &metaSets) {
    const vector<Probe *> &gcProbes = asFactory.getGcControlProbes();
    // If no meta probesets then just use the layout in the cdf/spf/pgf

    if (!haveSubset) {
        const vector<int> &layoutOrder = layout.getProbeLayoutOrder();
        desiredOrder.reserve(gcProbes.size() + layoutOrder.size());
        for (int i = 0; i < layoutOrder.size(); i++) {
            desiredOrder.push_back(layoutOrder[i]);
        }
    }
    // If meta sets then we want to order in the same order as the
    // meta probesets.
    else {
        for (int i = 0; i < metaSets.size(); i++) {
            for (int j = 0; j < metaSets[i]->probesets.size(); j++) {
                ProbeListPacked pl = layout.getProbeListByName(metaSets[i]->probesets[j]);
                if (!pl.isNull()) {
                    for (int pIx = 0; pIx < pl.probe_cnt(); pIx++) {
                        desiredOrder.push_back(pl.get_probeId(pIx));
                    }
                }
            }
        }
    }
    for (int i = 0; i < gcProbes.size(); i++) {
      desiredOrder.push_back(gcProbes[i]->id);
    }
}

void ProbesetSummarizeEngine::addProbeSet(ProbeSet *basicPs, ProbeSet *currentPs) {
    if (basicPs->atomsPerGroup.size() != 1) {
        Err::errAbort("addProbeset() - Error: Can't merge probesets with more than one block");
    }
    for (int i = 0; i < currentPs->atoms.size(); i++) {
        basicPs->atoms.push_back(currentPs->atoms[i]);
        ///@todo we should really make sure that PM-only and PM-MM probesets are merged correctly
            if (basicPs->hasMM() != currentPs->hasMM())
                Err::errAbort("Cannot merge probesets which have MM with those that do not.");
            basicPs->atomsPerGroup[0] += (currentPs->hasMM() ? 
                                          currentPs->atoms[i]->probes.size() / 2 : currentPs->atoms[i]->probes.size());
    }

}

void ProbesetSummarizeEngine::fillInAnalysisInfo(AnalysisInfo &info, AnalysisStream *as, std::string prefix) {
    ///@todo need to refactor to ensure consistency between meta info tracked in probeset-summarize
    ///      versus probeset-genotype
    ///@todo automatically pull state and options info?
    ///@todo sync up option name with label in chp header
    // FOR NOW -- PLEASE KEEP THIS IN SYNC WITH probeset-genotype
    vector<string> celFiles = getOptVector("cels");
    assert(as);
    QuantMethod *qMethod = as->getQuantMethod();
    Err::check(qMethod != NULL, "Options::fillInAnalysisInfo() - Must have quantification method.");
    info.m_AlgVersion = qMethod->getVersion();
    info.m_AlgName = as->getName();
    info.m_ProgramName = getOpt("program-name");
    info.m_ProgramVersion = getOpt("version-to-report");
    info.m_ProgramCompany = getOpt("program-company");
    info.m_ChipType = getOpt("chip-type");
    info.m_ProgID = "";
    info.m_ExecGuid = getOpt("exec-guid");
    info.m_AnalysisGuid = as->getGuid();

    // State and Execution info
    info.addParam("apt-engine", "ProbesetSummarizeEngine");
    info.addParam(prefix + "program-name", getOpt("program-name"));
    info.addParam(prefix + "command-line", getOpt("command-line"));
    info.addParam(prefix + "exec-guid", getOpt("exec-guid"));
    info.addParam(prefix + "analysis-guid", as->getGuid());
    info.addParam(prefix + "time-str", getOpt("time-start"));
    info.addParam(prefix + "version", getOpt("version-to-report"));
    info.addParam(prefix + "cvs-id", getOpt("program-cvs-id"));
    info.addParam(prefix + "free-mem", getOpt("free-mem-at-start"));
    ///@todo add other state:
    ///@todo time-end state

    // Engine Options
    std::string subprefix = "opt-";
    ///@todo probably remove do-ps-names -- no longer an option
    info.addParam(prefix + subprefix + "do-ps-names",   (getOpt("cdf-file") != "" || getOpt("spf-file")!="" ? "true" : "false"));
    info.addParam(prefix + subprefix + "chip-type", getOpt("chip-type"));
    info.addParam(prefix + subprefix + "probe-count", getOpt("probe-count"));
    info.addParam(prefix + subprefix + "force", getOpt("force"));
    info.addParam(prefix + subprefix + "precision", getOpt("precision"));
    info.addParam(prefix + subprefix + "out-dir", getOpt("out-dir"));
    info.addParam(prefix + subprefix + "cc-expr-chp-out-dir", getOpt("cc-chp-out-dir"));
    info.addParam(prefix + subprefix + "cc-md-chp-out-dir", getOpt("cc-md-chp-out-dir"));
    info.addParam(prefix + subprefix + "xda-chp-out-dir", getOpt("xda-chp-out-dir"));
    info.addParam(prefix + subprefix + "cdf-file",       Fs::basename(getOpt("cdf-file")));
    info.addParam(prefix + subprefix + "spf-file",       Fs::basename(getOpt("spf-file")));
    info.addParam(prefix + subprefix + "pgf-file",       Fs::basename(getOpt("pgf-file")));
    info.addParam(prefix + subprefix + "clf-file",       Fs::basename(getOpt("clf-file")));
    info.addParam(prefix + subprefix + "bgp-file",       Fs::basename(getOpt("bgp-file")));
    info.addParam(prefix + subprefix + "ps-list-file",   Fs::basename(getOpt("probeset-ids")));
    info.addParam(prefix + subprefix + "meta-ps-file",   Fs::basename(getOpt("meta-probesets")));
    info.addParam(prefix + subprefix + "qc-groups-file", Fs::basename(getOpt("qc-probesets")));
    info.addParam(prefix + subprefix + "kill-list",      Fs::basename(getOpt("kill-list")));
    info.addParam(prefix + subprefix + "temp-dir", getOpt("temp-dir"));
    info.addParam(prefix + subprefix + "use-disk", getOpt("use-disk"));
    info.addParam(prefix + subprefix + "diskCache", getOpt("disk-cache"));
    info.addParam(prefix + subprefix + "cel-count", ToStr(celFiles.size()));
    for (uint32_t i = 0; i < celFiles.size(); i++) {
      std::string paramName = prefix + subprefix + "cel-" + ToStr(i+1);
      info.addParam(paramName, Fs::basename(celFiles[i]));
    }

    // Algorithm Options
    info.addParam(prefix + subprefix + "analysis-name", as->getName());
    info.addParam(prefix + subprefix + "set-analysis-name", getOpt("set-analysis-name"));
    info.addParam(prefix + subprefix + "analysis-spec", as->getAnalysisSpec());
    info.addParam(prefix + subprefix + "feat-effect-file", Fs::basename(getOpt("use-feat-eff")));
    info.addParam(prefix + subprefix + "target-sketch-file", Fs::basename(getOpt("target-sketch")));
    info.addParam(prefix + subprefix + "do-residuals", getOpt("feat-details"));
    info.addParam(prefix + subprefix + "do-feature-effects", getOpt("feat-effects"));
    info.addParam(prefix + subprefix + "write-sketch", getOpt("write-sketch"));
    info.addParam(prefix + subprefix + "reference-profile-file",Fs::basename(getOpt("reference-profile")));
    info.addParam(prefix + subprefix + "write-profile", getOpt("write-profile"));

    // no set-gender-method
    // no read-priors-brlmm
    // no initial-calls
    // no output-summaries
    // no list-sample
    // no write-prior
    // no chrX-file
    // no chrX-snps
    // no no-gender-force
    // no special-snps
    // no compute-prior
    // no norm-size
    // no dm-thresh
    // no dm-het-mult
    // no read-genders
    // no model-file-brlmm
    // no model-file-brlmmp
    // no model-file-birdseed
    // no genotypes-file
    // no dm-calls-out
    // no qmethod-spec
    // no genotype-only

    // Quantification Related Meta Info
    info.addParam("quantification-name", qMethod->getType());
    info.addParam("quantification-version", qMethod->getVersion());
    info.addParam("quantification-scale", QuantMethod::scaleToTxt(qMethod->getScale()));
    info.addParam("quantification-type", QuantMethod::quantTypeToTxt(qMethod->getQuantType()));
}

void ProbesetSummarizeEngine::printStandardMethods(std::ostream &out) {
    map<string,string>::iterator iter;
    out << endl << "Standard Methods:" << endl;
    unsigned int maxLength = 0;
    for (iter = m_stdMethods.begin(); iter != m_stdMethods.end(); iter++) {
        if (iter->first.size() > maxLength)
            maxLength = iter->first.size();
    }
    for (iter = m_stdMethods.begin(); iter != m_stdMethods.end(); iter++) {
        unsigned int currentLength = 0;
        out << " '" << iter->first << "' ";
        currentLength = iter->first.size();
        while (currentLength < maxLength + 1) {
            out.put(' ');
            currentLength++;
        }
        out << iter->second << endl;
    }
}

void ProbesetSummarizeEngine::closeGlobalA5() {
    // A5 close global groups
    if (m_a5_global_output_group!=NULL) {
        m_a5_global_output_group->close();
        delete m_a5_global_output_group;
        m_a5_global_output_group=NULL;
    }
    if (m_a5_global_input_group!=NULL) {
        m_a5_global_input_group->close();
        delete m_a5_global_input_group;
        m_a5_global_input_group=NULL;
    }

    // A5 close global files
    if (m_a5_global_output_file != m_a5_global_input_file) {
        affx::TsvReport::closeA5File(m_a5_global_output_file);
        affx::TsvReport::closeA5File(m_a5_global_input_file);
    } 
    else {
        affx::TsvReport::closeA5File(m_a5_global_output_file);
        m_a5_global_input_file = NULL;
    }
}

void ProbesetSummarizeEngine::extraHelp() {
    ChipStreamFactory cFactory;
    PmAdjusterFactory aFactory;
    QuantMethodFactory qFactory(QuantMethodFactory::Expression);
    AnalysisStreamFactory asFactory;

    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", cFactory.getDocs());
    EngineUtil::printSelfDocs("Pm Intensity Adjustments:", aFactory.getDocs());
    EngineUtil::printSelfDocs("Quantification Methods:", qFactory.getDocs());
    EngineUtil::printSelfDocs("Analysis Streams:", asFactory.getDocs());
}

void ProbesetSummarizeEngine::explain() {
    EngineUtil::explainParameter(getOpt("explain"));
}
