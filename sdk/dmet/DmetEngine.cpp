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
 * @file   DmetEngine.cpp
 * @author Alan Williams
 * @date   Mon Jun 23 14:57:34 PDT 2008
 * 
 * @brief Analysis engine for DMET 3.0
 */

//
#include "dmet/DmetEngine.h"
//
#include "dmet/DmetCHPWriter.h"
#include "dmet/DmetCopyNumberEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-probeset-genotype/ProbesetGenotypeEngine.h"
#include "chipstream/apt-probeset-summarize/ProbesetSummarizeEngine.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
//
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

DmetEngine::Reg DmetEngine::reg;

DmetEngine * DmetEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == DmetEngine::EngineName())
		return (DmetEngine *)engine;
	return NULL;
}

/**
 * Constructor
 */
DmetEngine::DmetEngine() {
    defineOptions();

    //
    m_ArgvPosAPS = -1;
    m_ArgvPosCN = -1;
    m_ArgvPosAPG = -1;
}

int DmetEngine::parseArgv( const char * const * const argv, int start ){
    vector<string> argvStrings;
    for (const char* const * arg=argv;*arg!=NULL;arg++) {
        argvStrings.push_back(*arg);
    }
    int argc = argvStrings.size();

    // Parse DmetEngine Options
    int argvPos = Options::parseArgv(argv, start);

    vector<string> celFiles;
    for(vector<const char *>::size_type i = 0; i < getArgCount(); i++)
        celFiles.push_back(getArg(i));
    setOpt("cels",celFiles);

    // Allow user to override APS defaults
    if(argc > argvPos) {
        ProbesetSummarizeEngine pse;
        int newArgvPos = pse.parseArgv(argv, argvPos+1);
        Verbose::out(1,"Parsed " + ToStr(newArgvPos - argvPos - 1) + " extra options for ProbesetSummarizeEngine!");
        m_ArgvPosAPS = argvPos + 1;
        argvPos = newArgvPos;
    } else {
        m_ArgvPosAPS = -1;
    }

    // Allow user to override CN defaults
    if(argc > argvPos) {
        DmetCopyNumberEngine cde;
        int newArgvPos = cde.parseArgv(argv,argvPos+1);
        Verbose::out(1,"Parsed " + ToStr(newArgvPos - argvPos - 1) + " extra options for DmetCopyNumberEngine!");
        m_ArgvPosCN = argvPos + 1;
        argvPos = newArgvPos;
    } else {
        m_ArgvPosCN = -1;
    }
    
    // Allow user to override APG defaults
    if(argc > argvPos) {
        ProbesetGenotypeEngine pge;
        int newArgvPos = pge.parseArgv(argv,argvPos+1);
        Verbose::out(1,"Parsed " + ToStr(newArgvPos - argvPos - 1) + " extra options for ProbesetGenotypeEngine!");
        m_ArgvPosAPG = argvPos + 1;
        argvPos = newArgvPos;
    } else {
        m_ArgvPosAPG = -1;
    }

    m_argv = argv;
    return argvPos;
}

/**
 * Destructor
 */
DmetEngine::~DmetEngine() {
}

void DmetEngine::defineOptions() {

    defineOptionSection("Input Options");
        defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
        defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets. Use either --cdf-file or --spf-file. ",
                     "");
        defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "File defining probe sets in spf (simple probe format) which is like a text cdf file.",
                     "");
        defineOption("", "special-snps",PgOpt::STRING_OPT,
                    "File containing all snps of unusual copy (chrX,mito,Y)",
                    "");
        defineOption("", "chrX-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrX. "
                     "Used for copy number probe chrX/Y ratio gender calling. ",
                     "");
        defineOption("", "chrY-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrY. "
                     "Used for copy number probe chrX/Y ratio gender calling. ",
                     "");
        defineOption("", "reference-input", PgOpt::STRING_OPT,
                     "Reference file with cluster prior information. ",
                     "");
        defineOption("s", "probeset-ids", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to analyze.",
                     "");
        defineOption("", "probeset-ids-reported", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to report. "
                     "This should be a subset of those specified with --probeset-ids if that option is used.",
                     "");
        defOptMult("", "chip-type", PgOpt::STRING_OPT,
                     "Chip types to check library and CEL files against. "
                     "Can be specified multiple times. "
                     "The first one is propigated as the chip type in the output files. "
                     "Warning, use of this option will override the usual check between chip types "
                     "found in the library files and cel files. You should use this option instead "
                     "of --force when possible. ",
                     "");
        defineOption("", "region-model", PgOpt::STRING_OPT,
                     "Regions model parameters. ",
                     "");
        defineOption("", "probeset-model", PgOpt::STRING_OPT,
                     "Probeset model parameters. ",
                     "");
	defineOption("", "cn-region-gt-probeset-file", 
			       PgOpt::STRING_OPT,
			       "Tab delimited file mapping probeset ids to copynumber regions. ",
			       "");

    defineOptionSection("Output Options");
        defineOption("", "cc-chp-output", PgOpt::BOOL_OPT,
                    "Output resulting calls in directory called 'chp' under out-dir. "
                    "This makes one AGCC Multi Data CHP file per cel file analyzed.",
                    "false");
        defineOption("", "reference-output", PgOpt::STRING_OPT,
                     "File to write reference values to. Specifying this option will turn on "
                     "dynamic clustering. WARNING: Currently the "
                     "resulting reference file is not really usable as a reference. "
                     "See the manual for more info. ",
                     "");
        defineOption("", "batch-name", PgOpt::STRING_OPT,
                     "The name of the batch for the dynamic cluster analysis. ",
                     "");

    defineOptionSection("Analysis Options");
        defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
                  "Explicitly set the analysis name. This affects output file names (ie prefix) and various meta info.",
                    "dmet");
        defineOption("", "ps-analysis", PgOpt::STRING_OPT,
                  "Explicitly set the ProbesetSummarizeEngine analysis string.",
                  "");
        defineOption("", "gt-analysis", PgOpt::STRING_OPT,
                  "Explicitly set the ProbesetGenotypeEngine analysis string.",
                  "");
        defineOption("", "gt-qmethod-spec", PgOpt::STRING_OPT,
                  "Explicitly set the ProbesetGenotypeEngine quant spec.",
                  "");
        defineOption("", "sample-type", PgOpt::STRING_OPT,
                    "Set the type of samples being processed. eg genomic, plasmid.",
                    "unknown");
        defineOption("", "batch-info", PgOpt::BOOL_OPT,
                    "Indicates whether or not information about other cel files in the batch "
                    "should be reported in CHP headers.",
                    "false");
        defineOption("", "null-context", PgOpt::BOOL_OPT,
                    "Indicates whether or not context info should be populated in the CHP files.",
                    "true");
        defineOption("", "run-cn-engine", PgOpt::BOOL_OPT,
                  "Indicates if the CN engine should be run or not.",
                    "true");
        defineOption("", "pra-thresh", PgOpt::INT_OPT,
                  "The threshold for calling PRAs based on the cluster mean strength.",
                  "3");
        defineOption("", "geno-call-thresh", PgOpt::DOUBLE_OPT,
                  "The confidence threshold for reporting calls in the CHP file.",
                  "0.1");

    defineOptionSection("Gender Options");

        defineOption("", "female-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio method.",
                    "0.17");
        defineOption("", "male-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio method.",
                    "0.68");

    defineOptionSection("Advanced Options");
        defineOption("", "call-coder-max-alleles", PgOpt::INT_OPT,
                    "For encoding/decoding calls, the max number of alleles per marker to allow.",
                    "6");
        defineOption("", "call-coder-type", PgOpt::STRING_OPT,
                    "The data size used to encode the call.",
                    "UCHAR");
        defineOption("", "call-coder-version", PgOpt::STRING_OPT,
                    "The version of the encoder/decoder to use",
                    "1.0");

    defineOptionSection("Execution Control Options");
        defineOption("", "use-disk", PgOpt::BOOL_OPT,
                     "Store CEL intensities to be analyzed on disk.", "true");
        defineOption("", "disk-cache", PgOpt::INT_OPT,
                     "Size of intensity memory cache in millions of intensities (when --use-disk=true).",
                     "50");

    defineOptionSection("Engine Options (Not used on command line)");
        defOptMult("", "cels", PgOpt::STRING_OPT,
                     "Cel files to process.",
                     "");
        defOptMult("", "report", PgOpt::STRING_OPT,
                     "Probesets to report. eg consented.",
                     "");
}

void DmetEngine::defineStates() {
    defineOption("", "num-rows", PgOpt::INT_OPT,
                         "The number of rows on the chip.",
                         "-1");
    defineOption("", "num-cols", PgOpt::INT_OPT,
                         "The number of cols on the chip.",
                         "-1");
    defineOption("", "probe-count", PgOpt::INT_OPT,
                         "The number of probes on the chip.",
                         "-1");
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void DmetEngine::checkOptionsImp() {

    defineStates();

    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
    setLibFileOpt("special-snps");
    setLibFileOpt("chrX-probes");
    setLibFileOpt("chrY-probes");
    setLibFileOpt("reference-input");
    setLibFileOpt("probeset-ids");
    setLibFileOpt("probeset-ids-reported");
    setLibFileOpt("region-model");
    setLibFileOpt("probeset-model");
	setLibFileOpt("cn-region-gt-probeset-file");

    if (getOpt("out-dir") == "") {Err::errAbort("Must specify an output directory.");}
    if (getOpt("temp-dir") == "") { 
      setOpt("temp-dir", Fs::join(getOpt("out-dir"),"temp"));
    }

    string cdfFile = getOpt("cdf-file");
    string spfFile = getOpt("spf-file");
    string specialSnps = getOpt("special-snps");
    string chrXProbes = getOpt("chrX-probes");
    string chrYProbes = getOpt("chrY-probes");

	if(getOpt("sample-type") == "plasmid") { 
		setOpt("run-cn-engine","false"); 
	} else { 
		setOpt("run-cn-engine",getOpt("run-cn-engine"));
	}

    string refOut = getOpt("reference-output");
	string batchName = getOpt("batch-name");
    if(refOut == "") {
		if (batchName.empty() != true)
			Err::errAbort("You cannot provide a batch-name when running in single sample mode. batch-name is only valid when output-reference is specified.");
    } else {
		if (batchName.empty() == true)
			Err::errAbort("You must define the batch-name parameter");
        if(Fs::isReadable(refOut))
            if(Fs::rm(refOut, false) != APT_OK)
                Err::errAbort("Unable to remove existing reference-output file '" + refOut + "'");
	}

    ///@todo check chip type in reference file
    ///@todo check reference file version

    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if(celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);

    // Build a consent file if vector of markers was passed in
    vector<string> consented = getOptVector("report");
    if(consented.size() > 0) {
      FsPath probeset_path;
      probeset_path.setPath(getOpt("out-dir"),"probesets-reported","txt");
      //probeset_path.ensureWriteableDirPath();
	  Fs::ensureWriteableDirPath(getOpt("out-dir", false));
      writeProbesetList(probeset_path.asUnixPath(), consented);
      setOpt("probeset-ids-reported", probeset_path.asUnixPath());
    } else {
      setOpt("probeset-ids-reported", getOpt("probeset-ids-reported"));
    }

    if(cdfFile == "" && spfFile == "")
        Err::errAbort("Must specify either a cdf file or spf (simple probe format) file.");
    if (chrXProbes != "" && chrYProbes == "")
        Err::errAbort("Must provide a chrY Probe File when providing a chrX Probe File.");
    if (chrXProbes == "" && chrYProbes != "")
        Err::errAbort("Must provide a chrX Probe File when providing a chrY Probe File.");

    // Check chip types
    vector<string> chipTypesInLayout;

    /* Get the intial info about the chip and check cel files to make sure
       they match. */
    colrow_t numRows = 0, numCols = 0;
    int probeCount=0, probeSetCount=0;

    if(cdfFile != "")
        EngineUtil::getCdfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, cdfFile);
    else if(spfFile != "")
        EngineUtil::getSpfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, spfFile);
    else
        Err::errAbort("Must specify either a cdf file or spf (simple probe format) file.");

    setOpt("num-rows", ToStr(numRows));
    setOpt("num-cols", ToStr(numCols));
    setOpt("probe-count", ToStr(probeCount));

    if(chipTypesInLayout.empty() || chipTypesInLayout[0] == "" || probeCount == 0) 
        Err::errAbort("Problem determining ChipType in file: " + 
              ( cdfFile != "" ? cdfFile : spfFile));

    /* Did the user "force" a set of chip types via options? */
    vector<string> chipTypesSupplied = getOptVector("chip-type");

    /* Figure out what chip type to report */
    if(chipTypesSupplied.size() > 0) {
        setOpt("chip-type", chipTypesSupplied);
    } else if(chipTypesInLayout.size() > 0) {
        setOpt("chip-type", chipTypesInLayout);
    } else {
        Err::errAbort("Unable to figure out a chip type.");
    }

    /* Do Chip Type Check */
    if(!getOptBool("force")) {
        vector<string> chipTypeJustPrimary;
        vector<string> chipTypesToCheck;

        if(chipTypesSupplied.size() > 0) {
            chipTypesToCheck = chipTypesSupplied;
            EngineUtil::checkChipTypeVectors(chipTypesSupplied, chipTypesInLayout);
        } else {
            chipTypesToCheck = chipTypesInLayout;
        }

        chipTypeJustPrimary.push_back(chipTypesToCheck[0]);
        EngineUtil::checkCelChipTypes(chipTypesToCheck, probeCount, celFiles, numRows, numCols);

        // Check special SNPs files
        if (specialSnps != "") {
            EngineUtil::checkTsvFileChipType(specialSnps, chipTypeJustPrimary);
        }
        
        // And other files
        if (chrXProbes != "") {
            EngineUtil::checkTsvFileChipType(chrXProbes, chipTypesToCheck);
        }
        if (chrYProbes != "") {
            EngineUtil::checkTsvFileChipType(chrYProbes, chipTypesToCheck);
        }
    } // end if(!force)
}

bool DmetEngine::generateApsProbesetList(string requestedProbesets, string probesetsInRegion, string output, bool allowProbesetConsent) {

    map<string, string> probesetToRegion;
    map<string, vector<string> > regionProbesets;
    vector<string> probesetsToProcess;

    // Here are the list of probesets in the region
    {
        affx::TsvFile tsv;
        string probeSetName, regionName;
        tsv.bind(0,"probeset_id",&probeSetName,affx::TSV_BIND_REQUIRED);
        tsv.bind(0,"region",&regionName,affx::TSV_BIND_REQUIRED);
        if(tsv.open(probesetsInRegion) != affx::TSV_OK)
            Err::errAbort("Couldn't open file: '" + probesetsInRegion + "' to read.");
        while(tsv.nextLevel(0) == affx::TSV_OK) {
            probesetToRegion[probeSetName] = regionName;
            regionProbesets[regionName].push_back(probeSetName);
        }
    }

    // Here are the list of things requested
    if(requestedProbesets != "") {
        affx::TsvFile tsv;
        string probeSetName;
        tsv.bind(0,"probeset_id",&probeSetName,affx::TSV_BIND_REQUIRED);
        if(tsv.open(requestedProbesets) != affx::TSV_OK)
            Err::errAbort("Couldn't open file: '" + requestedProbesets + "' to read.");
        while(tsv.nextLevel(0) == affx::TSV_OK) {
            if(regionProbesets.find(probeSetName) != regionProbesets.end()) {
                Verbose::out(1,"Adding CN region " + probeSetName + " to the list of regions to process.");
                for(size_t i=0; i<regionProbesets[probeSetName].size(); i++) {
                    probesetsToProcess.push_back(regionProbesets[probeSetName][i]);
                }
            } else if(allowProbesetConsent) {
                Verbose::out(3,"Adding individual probeset " + probeSetName + " to list to process for CN.");
                if(probesetToRegion.find(probeSetName) != probesetToRegion.end())
                    probesetsToProcess.push_back(probeSetName);
            }
        }
    } else {
        // Assume everything is consented
        Verbose::out(1,"Adding all CN probesets to list to process for CN.");
        map<string,string>::iterator iter;
        for(iter=probesetToRegion.begin(); iter != probesetToRegion.end(); iter++)
            probesetsToProcess.push_back(iter->first);
    }

    if(probesetsToProcess.size() < 1) {
        return false;
    } else {
        writeProbesetList(output, probesetsToProcess);
        return true;
    }
}

/// @todo: put this in TsvFile::Util.
void DmetEngine::writeProbesetList(string filename, vector<string> &probesets) {
    affx::TsvFile out;
    out.addHeader("dmet-probeset-list","1");
    out.defineColumn(0, 0, "probeset_id");
    //
    if (out.writeTsv_v1(filename)!=affx::TSV_OK) {
      Err::errAbort("cant open file.");
    }
    //
    for(size_t i=0; i<probesets.size(); i++) {
        out.set(0,0,probesets[i]);
        out.writeLevel(0);
    }
    out.close();
}

/**
 * Compute probeset summaries values for using in computing CN state
 */
bool DmetEngine::computeCnSummaries() {
    ProbesetSummarizeEngine engine;

    vector<string> celFiles = getOptVector("cels");
    engine.setOpt("cels",celFiles);
    vector<string> chipTypes = getOptVector("chip-type");
    if(chipTypes.size() > 0)
        engine.setOpt("chip-type",chipTypes);
    engine.setOpt("cdf-file",getOpt("cdf-file"));
    engine.setOpt("spf-file",getOpt("spf-file"));
    ///@todo If we have a reference, set analysis based on what was used to build the reference
    if(getOpt("ps-analysis") != "")
        engine.setOpt("analysis",getOpt("ps-analysis"));
    else
        engine.setOpt("analysis",
            "quant-norm.sketch=50000,pm-only,plier.optmethod=1.FixFeatureEffect=true,expr.genotype=true");
    engine.setOpt("out-dir",Fs::join(getOpt("out-dir"),"aps"));
    engine.setOpt("verbose",getOpt("verbose"));
    engine.setOpt("force",getOpt("force"));
    engine.setOpt("cc-chp-output","false");

    // We only want to process intersection of probesets in consent file
    // and those in regions
    Fs::ensureWriteableDirPath(getOpt("out-dir", false));

    FsPath probesetFile;
    probesetFile.setDirName(getOpt("out-dir"));
    probesetFile.setFileNameExt("probesets","txt");
    if (!generateApsProbesetList(getOpt("probeset-ids-reported"),
                                 getOpt("probeset-model"),
                                 probesetFile.asUnixPath())) {
      return false;
    }

    engine.setOpt("probeset-ids",probesetFile.asUnixPath());
    engine.setOpt("temp-dir", getOpt("temp-dir"));
    engine.setOpt("use-disk", getOpt("use-disk"));
    engine.setOpt("disk-cache", getOpt("disk-cache"));
    engine.setOpt("set-analysis-name", getOpt("set-analysis-name"));
    engine.setOpt("command-line", getOpt("command-line"));
    engine.setOpt("exec-guid", getOpt("exec-guid"));
    engine.setOpt("program-name", getOpt("program-name"));
    engine.setOpt("program-company", getOpt("program-company"));
    engine.setOpt("program-version", getOpt("program-version"));
    engine.setOpt("program-cvs-id", getOpt("program-cvs-id"));
    engine.setOpt("version-to-report", getOpt("version-to-report"));
    engine.setOpt("summaries","false");

    engine.setOpt("a5-global-file",getOpt("reference-output"));
    //engine.setOpt("a5-global-file-no-replace","true");
    engine.setOpt("a5-group","ProbesetSummarizeEngine");
    engine.setOpt("a5-sketch","true");
    engine.setOpt("a5-summaries","true");
    engine.setOpt("a5-feature-effects","true");
    if(getOpt("reference-output") != "") {
        engine.setOpt("a5-feature-effects-use-global","true");
        engine.setOpt("a5-sketch-use-global","true");
    }
  
    engine.setOpt("a5-global-input-file",getOpt("reference-input"));
    engine.setOpt("a5-input-group","ProbesetSummarizeEngine");
    if(getOpt("reference-input") != "") {
        engine.setOpt("a5-sketch-input-global","true");
        engine.setOpt("a5-feature-effects-input-global","true");
    }

    if(m_ArgvPosAPS > 0) 
        engine.parseArgv(m_argv, m_ArgvPosAPS);
    engine.run();

    return true;
}

/**
 * Compute the CN state
 */
void DmetEngine::computeCnState() {
    DmetCopyNumberEngine engine;

    engine.setOpt("out-dir",Fs::join(getOpt("out-dir"),"adc"));
    engine.setOpt("force",getOpt("force"));
    engine.setOpt("verbose",getOpt("verbose"));
    engine.setOpt("set-analysis-name", getOpt("set-analysis-name"));
    engine.setOpt("command-line", getOpt("command-line"));
    engine.setOpt("exec-guid", getOpt("exec-guid"));
    engine.setOpt("program-name", getOpt("program-name"));
    engine.setOpt("program-company", getOpt("program-company"));
    engine.setOpt("program-version", getOpt("program-version"));
    engine.setOpt("program-cvs-id", getOpt("program-cvs-id"));
    engine.setOpt("version-to-report", getOpt("version-to-report"));

    engine.setOpt("a5-output", "true");
    engine.setOpt("text-output", "true");
    engine.setOpt("region-model", getOpt("region-model"));
    engine.setOpt("cn-region-gt-probeset-file", getOpt("cn-region-gt-probeset-file"));
    engine.setOpt("probeset-model", getOpt("probeset-model"));
    engine.setOpt("a5-summaries",
                  Fs::join(getOpt("out-dir"), "aps",
                               getOpt("set-analysis-name") + ".summary.a5"));

    if(m_ArgvPosCN > 0) 
        engine.parseArgv(m_argv, m_ArgvPosCN);
    engine.run();
}

string setPra(string analysis, int praThresh) {
    ///@todo check that pra-thresh is not already set
    return analysis + ".pra-thresh=" + ToStr(praThresh);
}

/**
 * Compute genotypes
 */
void DmetEngine::computeGenotypes() {
    ProbesetGenotypeEngine engine;

    vector<string> celFiles = getOptVector("cels");
    engine.setOpt("cels",celFiles);
    vector<string> chipTypes = getOptVector("chip-type");
    if(chipTypes.size() > 0)
        engine.setOpt("chip-type",chipTypes);
    engine.setOpt("cdf-file",getOpt("cdf-file"));
    engine.setOpt("spf-file",getOpt("spf-file"));
    string analysis;
    if(getOpt("gt-analysis") != "") {
        analysis = getOpt("gt-analysis");
    } else {
        ///@todo If we have a reference, set analysis and qmethod-spec based on 
        ///      what was used to build the reference
        if(getOpt("reference-output") != "") {
            analysis = "quant-norm.sketch=50000,pm-only,brlmm-p-multi.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.KX=0.3.KH=0.3.KXX=0.1.KAH=-0.1.KHB=-0.1.KYAH=-0.05.KYHB=-0.05.KYAB=-0.1.transform=MVA.AAM=2.8.BBM=-2.8.AAV=0.10.BBV=0.10.ABV=0.10.V=1.AAY=10.7.ABY=11.3.BBY=10.7.copyqc=0.00000.wobble=0.05.MS=0.1.copytype=-1.clustertype=2.CSepPen=0.5.ocean=0.0000001";
        } else {
            analysis = "quant-norm.sketch=50000,pm-only,brlmm-p-multi.CM=2.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.KX=0.3.KH=0.3.KXX=0.1.KAH=-0.1.KHB=-0.1.KYAH=-0.05.KYHB=-0.05.KYAB=-0.1.transform=MVA.AAM=2.8.BBM=-2.8.AAV=0.10.BBV=0.10.ABV=0.10.V=1.AAY=10.7.ABY=11.3.BBY=10.7.copyqc=0.00000.wobble=0.05.MS=0.1.copytype=-1.clustertype=2.CSepPen=0.5.ocean=0.0000001";
        }
        analysis += ".cc-alleles=6";
        analysis += ".cc-type=UCHAR";
        analysis += ".cc-version=1.0";
    }
    if(getOpt("gt-qmethod-spec") != "")
        engine.setOpt("qmethod-spec",getOpt("gt-qmethod-spec"));
    else
        engine.setOpt("qmethod-spec","plier.optmethod=1.FixFeatureEffect=true");

    analysis = setPra(analysis,getOptInt("pra-thresh"));
    engine.setOpt("analysis",analysis);
    engine.setOpt("out-dir",Fs::join(getOpt("out-dir"),"apg"));
    engine.setOpt("verbose",getOpt("verbose"));
    engine.setOpt("force",getOpt("force"));
    engine.setOpt("cc-chp-output","false");
    engine.setOpt("probeset-ids-reported",getOpt("probeset-ids-reported"));
    engine.setOpt("probeset-ids",getOpt("probeset-ids"));
    engine.setOpt("temp-dir", getOpt("temp-dir"));
    engine.setOpt("use-disk", getOpt("use-disk"));
    engine.setOpt("disk-cache", getOpt("disk-cache"));
    engine.setOpt("set-analysis-name", getOpt("set-analysis-name"));
    engine.setOpt("call-coder-max-alleles", getOpt("call-coder-max-alleles"));
    engine.setOpt("call-coder-type", getOpt("call-coder-type"));
    engine.setOpt("call-coder-version", getOpt("call-coder-version"));
    engine.setOpt("command-line", getOpt("command-line"));
    engine.setOpt("exec-guid", getOpt("exec-guid"));
    engine.setOpt("program-name", getOpt("program-name"));
    engine.setOpt("program-company", getOpt("program-company"));
    engine.setOpt("program-version", getOpt("program-version"));
    engine.setOpt("program-cvs-id", getOpt("program-cvs-id"));
    engine.setOpt("version-to-report", getOpt("version-to-report"));
    engine.setOpt("prior-size", "1");
    engine.setOpt("special-snps", getOpt("special-snps"));
    engine.setOpt("em-gender", "false");
    if(getOpt("chrX-probes") != "") {
        engine.setOpt("chrX-probes", getOpt("chrX-probes"));
        engine.setOpt("chrY-probes", getOpt("chrY-probes"));
        engine.setOpt("set-gender-method", "cn-probe-chrXY-ratio");
        engine.setOpt("no-gender-force", "false");
    } else {
        engine.setOpt("set-gender-method", "none");
        engine.setOpt("no-gender-force", "true");
    }
    engine.setOpt("female-thresh", getOpt("female-thresh"));
    engine.setOpt("male-thresh", getOpt("male-thresh"));
    engine.setOpt("table-output", "false");
    engine.setOpt("output-context", "true");
    engine.setOpt("output-forced-calls", "true");

    if (getOptBool("run-cn-engine") && getOpt("cn-region-gt-probeset-file") != "") {
      /// @todo perhaps we should use state to track the file name
      ///       as it is we need to keep this in sync with the cn engine call
      engine.setOpt("genotype-markers-cn-file", 
                    Fs::join(getOpt("out-dir"),
                                 "adc",
                                 getOpt("set-analysis-name") + ".gt.markers.cn.call.a5"));
    }

    engine.setOpt("a5-global-file",getOpt("reference-output"));
    engine.setOpt("a5-global-file-no-replace","true");
    engine.setOpt("a5-group","ProbesetGenotypeEngine");
    engine.setOpt("a5-calls","true");
    engine.setOpt("a5-summaries","true");
    engine.setOpt("a5-feature-effects","true");
    engine.setOpt("a5-sketch","true");
    engine.setOpt("a5-write-models","true");
    if(getOpt("reference-output") != "") {
        engine.setOpt("a5-sketch-use-global","true");
        engine.setOpt("a5-feature-effects-use-global","true");
        engine.setOpt("a5-write-models-use-global","true");
    }
  
    engine.setOpt("a5-global-input-file",getOpt("reference-input"));
    engine.setOpt("a5-input-group","ProbesetGenotypeEngine");
    if(getOpt("reference-input") != "") {
        engine.setOpt("a5-sketch-input-global","true");
        engine.setOpt("a5-feature-effects-input-global","true");
        engine.setOpt("a5-models-input-global","true");
    }

    if(m_ArgvPosAPG > 0) 
        engine.parseArgv(m_argv, m_ArgvPosAPG);
    engine.run();
}

string DmetEngine::makeVal(string name, int optionIndex) {
    vector<string> vals;
    vals = getOptVector(name, optionIndex);

    /*
    if(name.find("file") != string::npos)
        for(int i=0; i<vals.size(); i++)
            vals[i] = Fs::basename(vals[i]);
    */
    string val;
    if(vals.size() > 1) {
        val = "'" + vals[0] + "'";
        for(int i=1; i<vals.size(); i++) {
            val += ",'" + vals[i] + "'";
        }
    } else if(vals.size() == 1) {
        val = vals[0];
    } else {
        val = getOpt(name,optionIndex);
    }
    return val;
}

void DmetEngine::fillInAnalysisInfo(AnalysisInfo &info) {
    string prefix = "apt-";

    vector<string> celFiles = getOptVector("cels");
    info.m_AlgVersion = "3.0";
    info.m_AlgName = getOpt("set-analysis-name");
    info.m_ProgramName = getOpt("program-name");
    info.m_ProgramVersion = getOpt("version-to-report");
    info.m_ProgramCompany = getOpt("program-company");
    info.m_ChipType = getOpt("chip-type");
    info.m_ProgID = "";
    info.m_ExecGuid = getOpt("exec-guid");
    info.m_AnalysisGuid = getOpt("analysis-guid");

    // Execution info
    info.addParam("apt-engine", "DmetEngine");
    info.addParam(prefix + "program-name", getOpt("program-name"));
    info.addParam(prefix + "command-line", getOpt("command-line"));
    info.addParam(prefix + "exec-guid", getOpt("exec-guid"));
    info.addParam(prefix + "analysis-guid", getOpt("analysis-guid"));
    info.addParam(prefix + "time-str", getOpt("time-start"));
    info.addParam(prefix + "version", getOpt("version-to-report"));
    info.addParam(prefix + "cvs-id", getOpt("program-cvs-id"));

    // Engine Options
    string opt = "opt-";
    vector<string> initialOptionNames;
    getOptionNames(initialOptionNames,1);
    for(int i=0; i<initialOptionNames.size(); i++) {
        string name = initialOptionNames[i];
        if((name != "cels") || getOptBool("batch-info")){
            info.addParam(prefix + opt + name, makeVal(name,1));
        }
    }

    // Cel files in batch
    if(getOptBool("batch-info")){
        for(uint32_t i = 0; i < celFiles.size(); i++) {
            std::string paramName = prefix + opt + "cel-" + ToStr(i+1);
            info.addParam(paramName, Fs::basename(celFiles[i]));
        
            affymetrix_fusion_io::FusionCELData cel;
            try {
                cel.SetFileName(celFiles[i].c_str());
                if(!cel.ReadHeader()) {
                    Err::errAbort("Unable to read CEL file " + celFiles[i]);
                }
                affymetrix_calvin_io::GenericData *gdata = cel.GetGenericData();
                if (gdata != NULL)
                {
                    std::string paramName = prefix + opt + "cel-guid-" + ToStr(i+1);
                    info.addParam(paramName, gdata->Header().GetGenericDataHdr()->GetFileId());
                }
                cel.Close();
            }
            catch (...)
            {
                Err::errAbort("Unable to read CEL file " + celFiles[i]);
            }
        }
    }
    
    // Engine State
    string option = "state-";
    std::vector<std::string> optionNames;
    getOptionNames(optionNames);
    for(int i=0; i<optionNames.size(); i++) {
        string name = optionNames[i];
        if((name != "cels") || getOptBool("batch-info")){
            info.addParam(prefix + option + name, makeVal(name));
        }
    }

    // Sanity check
    Err::check(info.m_ParamValues.size() == info.m_ParamNames.size(),
             "AnalysisInfo - Names and values out of sync.");
}

/**
 * Build final CHP files
 */
void DmetEngine::generateChpFiles() {

  DmetCHPWriter engine;

  AnalysisInfo info;
  fillInAnalysisInfo(info);
  engine.setAnalysisInfo(info);

  vector<string> celFiles = getOptVector("cels");
  engine.setOpt("cels",celFiles);
  engine.setOpt("out-dir",Fs::join(getOpt("out-dir"),"chp"));
  engine.setOpt("verbose",getOpt("verbose"));
  engine.setOpt("set-analysis-name", getOpt("set-analysis-name"));
  engine.setOpt("exec-guid", getOpt("exec-guid"));
  engine.setOpt("program-name", getOpt("program-name"));
  engine.setOpt("program-company", getOpt("program-company"));
  engine.setOpt("program-version", getOpt("program-version"));
  engine.setOpt("batch-name", getOpt("batch-name"));

  if (getOptBool("run-cn-engine")) {
    engine.setOpt("a5-copynumber",
                  Fs::join(getOpt("out-dir"), "adc",
                               getOpt("set-analysis-name") + ".copynumber.a5"));
  }

  engine.setOpt("a5-summaries",
                Fs::join(getOpt("out-dir"), "apg",
                             getOpt("set-analysis-name") + ".summary.a5"));
  engine.setOpt("a5-calls", 
                Fs::join(getOpt("out-dir"), "apg",
                             getOpt("set-analysis-name") + ".calls.a5"));
  engine.setOpt("a5-forced-calls", 
                Fs::join(getOpt("out-dir"), "apg", getOpt("set-analysis-name") + ".forced-calls.a5"));
  engine.setOpt("a5-confidences", 
                Fs::join(getOpt("out-dir"), "apg",
                             getOpt("set-analysis-name") + ".confidences.a5"));
  engine.setOpt("a5-context", 
                Fs::join(getOpt("out-dir"), "apg",
                             getOpt("set-analysis-name") + ".context.a5"));
  engine.setOpt("report-file", 
                Fs::join(getOpt("out-dir"), "apg",
                             getOpt("set-analysis-name") + ".report.txt"));

  engine.setOpt("spf-file",getOpt("spf-file"));
  engine.setOpt("cdf-file",getOpt("cdf-file"));
  engine.setOpt("region-model", getOpt("region-model"));
  engine.setOpt("null-context", getOpt("null-context"));
  engine.setOpt("geno-call-thresh", getOpt("geno-call-thresh"));
  engine.setOpt("call-coder-max-alleles", getOpt("call-coder-max-alleles"));
  engine.setOpt("call-coder-type", getOpt("call-coder-type"));
  engine.setOpt("call-coder-version", getOpt("call-coder-version"));

  engine.run();
}

/**
 * compute CHP files
 */
void DmetEngine::runImp() {

    setOpt("analysis-guid", affxutil::Guid::GenerateNewGuid());

    string errMsg;

    if(!Fs::isWriteableDir(getOpt("out-dir")))
        if(Fs::mkdirPath(getOpt("out-dir"), false) != APT_OK)
            Err::errAbort("Can't make or write to directory: " + getOpt("out-dir"));
    if(getOptBool("run-cn-engine")) {
        Verbose::out(1,"");
        Verbose::out(1,"Step 1: Computing probeset summaries for copy number state calling");
        if(computeCnSummaries()) {
            Verbose::out(1,"");
            Verbose::out(1,"Step 2: Computing copy number states");
            computeCnState();
            setOpt("run-cn-engine","true");
        } else {
            Verbose::out(1,"");
            Verbose::out(1,"No CN regions to compute. Skipping step 2.");
            setOpt("run-cn-engine","false");
        }
    } else {
        Verbose::out(1,"");
        Verbose::out(1,"Not computing CN state. Skipping steps 1 and 2.");
        setOpt("run-cn-engine","false");
    }
    Verbose::out(1,"");
    Verbose::out(1,"Step 3: Computing genotypes");
    computeGenotypes();
    if(getOptBool("cc-chp-output")) {
        Verbose::out(1,"");
        Verbose::out(1,"Step 4: Generating CHP files");
        generateChpFiles();
    }
    Verbose::out(1,"");
    Verbose::out(1, "Done.");
}

