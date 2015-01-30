////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
#include "canary/CanaryEngine.h"
//
#include "canary/CanaryOptions.h"
#include "canary/CanaryWithPrior.h"
#include "canary/canary_utils.h"
//
#include "calvin_files/data/src/ProbeSetMultiDataData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/SparseMart.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include "newmat.h"
//
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <valarray>
#include <vector>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CanaryEngine::Reg CanaryEngine::reg;

CanaryEngine* CanaryEngine::FromBase(BaseEngine *engine) {
    if (engine != NULL && engine->getEngineName() == CanaryEngine::EngineName())
        return (CanaryEngine *)engine;
    return NULL;
}

CanaryEngine::CanaryEngine() {
    defineOptions();
}

CanaryEngine::~CanaryEngine() {
}

/** Fill in options given command line arguments. */
void CanaryEngine::fillInOptions(CanaryOptions &o) {
    o.version = getOpt("program-version");
    o.cvsId = getOpt("program-cvs-id");
    o.progName = getOpt("program-name");
    o.progName = Fs::basename(o.progName);
    o.execGuid = getOpt("exec-guid");
    o.timeStr = getOpt("time-start");

    o.genomeVersion = getOpt("genome-version");
    o.mapName = getOpt("map-name");
    o.mapVersion = getOpt("map-version");

    o.force = getOptBool("force");
    o.precision = getOptInt("precision");
    o.verbosity = getOptInt("verbose");
    o.cdfFile = getOpt("cdf-file");
    o.spfFile = getOpt("spf-file");
    o.outDir = getOpt("out-dir");
    o.cnvRegionFile = getOpt("cnv-region-file");
    o.cnvMapFile = getOpt("cnv-map-file");
    o.canaryPriorFile = getOpt("cnv-prior-file");
    o.cnvNormFile = getOpt("cnv-normalization-file");
    o.setAnalysisName = getOpt("analysis-name");
    o.ccMDChpOutput = getOptBool("cc-chp-output");
    o.tableOutput = getOptBool("table-output");

    o.chipType = getOpt("chip-type");

    Util::chompLastIfSep(o.outDir); // for windows

    if (getOpt("apt-summarize-analysis") != "")
        o.aptSummarizeAnalysis = getOpt("apt-summarize-analysis");

    if (getOpt("apt-canary-analysis") != "")
        o.aptCanaryAnalysis = getOpt("apt-canary-analysis");

    o.chpFiles = getOptVector("result-files");
    o.celFiles = getOptVector("cels");

    // Get CEL file GUIDs
    o.celGuids.resize(o.celFiles.size());
    for (int chip=0; chip<o.celFiles.size(); chip++) {
        FusionCELData cel;
        try {
            cel.SetFileName(o.celFiles[chip].c_str());
            if (!cel.ReadHeader()) {
                Err::errAbort("Unable to read CEL file " + o.celFiles[chip]);
            }
            GenericData *gdata = cel.GetGenericData();
            if (gdata != NULL) {
                o.celGuids[chip] = gdata->Header().GetGenericDataHdr()->GetFileId();
            }
            cel.Close();
        }
        catch (...) {
            Err::errAbort("Unable to read CEL file " + o.celFiles[chip]);
        }
    }

    return;
}


void CanaryEngine::extraHelp() {
    cout << '\n';
    cout << "Canary Parameters:\n";
    cout << '\n';
    cout << "    The canary algorithm does copy number calling on defined\n";
    cout << "    regions using priors. The following parameters are \n";
    cout << "    accessible using the --apt-canary-analysis option.\n";
    cout << "    Use key1=val1,key2=val2,... string format.\n";
    cout << '\n';
    cout << "        af-weight \n";
    cout << "        TOL \n";
    cout << "        hwe_tol \n";
    cout << "        hwe_tol2 \n";
    cout << "        fraction-giveaway-0 \n";
    cout << "        fraction-giveaway-1 \n";
    cout << "        fraction-giveaway-2 \n";
    cout << "        fraction-giveaway-3 \n";
    cout << "        fraction-giveaway-4 \n";
    cout << "        min-fill-prop \n";
    cout << "        conf-interval-half-width \n";
    cout << "        inflation \n";
    cout << "        min-cluster-variance \n";
    cout << "        pseudopoint-factor \n";
    cout << "        regularize_variance_factor \n";
    cout << '\n';
}

void CanaryEngine::defineOptions() {
    // Input options
    defineOptionSection("Input Options");

    defineOption("", "cel-files", PgOpt::STRING_OPT, 
                 "Text file specifying cel files to process, one per line with the first line being 'cel_files'.", 
                 "");
    defineOption("", "cdf-file", PgOpt::STRING_OPT, 
                 "File defining probe sets. Use either --cdf-file or --spf-file ", 
                 "");
    defineOption("", "spf-file", PgOpt::STRING_OPT,
                 "File defining probe sets in spf (simple probe format) which is like a text cdf file.",
                 "");
    defineOption("", "cnv-region-file", PgOpt::STRING_OPT,
                 "File defining CNV regions and what probesets to use for each CNV region.", 
                 "");
    defineOption("", "cnv-prior-file", PgOpt::STRING_OPT,
                 "File defining the canary priors for a given CNV regions file.", 
                 "");
    defineOption("", "cnv-map-file", PgOpt::STRING_OPT,
                 "File (bed format) used for visualizing CNV regions in other applications. This arg causes the map file name to be included in the CHP meta info.", 
                 "");
    defineOption("", "cnv-normalization-file", PgOpt::STRING_OPT,
                 "File containing probesets to use (restricted to) for doing probe level normalization.", 
                 "");
    defOptMult("", "chip-type", PgOpt::STRING_OPT,
               "Chip types to check library and CEL files against. "
               "Can be specified multiple times. "
               "The first one is propigated as the chip type in the output files. "
               "Warning, use of this option will override the usual check between chip types "
               "found in the library files and cel files. You should use this option instead "
               "of --force when possible. ",
               "");


    // Output options
    defineOptionSection("Output Options");

    defineOption("", "table-output", PgOpt::BOOL_OPT,
                 "Output matching matrices of tab delimited genotype calls and confidences.", "true");
    defineOption("", "cc-chp-output", PgOpt::BOOL_OPT,
                 "Output resulting calls in binary CHP format. "
                 "This makes one AGCC Multi Data CHP file per cel file analyzed.", "false");

    // Analysis options
    defineOptionSection("Analysis Options");

    defineOption("", "apt-summarize-analysis", PgOpt::STRING_OPT,
                 "String representing analysis parameters for the "
                 "apt-probeset-summarize step a.k.a. pre-canary",
                 "");
    defineOption("", "apt-canary-analysis", PgOpt::STRING_OPT,
                 "String representing analysis parameters for canary.",
                 "");

    // Execution control options
    defineOptionSection("Execution Control Options");

    defineOption("", "precision", PgOpt::INT_OPT,
                 "Precision after decimal place", "4");

    defineOption("", "analysis-name", PgOpt::STRING_OPT,
                 "Set the name of the analysis.", "");
    defineOption("", "use-disk", PgOpt::BOOL_OPT,
                 "Store CEL intensities to be analyzed on disk.", "true");
    defineOption("", "disk-cache", PgOpt::INT_OPT,
                 "Size of intensity memory cache in millions of intensities (when --use-disk=true).",
                 "50");

    defineOptionSection("Engine Options (Not used on command line)");
    defOptMult("", "cels", PgOpt::STRING_OPT,
               "Cel files to process.",
               "");
    defOptMult("", "result-files", PgOpt::STRING_OPT,
               "CHP file names to output. Must be paired with cels.",
               "");
}

void CanaryEngine::defineStates() {
    defineOption("", "num-rows", PgOpt::INT_OPT,
                 "The number of rows on the chip.",
                 "-1");
    defineOption("", "num-cols", PgOpt::INT_OPT,
                 "The number of cols on the chip.",
                 "-1");
    defineOption("", "probe-count", PgOpt::INT_OPT,
                 "The number of probes on the chip.",
                 "-1");
    defineOption("", "channel-count", PgOpt::INT_OPT,
                 "The number of data channels on the chip.",
                 "-1");
    defineOption("", "map-name", PgOpt::STRING_OPT,
                 "The map label",
                 "");
    defineOption("", "map-version", PgOpt::STRING_OPT,
                 "Version of the map file.",
                 "");
    defineOption("","genome-version", PgOpt::STRING_OPT,
                 "The version of the genome.",
                 "");
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 * @param opts - Options to be checked.
 */

void CanaryEngine::checkOptionsImp() {

    defineStates();

    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
    setLibFileOpt("cnv-region-file");
    setLibFileOpt("cnv-prior-file");
    setLibFileOpt("cnv-normalization-file");

    if (getOpt("out-dir") == "") {
        Err::errAbort("Must specify an output directory.");
    }
    if (getOpt("temp-dir") == "") { 
      setOpt("temp-dir",Fs::join(getOpt("out-dir"),"temp"));
    }

    if (getOpt("cdf-file") == "" && getOpt("spf-file") == "")
        Err::errAbort("Must specify either a cdf file or spf file.");

    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if (celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);

    vector<string> resultFiles = getOptVector("result-files");
    if (resultFiles.size() > 0 && celFiles.size() != resultFiles.size()) {
        Err::errAbort("result-files option used but is not the same size as cel file listing");
    }

    if (getOpt("cnv-region-file") == "")
        Err::errAbort("Must specify a CNV region file");
    if (getOpt("cnv-normalization-file") == "")
        Err::errAbort("Must specify normalization file");

    // Check chip types
    vector<string> chipTypesInLayout;

    /* Get the intial info about the chip and check cel files to make sure they match. */
    colrow_t numRows = 0, numCols = 0;
    int probeCount = 0, probeSetCount = 0, channelCount = 1;

    string cdfFile = getOpt("cdf-file");
    string spfFile = getOpt("spf-file");
    if (cdfFile!="")
        EngineUtil::getCdfChipType(chipTypesInLayout, numRows, numCols, 
                                   probeCount, probeSetCount, cdfFile);
    else if (spfFile!="")
        EngineUtil::getSpfChipType(chipTypesInLayout, numRows, numCols, 
                                   probeCount, probeSetCount, spfFile);
    else
        Err::errAbort("Must specify a cdf file, spf file, or PGF and CLF files.");

    setOpt("num-rows", ToStr(numRows));
    setOpt("num-cols", ToStr(numCols));
    setOpt("probe-count", ToStr(probeCount));
    setOpt("channel-count", ToStr(channelCount));

    if (chipTypesInLayout.empty() || chipTypesInLayout[0] == "" || 
        probeCount == 0) 
        Err::errAbort("Problem determining ChipType in file: " + 
                      ( cdfFile != "" ? cdfFile : spfFile));

    /* Did the user "force" a set of chip types via options? */
    vector<string> chipTypesSupplied = getOptVector("chip-type");

    /* Figure out what chip type to report */
    if (chipTypesSupplied.size() > 0) {
        setOpt("chip-type", chipTypesSupplied[0]);
    } else if (chipTypesInLayout.size() > 0) {
        setOpt("chip-type", chipTypesInLayout[0]);
    } else {
        Err::errAbort("Unable to figure out a chip type.");
    }

    vector<string> chipTypesToCheck;
    if (chipTypesSupplied.size() > 0) {
        chipTypesToCheck = chipTypesSupplied;
    } else {
        chipTypesToCheck = chipTypesInLayout;
    }

    /* Do Chip Type Check */
    if (!getOptBool("force")) {
        if (chipTypesSupplied.size() > 0) {
            EngineUtil::checkChipTypeVectors(chipTypesSupplied, chipTypesInLayout);
        }
        EngineUtil::checkCelChipTypes(chipTypesToCheck, probeCount, celFiles, numRows, numCols);
    }

    // check the map files and get the map name and version
    vector<string> files;
    files.push_back(getOpt("cnv-region-file"));
    files.push_back(getOpt("cnv-map-file"));
    files.push_back(getOpt("cnv-normalization-file"));
    files.push_back(getOpt("cnv-prior-file"));
    string mapName, mapVersion;
    check_input_file_headers(chipTypesToCheck, files, mapName, mapVersion, 
                             getOptBool("force"));
    setOpt("map-name", mapName);
    setOpt("map-version",mapVersion);

    // default analysis name is pulled from the regions file
    if (getOpt("analysis-name") == "")
        setOpt("analysis-name", mapName + "-" + mapVersion);
    else
        setOpt("analysis-name",getOpt("analysis-name"));

    if (getOpt("exec-guid") == "")
        setOpt("exec-guid", affxutil::Guid::GenerateNewGuid());
    else
        setOpt("exec-guid", getOpt("exec-guid"));

    // Always set this -- as we own it
    setOpt("analysis-guid", affxutil::Guid::GenerateNewGuid()) ;
    setOpt("genome-version", get_genome_version(getOpt("cnv-region-file")));

    if (!Fs::isWriteableDir(getOpt("out-dir"))) 	 
        if (Fs::mkdirPath(getOpt("out-dir"), false) != APT_OK) 	 
            Err::errAbort("Can't make or write to directory: " + 
                          getOpt("out-dir")); 	 
    makeTempDir(getOpt("temp-dir"));
}

/**
 * @brief checks if there is enough disk space for diskmarts in the
 * temp-dir and for chp files in out-dir
 */
void CanaryEngine::checkDiskSpaceImp() {
    int64_t out_disk_available = -1, scratch_disk_available = -1;
  
    int cel_count = getOptVector("cels").size();
    int row_count = getOptInt("num-rows");
    int col_count = getOptInt("num-cols");
    int channel_count = getOptInt("channel-count");
    // @todo get analysis_stream_count, if necessary
    // int analysis_stream_count = getOptInt("analysis-count");
    int analysis_stream_count = 1;

    std::vector<std::string> result_files = getOptVector("result-files");
    std::string out_dir;
    if (!result_files.empty()) {
      out_dir = Fs::dirname(result_files[0]);
      ///@todo iterate through all of the result dirs to check available
      // disk space
    }
    
    ///@todo get name of parameter that specifies chp-out-dir, if any
//   if (out_dir.empty()) {
//     out_dir = getOpt("cc-chp-out-dir");
//   }
   
    if (out_dir.empty()) {
        out_dir = getOpt("out-dir");
    }
        
    out_disk_available = Fs::getFreeDiskSpace(out_dir, false);
        
    uint64_t out_disk_needed = 0;
    std::string region_file = getOpt("cnv-region-file");
    if (!region_file.empty()) {
        TsvFile tsv;
        if (tsv.open(region_file) != TSV_OK) {
            Err::errAbort("Couldn't open file: '" + region_file + "' to read.");
        }
        int num_regions = tsv.countTotalDataLines();
        out_disk_needed = static_cast<uint64_t>(num_regions) * cel_count * 38.39;
    }
        
    // (1 + analysis_stream_count): one initial raw diskmart + one
    // diskmarts for each analysis stream
    uint64_t scratch_disk_needed = 0;
    if (getOptBool("use-disk")) {
        scratch_disk_needed = static_cast<uint64_t>(row_count) * col_count * cel_count * channel_count * 4.04 * (1 + analysis_stream_count);
    }

    std::string temp_dir = getOpt("temp-dir");
    AptErr_t aptErr = APT_OK;
    bool same_disk = Fs::isSameVolume(temp_dir, out_dir, aptErr, false);
    if (temp_dir.empty() ||
        same_disk || 
        ((aptErr != APT_OK) && (temp_dir == out_dir))
        ) {
        out_disk_needed += scratch_disk_needed;
    }
    else {
        scratch_disk_available = Fs::getFreeDiskSpace(temp_dir, false);
        if (scratch_disk_needed > 0 && 
            scratch_disk_available >= 0 &&
            scratch_disk_needed >= scratch_disk_available) {
            // format in kb
            scratch_disk_needed = (scratch_disk_needed + 512) / 1024;
            scratch_disk_available = (scratch_disk_available + 512) / 1024;
            Err::errAbort("In " + temp_dir + ", need " + 
                          ToStr(scratch_disk_needed) + 
                          "kb of scratch diskspace.  Only " + 
                          ToStr(scratch_disk_available) + "kb available.");
        }
    }
        
    if (out_disk_needed > 0 &&
        out_disk_available >= 0 &&
        out_disk_needed >= out_disk_available) {
        // format in kb
        out_disk_needed = (out_disk_needed + 512) / 1024;
        out_disk_available = (out_disk_available + 512) / 1024;
        Err::errAbort("In " + out_dir + ", need " + ToStr(out_disk_needed) + 
                      "kb of diskspace.  Only " + ToStr(out_disk_available) + 
                      "kb available.");
    }
}

void CanaryEngine::runImp() {
    CanaryOptions myOpts;
    fillInOptions(myOpts);
    myOpts.setTuneableStr();
    CanaryCompute(myOpts);
    removeTempDir(getOpt("temp-dir"));
}

/*
 * Read the design of the chip from the cdf or spf file.
 *
 * @param layout       - Object to be filled in.
 * @param opts         - Supplied engine options
 * @param libpsets     - Set of probeset names to be filled in
 */

void CanaryEngine::loadCdfLayout_canary(ChipLayout & layout, 
                                        CanaryOptions & opts,
                                        set<string> & libpsets) {
    set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    set<const char *, Util::ltstr> probeSetsToLoad;
    vector<const char *> probesetNames;
    vector<bool> probeSubset;

    if (opts.cdfFile != "") {
        Verbose::out(1,"Loading chip layout information from CDF file.");
        probeidmap_t dummy_killList;
        if (!layout.openCdf(opts.cdfFile, probeSetsToLoad,
                            &probesetNames, probeSubset, "", dummy_killList, 
                            false,
                            psTypesToLoad))
            Err::errAbort("Unable to open CDF file.");
    }

    else if (opts.spfFile != "") {
        Verbose::out(1,"Loading chip layout information from SPF file.");
        layout.openSpf(opts.spfFile, probeSetsToLoad, &probesetNames,
                       probeSubset, "", false, psTypesToLoad);
    }

    else Err::errAbort("Must have either a cdf file or spf file.");

    libpsets.clear();
    for (unsigned int k=0; k<probesetNames.size(); k++) {
        libpsets.insert(probesetNames[k]);
        delete[] probesetNames[k];
        probesetNames[k] = NULL;
    }
}


/*
 * Return the list of cnv probes to run through pre-canary
 */

set<string> CanaryEngine::get_cnv_psets(map<string,vector<string> > &cnv_region_map) {
    set<string>probe_names;

    // walk the map and load up probe_names
    map<string,vector<string> >::iterator iter;
    for (iter=cnv_region_map.begin(); iter != cnv_region_map.end(); iter++) {
        for (int i=0; i<iter->second.size(); i++)
            probe_names.insert(iter->second[i]);
    }

    return probe_names;
}


/*
 * Return the list of normalization to use
 */

set<string> CanaryEngine::get_normalization_psets(string ifname) {
    set<string>probe_names;
  
    TsvFile cn_normfile;
    string the_probe;
    cn_normfile.bind(0,"probeset_id",&the_probe,TSV_BIND_REQUIRED);
    if (cn_normfile.open(ifname) != TSV_OK)
        Err::errAbort("Couldn't open file: " + ifname + " to read.");

    while (cn_normfile.nextLevel(0)==TSV_OK) {
	cn_normfile.get(cn_normfile.lineLevel(),0,the_probe);
	probe_names.insert(the_probe);
    }

    return probe_names;
}

/*
 * Returns a boolean vector of probes to be analysed according to
 * the set of probe set names
 */

vector<string> CanaryEngine::get_available_psets(set<string> & cnv_probes,
                                                 set<string> & normalization_probes, 
                                                 set<string> & all_probes,
                                                 CanaryOptions & opts) {
    set<string>psets;
    set<string>not_found;

    for (set<string>::iterator iter=cnv_probes.begin();
         iter!=cnv_probes.end(); iter++) {
        if (all_probes.find(*iter) != all_probes.end()) psets.insert(*iter);
        else not_found.insert(*iter);
    }

    for (set<string>::iterator iter=normalization_probes.begin();
         iter!=normalization_probes.end(); iter++) {
        if (all_probes.find(*iter) != all_probes.end()) psets.insert(*iter);
        else not_found.insert(*iter);
    }

    vector<string>pset_vec;
    for (set<string>::iterator iter=psets.begin(); iter!=psets.end(); iter++) {
        pset_vec.push_back(*iter);
    }

    /*
      Used to be able to put missing probes in a file.  To conform with
      GTC just write them to verbose output.
    */
    for (set<string>::iterator iter=not_found.begin();
         iter!=not_found.end(); iter++) {
        string vstr = *iter + ": probe not found in chip definition";
        Verbose::out(2,vstr.c_str());
    }

    return pset_vec;
}


void CanaryEngine::averageAlleles(CanaryOptions &opts) {
    // Input and output file names
    ///@todo Next time maybe not hard code file name extensions in 
    // multiple places.
    string ifname = Fs::join(opts.outDir, opts.setAnalysisName + ".summary.txt");
    string ofname = Fs::join(opts.outDir, opts.setAnalysisName + ".allele_average_summary.txt");

    // Open the input file and get the number of columns
    TsvFile infile;
    if (infile.open(ifname) != TSV_OK)
        Err::errAbort("Unable to open allele signal input file '" + ifname + "'");
    int ncol = infile.getColumnCount(0);

    // Declare the output TsvFile
    TsvFile outfile;
//	outfile.setPrecision(opts.precision);  // kicks out an error on compile

    // Define the column names from the old file
    for (int i=0; i<ncol; i++) {
        string cname = infile.getColumnName(0,i);
        outfile.defineColumn(0,i,cname);
    }

    outfile.writeTsv_v1(ofname);

    while (infile.nextLevel(0)==TSV_OK) {
        string probe_name;
        double value;
        infile.get(0,0,probe_name);

        // Just write out the values again
        if (probe_name[0] == 'C') {
            outfile.set(0,0,probe_name);
            for (int i=1; i<ncol; i++) {
                infile.get(0,i,value);
                outfile.set(0,i,value);
            }
        }

        else if (probe_name[0] == 'S') {
            probe_name.erase(probe_name.size()-2);
            outfile.set(0,0,probe_name);
            valarray<double> valuesA(ncol-1);
            valarray<double> valuesB(ncol-1);
            for (int i=0; i<ncol-1; i++) infile.get(0,i+1,valuesA[i]);
            infile.nextLevel(0);
            for (int i=0; i<ncol-1; i++) infile.get(0,i+1,valuesB[i]);

            for (int i=0; i<ncol-1; i++) {
                double ABavg = 0.5*(valuesA[i] + valuesB[i]);
                outfile.set(0,i+1,ABavg);
            }
        }

        else
        {
            string vstr = probe_name +
                " not recognized as copy number or SNP when averaging alleles";
            Verbose::out(2,vstr);
        }
        outfile.writeLevel(0);
    }

    infile.close();
    outfile.close();

    // get rid of the apt-probeset-summarize input
    if (!opts.tableOutput)
        Fs::rm(ifname, false);
}

void CanaryEngine::get_meta_tags(string file, vector<string> &keys, 
                                 vector<string> &vals) {
    TsvFile tsv;
    string key, val;
    if (tsv.open(file) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + file + "' to read.");
    tsv.headersBegin();
    while (tsv.headersNext(key,val)==affx::TSV_OK) {
        keys.push_back(key);
        vals.push_back(val);
    }
    tsv.close();
}

string CanaryEngine::get_genome_version(string cnvRegionFile) {
    TsvFile tsv;
    string genomeVersion;
    if (tsv.open(cnvRegionFile) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + cnvRegionFile + "' to read.");
    if (tsv.headersFindNext("genome-version", genomeVersion) != TSV_OK)
        Err::errAbort("Region file is missing 'genome-version' meta info.");
    tsv.close();
    return genomeVersion;
}

void CanaryEngine::check_input_file_headers(vector<string> chipTypes, 
                                            vector<string> files, 
                                            string &mapName, 
                                            string &mapVersion, bool force) {
    TsvFile tsv;
    if (tsv.open(files[0]) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + files[0] + "' to read.");
    if (tsv.headersFindNext("map-name", mapName) != TSV_OK)
        Err::errAbort("File '" + files[0] + "' is missing 'map-name' meta info.");
    if (tsv.headersFindNext("map-version", mapVersion) != TSV_OK)
        Err::errAbort("File '" + files[0] + "' is missing 'map-version' meta info.");
    tsv.close();

    if (force)
        return;

    for (int i=1; i < files.size(); i++) {
        TsvFile tsv2;
        string val;
        if (tsv2.open(files[i]) != TSV_OK)
            Err::errAbort("Couldn't open file: '" + files[i] + "' to read.");
        if (tsv2.headersFindNext("map-name", val) != TSV_OK)
            Err::errAbort("File '" + files[i] + 
                          "' is missing 'map-name' meta info.");
        if (val != mapName)
            Err::errAbort("File '" + files[i] + "' has map-name '" + val + 
                          "'. Expecting '" + mapName + "'.");
        if (tsv2.headersFindNext("map-version", val) != TSV_OK)
            Err::errAbort("File '" + files[i] + 
                          "' is missing 'map-version' meta info.");
        if (val != mapVersion)
            Err::errAbort("File '" + files[i] + "' has map-version '" + val + 
                          "'. Expecting '" + mapVersion + "'.");
        tsv2.close();

    }

    for (int i=0; i < files.size(); i++)
        EngineUtil::checkTsvFileChipType(files[i],chipTypes,true);
}

void CanaryEngine::read_region_file(string cnvRegionFile, 
                                    map<string,vector<string> > &cnv_region_map,
                                    vector<string> &cnv_region_names) {
  
    TsvFile regionfile;

    string cnv_region;

    regionfile.bind(0,"cnv_region",&cnv_region, TSV_BIND_REQUIRED);

    if (regionfile.open(cnvRegionFile) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + cnvRegionFile + "' to read.");


    if (regionfile.cname2cidx(0,"cn_probeset_ids") == TSV_ERR_NOTFOUND)
        Err::errAbort("cn_probeset_ids column is missing from regions file");
    if (regionfile.cname2cidx(0,"snp_probeset_ids") == TSV_ERR_NOTFOUND)
        Err::errAbort("snp_probeset_ids column is missing from regions file");

    // Run through the region file
    while (regionfile.nextLevel(0)==TSV_OK) {
        vector<string>cn_probes;
        vector<string>snp_probes;

        // split up probesets
        regionfile.get(0,"cn_probeset_ids",&cn_probes);
        regionfile.get(0,"snp_probeset_ids",&snp_probes);

        cnv_region_names.push_back(cnv_region);

        ///@todo need to check that the entries in regions file are
        //unique given our use of maps.
        for (int k=0; k<cn_probes.size(); k++)
            cnv_region_map[cnv_region].push_back(cn_probes[k]);
        
        for (int k=0; k<snp_probes.size(); k++)
            cnv_region_map[cnv_region].push_back(snp_probes[k]);
    }
    regionfile.close();
}

void CanaryEngine::writeCnvSummaries(CanaryOptions &opts, 
                                     map<string,vector<string> > &cnv_region_map) {
    // Now it is time to read the probe summaries back in

    //OPTIONS
    string probedatafilename = Fs::join(opts.outDir, opts.setAnalysisName + ".allele_average_summary.txt");

    TsvFile probedatafile;
    if (probedatafile.open(probedatafilename) != TSV_OK)
        Err::errAbort("Unable to open probe data file '" + 
                      probedatafilename + "'");
    int ncol = probedatafile.getColumnCount(0);

    // save for the header of the output file
    vector<string>celfnames;
    for (int i=1; i<ncol; i++) {
        string cname = probedatafile.getColumnName(0,i);
        celfnames.push_back(cname);
    }

    map<string,valarray<double> > probe_data_map; 
    while (probedatafile.nextLevel(0)==TSV_OK) {
        string probe_name;
        valarray<double>values(ncol-1);
        probedatafile.get(0,0,probe_name);
        for (int i=1; i<ncol; i++) probedatafile.get(0,i,values[i-1]);
        probe_data_map[probe_name].resize(ncol-1);
        probe_data_map[probe_name] = values;
    }
    probedatafile.close();

    // In the future this should not be the spot to make a temporary
    // file disappear.  Leave it here until the code is reorganised to
    // flow.
    if (!opts.tableOutput)
        Fs::rm(probedatafilename, false);


    // get chip intensities for normalization
    set<string>normalization_psets =
        get_normalization_psets(opts.cnvNormFile);

    Matrix normdata;
    normdata.ReSize(normalization_psets.size(),ncol-1);


    int ncount=0;
    for (set<string>::iterator iter=normalization_psets.begin();
         iter!=normalization_psets.end(); iter++) {
        valarray<double>vec = probe_data_map[*iter];
        for (unsigned int k=0; k<vec.size(); k++)
            normdata.element(ncount,k) = vec[k];
        ncount++;
    }

    // Need to find a consistent way of assigning sizes
    valarray<double>nvalues(normdata.Nrows());
    valarray<double>norm_factor(normdata.Ncols());

    for (int j=0; j<normdata.Ncols(); j++) {
        for (int i=0; i<normdata.Nrows(); i++)
            nvalues[i] = normdata.element(i,j);

        norm_factor[j] = compute_median(nvalues);
    }


    //OPTIONS
    string cnvdatafilename = Fs::join(opts.outDir, opts.setAnalysisName + ".cnv_region_summary.txt");

    TsvFile cnvdatafile;

    // Define the column names from the old file
    cnvdatafile.defineColumn(0,0,"cnv_region_id"); // just to have it defined
    for (int i=1; i<ncol; i++) cnvdatafile.defineColumn(0,i,celfnames[i-1]);
    cnvdatafile.writeTsv_v1(cnvdatafilename);

    // Build and median polish a matrix for each CNV region.
    ///@todo output in same order as the regions file -- rather than the order returned by map
    for (map<string,vector<string> >::iterator iter=cnv_region_map.begin();
         iter!=cnv_region_map.end(); iter++) {
        string cnv_region = iter->first;
        vector<string> probe_names = iter->second;
        
        int N = ncol - 1;
        int P = int(probe_names.size());
	
        
        Matrix the_data;
        the_data.ReSize(N,P);
        
        for (int j=0; j<P; j++) {
            valarray<double> vals = probe_data_map[probe_names[j]];
            for (int i=0; i<N; i++)
                the_data.element(i,j) = log(vals[i]/norm_factor[i]);
        }
        
        valarray<double>cnv_summaries = median_polish(the_data);
        
        cnvdatafile.set(0,0,cnv_region);
        for (int i=0; i<N; i++)
            cnvdatafile.set(0,i+1,exp(cnv_summaries[i]));
        cnvdatafile.writeLevel(0);
    }
    
    cnvdatafile.close();
}

/*
 *  The main function for computing canary.
 */

//----------------------------------------------


void CanaryEngine::mainCanary(CanaryOptions &opts, 
                              map<string,vector<string> > &cnv_region_map, 
                              vector<string> &cnv_region_names) {
    string cnvdatafilename = Fs::join(opts.outDir,opts.setAnalysisName + ".cnv_region_summary.txt");

    Verbose::out(1,"Loading up CNV region signals");
    pair<list<string>,map<string,valarray<double> > >
        intensity_pair = load_broad_intensity_map(cnvdatafilename);

    // The CNV intensity data
    map<string,valarray<double> > intensity_map = intensity_pair.second;

    Verbose::out(1,"Loading up CNV region priors");
    map<string,CanaryPrior> prior_map =
        load_broad_prior_map(opts.canaryPriorFile);

    map<string,vector<int> > call_map;
    map<string,vector<double> > confidence_map;

    // Go through all of the CNV regions
    Verbose::out(1,"Processing CNV regions");
    ///@todo process in same order as the regions file -- rather than the order returned by map
    for (map<string,valarray<double> >::iterator iter=intensity_map.begin();
         iter!=intensity_map.end(); iter++) {
        string cnv_region = iter->first;
//herehere
//		if (cnv_region != "CNP1682") continue; // handy for isolating a CNV
//		if (cnv_region != "CNP12034") continue; // handy for isolating a CNV
//		if (cnv_region != "CNP10276") continue; // handy for isolating a CNV

        CanaryPrior CP = prior_map[cnv_region];
        valarray<double>vals = intensity_map[cnv_region];
        CanaryWithPrior CWP(opts,vals,CP);

        string best_model = CWP.best_fit();

        BroadEstepper1 best_model_BE = CWP.fitted_values(best_model);
        valarray<double>bmeans = best_model_BE.mean();

        vector<double>confidences =
            CWP.fitted_values(best_model).confidences();

        pair<string,BroadEstepper1> model_pair(best_model,best_model_BE);

        BroadEstepper1 newBE = directional_impute(opts,model_pair,CP);

        // tweak the variances

        valarray<double>newbevars = newBE.var();
        //OPTIONS
        double vartweak = newbevars.sum()/newbevars.size();
        double varshare = opts.tune_regularize_variance_factor;

//		newbevars = 0.6*newbevars + 0.4*vartweak;
        newbevars = (1.0 - varshare)*newbevars + varshare*vartweak;
        newBE.var(newbevars);
        newBE.update_prob();

//		confidences = newBE.confidences();
        confidence_map[cnv_region] = newBE.confidences();
        vector<int> assignments = newBE.assignments();
        call_map[cnv_region] = assignments;
    }

    if (opts.ccMDChpOutput || opts.tableOutput) {
        Verbose::out(1,"Generating output file(s).");

        int nregions = 0;
        for (int i=0; i<cnv_region_names.size(); i++)
            if (intensity_map.find(cnv_region_names[i]) != intensity_map.end()) 
                nregions++;

        int max_rlen = 0;
        for (map<string,valarray<double> >::iterator iter=intensity_map.begin();
             iter!=intensity_map.end(); iter++) {
            if (iter->first.size() > max_rlen) 
                max_rlen = (int)iter->first.size();
        }
   

        // What meta info do we want to dump
        vector<string> hNames;
        vector<string> hVals;
        hNames.push_back("apt-engine");                 hVals.push_back("CanaryEngine");
        hNames.push_back("apt-program-name");           hVals.push_back(opts.progName);
        hNames.push_back("apt-command-line");           hVals.push_back(opts.commandLine);
        hNames.push_back("apt-exec-guid");              hVals.push_back(opts.execGuid);
        hNames.push_back("apt-analysis-guid");          hVals.push_back(getOpt("analysis-guid"));
        hNames.push_back("apt-time-str");               hVals.push_back(opts.timeStr);
        hNames.push_back("apt-version");                hVals.push_back(opts.version);
        hNames.push_back("apt-cvs-id");                 hVals.push_back(opts.cvsId);
        hNames.push_back("apt-free-mem");               hVals.push_back(ToStr(opts.freeMemAtStart));
        hNames.push_back("apt-opt-genome-version");     hVals.push_back(opts.genomeVersion);
        hNames.push_back("apt-opt-force");              hVals.push_back(ToStr(opts.force));
        hNames.push_back("apt-opt-chip-type");          hVals.push_back(opts.chipType);
        hNames.push_back("apt-opt-out-dir");            hVals.push_back(opts.outDir);
        hNames.push_back("apt-opt-cdf-file");           hVals.push_back(Fs::basename(opts.cdfFile));
        hNames.push_back("apt-opt-spf-file");           hVals.push_back(Fs::basename(opts.spfFile));
        hNames.push_back("apt-opt-map-file");           hVals.push_back(Fs::basename(opts.cnvMapFile));
        hNames.push_back("apt-opt-region-file");        hVals.push_back(Fs::basename(opts.cnvRegionFile));
        hNames.push_back("apt-opt-canary-norm-file");   hVals.push_back(Fs::basename(opts.cnvNormFile));
        hNames.push_back("apt-opt-canary-prior-file");  hVals.push_back(Fs::basename(opts.canaryPriorFile));
        hNames.push_back("apt-opt-analysis-name");      hVals.push_back(opts.setAnalysisName);
        hNames.push_back("apt-opt-map-name");           hVals.push_back(opts.mapName);
        hNames.push_back("apt-opt-map-version");        hVals.push_back(opts.mapVersion);

        // other headers from regions file meta data
        vector<string> keys;
        vector<string> vals;
        get_meta_tags(opts.cnvRegionFile, keys, vals);
        for (int i=0; i < keys.size(); i++) {
            hNames.push_back("apt-opt-region-file-" + keys[i]);
            hVals.push_back(vals[i]);
        }

        ///@todo these should be the fully enumerated specs, not just what the user asked for -- ie use SelfDoc
        hNames.push_back("apt-opt-analysis-spec");  hVals.push_back(opts.aptSummarizeAnalysis);
        hNames.push_back("apt-opt-canary-spec");    hVals.push_back(opts.aptCanaryAnalysis);

        hNames.push_back("apt-opt-cel-count");      hVals.push_back(ToStr(opts.celFiles.size()));
        for (int i=0; i<opts.celFiles.size(); i++) {
            hNames.push_back("apt-opt-cel-"+ToStr(i+1)); hVals.push_back(Fs::basename(opts.celFiles[i]));
            if (opts.celGuids[i].empty() == false)
                hNames.push_back("apt-opt-cel-guid-"+ToStr(i+1)); hVals.push_back(opts.celGuids[i]);
        }

        if (opts.tableOutput) {
            TsvFile tsvCall,tsvConf;
            int cidx = 0;
            // first column labels the CNP

            tsvCall.defineColumn(0,cidx,"cnv_region_id",affx::TSV_TYPE_UNKNOWN);
            tsvConf.defineColumn(0,cidx,"cnv_region_id",affx::TSV_TYPE_UNKNOWN);
            cidx++;

            // remaining columns labelled by cel files
            for (list<string>::iterator iter=intensity_pair.first.begin();      
                 iter!=intensity_pair.first.end(); iter++){
                tsvCall.defineColumn(0,cidx,*iter,affx::TSV_TYPE_UNKNOWN);
                tsvConf.defineColumn(0,cidx,*iter,affx::TSV_TYPE_UNKNOWN);
                tsvConf.setPrecision(0,cidx,opts.precision);
                cidx++;
            }

            // header info
            for (int i=0; i<hNames.size(); i++) {
                tsvCall.addHeader(hNames[i],hVals[i]);
                tsvConf.addHeader(hNames[i],hVals[i]);
            }

            // @todo - Fix this up so it uses AnalysisInfo rather than writing directly to tsv
            vector<pair<string, string> > metaData = getMetaDataDescription();
            for (int i = 0; i < metaData.size(); i++) {
              pair<string,string> p = metaData[i];
              string key = "affymetrix-application-meta-data-info-" + p.first;
              tsvCall.addHeader(key, p.second);
              tsvConf.addHeader(key, p.second);
            }

            if (tsvCall.writeTsv_v1(Fs::join(opts.outDir, opts.setAnalysisName + ".calls.txt")) != TSV_OK) {
                Err::errAbort("Unable to open call file");
            }

            if (tsvConf.writeTsv_v1(Fs::join(opts.outDir, opts.setAnalysisName + ".confidences.txt")) != TSV_OK) {
                Err::errAbort("Unable to open confidences file");
            }

            ///@todo process in same order as the regions file -- rather than the order returned by map
            for (map<string,valarray<double> >::iterator iter=intensity_map.begin();
                 iter!=intensity_map.end(); iter++) {
                // give the CNP name
                tsvCall.set(0,0,iter->first);
                tsvConf.set(0,0,iter->first);
                
                // give the calls or confidences for the CNP
                for (int i=0; i<opts.celFiles.size(); i++) {
                    tsvCall.set(0,i+1,(unsigned char)call_map[iter->first][i]);
                    tsvConf.set(0,i+1,(float)confidence_map[iter->first][i]);
                }
                
                // write them out
                tsvCall.writeLevel(0);
                tsvConf.writeLevel(0);
                
            }
            tsvCall.close();
            tsvConf.close();
        }

        if (opts.ccMDChpOutput) {
            for (unsigned int k=0; k<opts.chpFiles.size(); k++) {
                try {
                    ParameterNameValueTypeList params;
                    ParameterNameValueType param;

                    CHPMultiDataData CHPMDD(opts.chpFiles[k]);
    
                    CHPMDD.SetArrayType(StringUtils::ConvertMBSToWCS(opts.chipType));
                    CHPMDD.SetAlgName(L"canary");
                    CHPMDD.SetAlgVersion(L"1.0");

                    param.SetName(L"program-name");
                    param.SetValueText(StringUtils::ConvertMBSToWCS(opts.progName));
                    CHPMDD.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
                    param.SetName(L"program-version");
                    param.SetValueText(StringUtils::ConvertMBSToWCS(opts.version));
                    CHPMDD.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
                    param.SetName(L"program-company");
                    param.SetValueText(StringUtils::ConvertMBSToWCS("Affymetrix"));
                    CHPMDD.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

                    CHPMDD.SetFilename(opts.chpFiles[k]);
    
                    CHPMDD.SetEntryCount(CopyNumberVariationMultiDataType, nregions, max_rlen);
    
                    FusionCELData celfile;
                    celfile.SetFileName(opts.celFiles[k].c_str());
                    if (!celfile.ReadHeader()) 
                        Err::errAbort("CEL can't read header\n");
                    GenericData * gdata = celfile.GetGenericData();
        
                    if (gdata)
                        CHPMDD.GetFileHeader()->GetGenericDataHdr()->
                            AddParent(*(gdata->Header().GetGenericDataHdr()));
    
                    for (int i=0; i<hNames.size(); i++) {
                        param.SetName(StringUtils::ConvertMBSToWCS(hNames[i]));
                        param.SetValueText(StringUtils::ConvertMBSToWCS(hVals[i]));
                        params.push_back(param);
                    }

                    param.SetName(StringUtils::ConvertMBSToWCS("apt-opt-cel-file"));
                    param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(opts.celFiles[k])));
                    params.push_back(param);
    
                    CHPMDD.AddAlgParams(params);
    
                    vector<pair<string, string> > metaData = getMetaDataDescription();
                    ParameterNameValueTypeList appMetaInfoList;
                    for (int i = 0; i < metaData.size(); i++) {
                      pair<string,string> p = metaData[i];
                      param.SetName(StringUtils::ConvertMBSToWCS(p.first));
                      param.SetValueText(StringUtils::ConvertMBSToWCS(p.second));
                      appMetaInfoList.push_back(param);
                    }
                    CHPMDD.AddAppMetaInfo(appMetaInfoList);

//                  this is how a generic parameter is addedwhere the key is not altered
//		            CHPMDD.GetGenericData().Header().GetGenericDataHdr()->
//			        AddNameValParam(param);

                    CHPMultiDataFileWriter writer(CHPMDD);
                    writer.SeekToDataSet(CopyNumberVariationMultiDataType);
   
                    for (int i=0; i<cnv_region_names.size(); i++) {
                        if (intensity_map.find(cnv_region_names[i]) != intensity_map.end()) {
                            ProbeSetMultiDataCopyNumberVariationRegionData entry;
                            entry.name = cnv_region_names[i];
                            entry.signal = intensity_map[cnv_region_names[i]][k];
                            entry.call = (unsigned char)call_map[cnv_region_names[i]][k];
                            entry.confidenceScore = (float)confidence_map[cnv_region_names[i]][k];
    
                            writer.WriteEntry(entry);
                        }
                    }
                }
                catch (...) {
                    Err::errAbort("Error writing out chp file " + opts.chpFiles[k]);
                }
            }// end for loop over cel files
	    
            if (!opts.tableOutput)
                Fs::rm(cnvdatafilename, false);

        }// end if (opts.ccMDChpOutput)
    }// end if (output)
}

//----------------------------------------------

/*****************************************************************
 *
 * The Engine for Canary.
 *
 *****************************************************************/

void CanaryEngine::CanaryCompute(CanaryOptions &opts) {
    AnalysisStream * analysis = NULL;
    map<string,vector<string> > cnv_region_map;
    vector<string> cnv_region_names;
    try {
        // Description of the chip
        ChipLayout layout;

        // Probeset names in the cdf or spf file.
        set<string> lib_psets;

        loadCdfLayout_canary(layout, opts, lib_psets);

        // Region file info
        read_region_file(opts.cnvRegionFile, cnv_region_map, cnv_region_names);
        set<string>cnv_psets = get_cnv_psets(cnv_region_map);

        set<string>normalization_psets =
            get_normalization_psets(opts.cnvNormFile);

        vector<string>pset_vec = get_available_psets(cnv_psets,
                                                     normalization_psets,
                                                     lib_psets,opts);

        // Make our analysis from our specifications.
        AnalysisStreamFactory asFactory(QuantMethodFactory::Expression);
        analysis = asFactory.constructExpressionAnalysisStream(
            opts.aptSummarizeAnalysis, layout, opts.stdMethods,
            opts.setAnalysisName);


        // figure out ordering of probes for disk mart
        vector<int> desiredOrder;
        for (int i=0; i<pset_vec.size(); i++) {
            ProbeListPacked pl = layout.getProbeListByName(pset_vec[i]);
            if (!pl.isNull()) {
                for (int pIx = 0; pIx < pl.probe_cnt(); pIx++) {
                    desiredOrder.push_back(pl.get_probeId(pIx));
                }
            }
        }
        IntensityMart* iMart = NULL;
        // RAM cache for data we're computing with.
        int intensityCacheSize = getOptInt("disk-cache") * 1048576;
        if (getOptBool("use-disk")) {
            DiskIntensityMart* diskMart = 
                new DiskIntensityMart(desiredOrder, 
                                      opts.celFiles,
                                      intensityCacheSize,
                                      getOpt("temp-dir"),
                                      "apt-canary.tmp",
                                      true);
            iMart = diskMart;
        }
        else {
            SparseMart* sparseMart = new SparseMart(desiredOrder, 
                                                    opts.celFiles, 
                                                    true);
            iMart = sparseMart;
        }

        CelReader cReader;
        cReader.registerStream(analysis->getChipStreamHead());
        cReader.setFiles(opts.celFiles);
        cReader.registerIntensityMart(iMart);

        analysis->addTextReporter(opts.outDir,opts.precision,false,false);

        // Read CEL files
        cReader.readFiles();
        analysis->prepare(*iMart);
        AnalysisInfo info = analysis->getInfo();
        vector<pair<string, string> > metaData = getMetaDataDescription();
        for (int i = 0; i < metaData.size(); i++) {
          pair<string,string> p = metaData[i];
          info.addClientMetaInfo(p.first, p.second);
        }
        analysis->setInfo(info);
        analysis->addStdHeaders(layout, opts.execGuid,
                                opts.timeStr, opts.commandLine,
                                opts.cvsId + "-" +  opts.version, 
                                analysis->getInfo());

        // Actual analysis computation.
        for (int k=0; k<pset_vec.size(); k++) {
            ProbeListPacked pList = layout.getProbeListByName(pset_vec[k]);
            ProbeSet* ps = ProbeListFactory::asProbeSet(pList);
            ProbeSetGroup psGroup(ps); // psGroup will free ps
            analysis->doAnalysis(psGroup, *iMart, true);
        }
    
        // Tell analyses that we are done computing.
        analysis->finish();
        Freez(iMart);
    } // end try
    catch(Except &e) {
        Freez(analysis);
        Err::errAbort(e.what());
    }
    catch(const BaseException &e) {
        Freez(analysis);
        Err::errAbort("newmat issue: " + ToStr(e.what()));
    }
    catch(CalvinException &ce) {
        Freez(analysis);
        Err::errAbort("Affymetrix GeneChip Command Console library has thrown an exception. Description: '" + 
                      StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
    }
    catch(const std::bad_alloc &e) {
        Freez(analysis);
        Err::errAbort("Ran out of memory.");
    }
    catch(const std::exception &e) {
        Freez(analysis);
        Err::errAbort("standard exception: " + ToStr(e.what()));
    }
    catch(...) {
        Freez(analysis);
        Err::errAbort("Uncaught Exception.");
    }
    /* Cleanup. */
    Freez(analysis);

    // Reading and writing temporary files here because of the merging
    averageAlleles(opts);
    writeCnvSummaries(opts,cnv_region_map);
	
    // Check to see if the output cnv file names are there
    if (!opts.chpFiles.size()) {
        for (unsigned int k=0; k<opts.celFiles.size(); k++) {
            string chpfname = Fs::basename(opts.celFiles[k]);
            chpfname = swap_suffix(chpfname,"." + opts.setAnalysisName + 
                                   ".cnvchp");
            chpfname = Fs::join(opts.outDir,chpfname);
            opts.chpFiles.push_back(chpfname);
        }
    }

    mainCanary(opts,cnv_region_map,cnv_region_names);
}
