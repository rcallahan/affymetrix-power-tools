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
 * @file   multiChannelExample.cpp
 * @author Alan Williams
 * @date   Fri Jan  4 13:53:19 PST 2008
 * 
 * @brief A simple example of multi channel analysis
 * Please do not use this example for real life data
 * analysis.
 */

//
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/BioTypes.h"
#include "chipstream/CelReader.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/FileGenders.h"
#include "chipstream/GenoSeedTxt.h"
#include "chipstream/MultiQuantMethodListener.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/QuantMethodGTypeReport.h"
#include "chipstream/QuantMethodReportListener.h"
#include "chipstream/SpecialSnps.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace affx;

using namespace std;

class AOptions {
public:
  bool help;
  int verbose;
  string algoParams;
  string outDir;
  string exprAnalysisString;
  string chrXSnps;
  string specialSnps;
  string readModelsBrlmmp;
  bool writeModels;
  bool summaries;
  bool doFeatureEffects;
  bool doResiduals;
  bool selectProbes;
  string genderFileIn;
  string quantMethod;
  std::string genoTypesFile;     ///< File to input seed genotypes
  vector<string> sketchInFile;
  bool writeSketch;

  string spfFile;
  string cdfFile;

  int channelCount;
  vector<string> sampleNames;
  vector<vector<string> > celFiles; ///< first dim channel, second is samples

  string channelFile;
};

void define_options(PgOptions* opts)
{
  opts->setUsage("multiChannelExample - example of feeding multiple channels into a\n"
                 "ChipStream for genotype calling.\n"
                 "\n"
                 "usage:\n"
                 "   ./multiChannelExample -p CM=1.bins=100.K=2.SB=0.003.MS=0.05 \\\n"
                 "       -o results \\\n"
                 "       --read-models-brlmmp GenomeWideSNP_5.models \\\n"
                 "       --chrX-snps GenomeWideSNP_5.chrx \\\n"
                 "       --read-genders genders.txt\n");
                 
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("", "explain", PgOpt::STRING_OPT,
                     "Display documentation for algo method parameters.",
                     "");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "File listing CEL files to process. There must be a 'sample_name' column, and 2 or more 'channel_*'. Current channels are ordered (starting with 0) by the order in this file.",
                     "");
  opts->defineOption("d", "cdf-file", PgOpt::STRING_OPT,
                     "CDF file providing chip layout information.",
                     "");
  opts->defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "SPF file providing chip layout information.",
                     "");
  opts->defineOption("", "channel-file", PgOpt::STRING_OPT,
                     "Channel file providing info on what alleles are in what channel for each SNP. There should be a row for each SNP. There must be a 'probeset_id', 'allele_A_channel', and an 'allele_B_channel' columns. The value corresponds to the channels in the --cel-files specified file. Valid values are integers >= 0.",
                     "");
  opts->defineOption("p", "algo-params", PgOpt::STRING_OPT,
                     "A string specifying algo parameters. See --explain option for acceptable parameters.",
                     "");
  opts->defineOption("a", "expr-analysis", PgOpt::STRING_OPT,
                     "Expression analysis string.",
                     "");
  opts->defineOption("o", "out-dir", PgOpt::STRING_OPT,
                     "Directory to write result files into.",
                     ".");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("", "chrX-snps", PgOpt::STRING_OPT,
                     "File containing snps on chrX (non-pseudoautosomal region).",
                     "");
  opts->defineOption("", "special-snps",PgOpt::STRING_OPT,
                     "File containing all snps of unusual copy (chrX,mito,Y)",
                     "");
  opts->defineOption("", "write-models", PgOpt::BOOL_OPT,
                     "Should we write snp specific models out for analysis? [experimental]",
                     "false");
  opts->defineOption("", "summaries", PgOpt::BOOL_OPT,
                     "Should we write out the allele summaries",
                     "false");
  opts->defineOption("", "feature-effects", PgOpt::BOOL_OPT,
                     "Output feature effects when available.",
                     "false");
  opts->defineOption("", "feature-details", PgOpt::BOOL_OPT,
                     "Output probe by chip specific details (often residuals) when available.",
                     "false");
  opts->defineOption("", "read-models-brlmmp", PgOpt::STRING_OPT,
                     "File to read precomputed BRLMM-P snp specific models from.",
                     "");
  opts->defineOption("", "read-genders", PgOpt::STRING_OPT,
                     "Explicitly read genders from a file.",
                     "");

  opts->defineOption("", "genotypes", PgOpt::STRING_OPT,
                     "File to read seed genotypes from.",
                     "");
  opts->defineOption("","select-probes",PgOpt::BOOL_OPT,
					  "evaluate probe/probe-set quality",
			"false");
  opts->defOptMult("", "target-sketch", PgOpt::STRING_OPT,
                    "File specifying a target distribution to use for quantile normalization.",
                    "");
  opts->defineOption("", "write-sketch", PgOpt::BOOL_OPT,
                    "Write the quantile normalization distribution (or sketch) to a file for reuse with target-sketch option.",
                    "false");
}

void fillInOptions(PgOptions *opts, AOptions &o) {
  PgOpt *asOpt = NULL;
  o.help = opts->getBool("help");
  o.verbose = opts->getInt("verbose");
  o.outDir = opts->get("out-dir");
  Util::chompLastIfSep(o.outDir);

  o.spfFile = opts->get("spf-file");
  o.cdfFile = opts->get("cdf-file");

  o.channelFile = opts->get("channel-file");

  o.quantMethod = "brlmm-p";
  o.algoParams = opts->get("algo-params");
  if(o.algoParams == "")
    o.algoParams = "CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.000001.wobble=0.05.MS=0.05";
  o.exprAnalysisString = opts->get("expr-analysis");
  if(o.exprAnalysisString == "")
    o.exprAnalysisString = "med-norm.target=1000,pm-only,med-polish.expon=true";

  o.specialSnps = opts->get("special-snps");
  o.chrXSnps = opts->get("chrX-snps");
  o.readModelsBrlmmp=opts->get("read-models-brlmmp");
  o.writeModels = opts->getBool("write-models");
  o.summaries = opts->getBool("summaries");
  o.doFeatureEffects = opts->getBool("feature-effects");
  o.doResiduals = opts->getBool("feature-details");
  o.genderFileIn = opts->get("read-genders");
  o.genoTypesFile = opts->get("genotypes");
  o.selectProbes = opts->getBool("select-probes");
  o.genoTypesFile = opts->get("genotypes");
  o.writeSketch = opts->getBool("write-sketch");
  asOpt = opts->findOpt("target-sketch");
  if (asOpt!=NULL) {
    asOpt->push_user_values_into(o.sketchInFile);
  }

  /* Sanity checks... */
  if(o.outDir == "") 
    Err::errAbort("Must specify an output directory.");

  /* Make our directory. */
  if(!Fs::isWriteableDir(o.outDir.c_str())) {
    if(Fs::mkdirPath(o.outDir, false) != APT_OK) {
      APT_ERR_ABORT("Can't make or write to directory: " + ToStr(o.outDir));
    }
  }

  /* load up cel file info */
  if (opts->get("cel-files")!="") {
    affx::TsvFile tsv;
#ifdef WIN32
    tsv.m_optEscapeOk = false;
#endif
    std::string celFiles = opts->get("cel-files");
    std::string sample;
    tsv.bind(0, "sample_name", &sample, TSV_BIND_REQUIRED);
    if(tsv.open(celFiles) != TSV_OK) {
      Err::errAbort("Couldn't open cell-files file: '" + celFiles + "'");
    }
    int numCol = tsv.getColumnCount(0);
    vector<int> colIndexes;
    for(int cidx=0; cidx < numCol; cidx++) {
        string name;
        tsv.cidx2cname(0,cidx,name);
        if(name.find("channel_") == 0) {
            ///@todo we order the channels by order seen, not the label
            ///      would be nice to track the label and associate via
            ///      label in lib/panel file
            colIndexes.push_back(cidx);
        }
    }
    o.channelCount = colIndexes.size();
    o.celFiles.resize(o.channelCount);
    tsv.rewind();
    while(tsv.nextLevel(0) == TSV_OK) {
        o.sampleNames.push_back(sample);
        for(int i=0; i<colIndexes.size(); i++) {
            string val;
            if(tsv.get(0,i+1,val)==TSV_OK) {
                o.celFiles[i].push_back(val);
            } else
                Err::errAbort("Unable to read column " + ToStr(i+1) + " from " + celFiles);
        }
    }
    tsv.close();
    Verbose::out(1, "Read " + ToStr(o.celFiles[0].size()) + " sample(s) and " + ToStr(o.celFiles.size()) + " channel(s) from: " + Fs::basename(celFiles));
  } else {
      Err::errAbort("Must provide a list of cel files for each channel.");
  }
}

void callGenotypes(AOptions &o) {

  /* Get the layout for the chip. */
  ChipLayout layout;
  if(o.cdfFile != "") {
    Verbose::out(1,"Opening cdf file: " + ToStr(o.cdfFile));
    layout.openCdfAll(o.cdfFile);
  } else if(o.spfFile != "") {
    Verbose::out(1,"Opening spf file: " + ToStr(o.spfFile));
    layout.openSpfAll(o.spfFile);
  } else {
    Err::errAbort("Must provide a CDF or SPF file.");
  }

  /* Get the channel file */
  map<string,vector<int> > sampleAlleleMap; ///@todo not particularly space efficient
  if (o.channelFile!="") {
    affx::TsvFile tsv;
#ifdef WIN32
    tsv.m_optEscapeOk = false;
#endif
    std::string probeset;
    int alleleA, alleleB;
    tsv.bind(0, "probeset_name", &probeset, TSV_BIND_REQUIRED);
    tsv.bind(0, "allele_A_channel", &alleleA, TSV_BIND_REQUIRED);
    tsv.bind(0, "allele_B_channel", &alleleB, TSV_BIND_REQUIRED);
    if(tsv.open(o.channelFile) != TSV_OK) {
      Err::errAbort("Couldn't open channel file: '" + o.channelFile + "'");
    }
    tsv.rewind();
    while(tsv.nextLevel(0) == TSV_OK) {
        sampleAlleleMap[probeset].push_back(alleleA);
        sampleAlleleMap[probeset].push_back(alleleB);
    }
    tsv.close();
  } else {
      Verbose::out(1,"No channel file. Assuming allele A on channel A and allele B on channel B.");
  }

  /* Construct the analysis streams for each factory */
  AnalysisStreamFactory asFactory; 
  vector<AnalysisStream *> streams;
  Verbose::out(1,"Allele estimates computed using: " + ToStr(o.exprAnalysisString));
  Verbose::progressBegin(1,"Construction ChipStreams for computing allele estimates on each channel",
    o.channelCount, 0, o.channelCount);
  for(int channel = 0; channel < o.channelCount; channel++) {
    Verbose::progressStep(1);
    // do we need to write out the sketch
    if(o.writeSketch) {
      string prefix = Fs::join(o.outDir,"ch" + ToStr(channel) + ".");
      asFactory.setWriteSketchDir(prefix);
    }
    // do we have a sketch to read in
    if(o.sketchInFile.size() == 1) {
      Verbose::out(1, "Opening target normalization file for channel " + ToStr(channel) + 
              ": " + Fs::basename(o.sketchInFile[0].c_str()));
      asFactory.readTargetSketchFromFile(o.sketchInFile[0].c_str());
    }else if(o.sketchInFile.size() > 1) {
      if(o.sketchInFile.size() != o.channelCount) {
          Verbose::out(1,"When specifying more than 1 sketch, number of sketches must match number of channels");
      }
      Verbose::out(1, "Opening target normalization file for channel " + ToStr(channel) + 
              ": " + Fs::basename(o.sketchInFile[channel].c_str()));
      asFactory.readTargetSketchFromFile(o.sketchInFile[channel].c_str());
    }
    streams.push_back(asFactory.constructAnalysisStream(o.exprAnalysisString.c_str(),layout));
  }
  Verbose::progressEnd(1, "Done.");
  
  /* Setup Intensity Mart for each channel */
  Verbose::progressBegin(1,"Setting up intensity marts for each channel",
    streams.size(), 0, streams.size());
  vector<DiskIntensityMart *> iMarts;
  for(int i=0; i<streams.size(); i++) {
      Verbose::progressStep(1);
      vector<int> desiredOrder;
      for(unsigned int psIx = 0; psIx < layout.getProbeSetCount(); psIx++) {
          ProbeListPacked pl = layout.getProbeListAtIndex(psIx);
          for(int pIx = 0; pIx < pl.probe_cnt(); pIx++) {
              desiredOrder.push_back(pl.get_probeId(pIx));
          }
      }
      DiskIntensityMart * ptr = new DiskIntensityMart(
                desiredOrder,
                o.celFiles[i], 
                50 * 1048576, 
                o.outDir, "multiChannelExample.tmp.ch" + ToStr(i), true);
      iMarts.push_back(ptr);
  }
  Verbose::progressEnd(1, "Done.");

  /* Setup CelReaders for each Channel */
  Verbose::progressBegin(1,"Setting up CEL readers for each channel",
    streams.size(), 0, streams.size());
  vector<CelReader> cReaders(streams.size());
  for(int i=0; i<streams.size(); i++) {
    Verbose::progressStep(1);
    /* As cel files are read the chipstream wants to see them and learn parameters. */
    cReaders[i].registerStream(streams[i]->getChipStreamHead());

    /* As cel files are read the data cache needs to store intensities. */
    cReaders[i].registerIntensityMart(iMarts[i]);
  }
  Verbose::progressEnd(1, "Done.");

  /* Setup the MultiChannelQuantMethod */
  MultiQuantMethodListener multiQuant;
  for(int i=0; i<streams.size(); i++) {
      // Because the stream owns the reporter we
      // use the listener rather than adding the multi quant directly
      QuantMethodReportListener *listener = new QuantMethodReportListener();
      listener->registerListener(&multiQuant);
      streams[i]->addReporter(listener);
  }

  // Create Quant Method
  Verbose::out(1,"Creating QuantMethod for genotyping from: " + o.quantMethod + "." + o.algoParams);
  QuantMethodFactory quantMethodFactory(QuantMethodFactory::GenoType);
  string callingString = o.quantMethod + "." + o.algoParams;
  QuantGTypeMethod *qMethod = 
      quantMethodFactory.quantGTypeMethodForString(callingString, layout, QuantMethodFactory::GenoType);
  QuantLabelZ *qBrlmmp = dynamic_cast<QuantLabelZ *>(qMethod);

  // setup genders
  vector<Gender> genders;
  if(o.genderFileIn != "") {
      Verbose::out(1,"Reading genders from " + o.genderFileIn);
      FileGenders genderReader(o.genderFileIn,o.sampleNames);
      genders = genderReader.getGenders();
  } else {
      Verbose::out(1,"No gender file. Assuming all samples are unknown.");
      for(int i = 0; i < o.sampleNames.size(); i++)
          genders.push_back(UnknownGender);
  }
  qBrlmmp->setGenders(genders);

  // setup special SNPs
  SpecialSnpMap specialSnps;
  map<string,bool> haploidSnps;
  if(o.specialSnps != "") {
      Verbose::out(1,"Reading in special SNPs from o.specialSnps");
      string chip, specialType;
      readSpecialSnpsFile(&specialSnps, &chip, &specialType, o.specialSnps);
      haploidFromSpecial(specialSnps,haploidSnps);
  } else if(o.chrXSnps != "") {
      Verbose::out(1,"Reading in special SNPs from o.chrXSnps");
      fillInHaploidSnps(o.chrXSnps, haploidSnps, NULL, true);
      specialFromHaploid(specialSnps, haploidSnps);
  } else {
      Verbose::out(1,"No special SNPs file provided. Assuming everything is bi-allelic.");
  }

  if(!specialSnps.empty())
      qBrlmmp->setSpecialSnps(specialSnps);
  else if(!haploidSnps.empty())
      qBrlmmp->setHaploidSnps(haploidSnps);

  // read in SNP models
  if(o.readModelsBrlmmp!="") {
      Verbose::out(1,"Reading in SNP models from " + o.readModelsBrlmmp);
      qBrlmmp->readSnpPriorMap(o.readModelsBrlmmp.c_str());
  } else {
      Verbose::out(1,"No SNP models provided. Assuming generic model.");
  }

  // write SNP models?
  if(o.writeModels) {
    string outfile = Fs::join(o.outDir,o.quantMethod + ".snp-posteriors");
    qBrlmmp->writeSnpPosteriorTsv(outfile,affx::TsvReport::FMT_TSV);
  }
  if (o.selectProbes){
	  string outfile = Fs::join(o.outDir,o.quantMethod + ".select-probes");
    qBrlmmp->writeSnpProbeTsv(outfile,affx::TsvReport::FMT_TSV);
  }
  
  // read in seed genotypes
  std::map<const char *,std::vector<GType>, Util::ltstr > knownGenotypes;
  if(!o.genoTypesFile.empty()) {
    Verbose::out(1,"Reading in SNP genotype calls from " + o.genoTypesFile);
    std::vector<std::string> probesets;
    GenoSeedTxt seeds(o.genoTypesFile.c_str(), layout, o.sampleNames, haploidSnps);
    seeds.fillProbesets(probesets);

    for(int i=0; i<probesets.size(); i++) {
      const char *name = Util::cloneString(probesets[i].c_str());
      vector<GType> callVec = seeds.getGenoCalls(name);
      knownGenotypes[name] = callVec;
      //Verbose::out(1,"Setting known genotypes for: " + ToStr(name));
    }
    qBrlmmp->setKnownGenoTypes(&knownGenotypes);
  }

  // Dummy Data Mart for Reporters
  // Used to spoof sample names
  DiskIntensityMart *dummyMart = iMarts[0];
    
  // Create Report
  Verbose::out(1,"Creating main genotype reporter.");
  QuantMethodGTypeReport *report = new QuantMethodGTypeReport();
  report->setDirPath(o.outDir);
  report->setFileprefix(o.quantMethod);
  report->setFormat(affx::TsvReport::FMT_TSV);
  report->setFilename("QuantMethodGTypeReport-DEBUG"); // filename not used -- for debugging
  report->setPrecision(5);
  
  string qMethodGuid = affxutil::Guid::GenerateNewGuid();
  string execGuid = affxutil::Guid::GenerateNewGuid();
  string timeStr = ""; ///@todo set time string
  string commandLine = ""; ///@todo set command line string
  string execVersion = ""; ///@todo set command line string
  AnalysisInfo info; ///@todo populate analysis info object
  
  report->addStdHeaders((QuantMethodReport *)report, 
                        execGuid,
                        qMethodGuid,
                        timeStr,
                        commandLine,
                        //o.sampleNames,
                        execVersion, 
                        info);
  //
  report->prepare(*qMethod, *dummyMart);

  // Create allele summaries report
  QuantMethodExprReport *summaryReport = NULL;
  QuantExprMethod *eMethod = NULL;
  if(o.summaries) {
        eMethod = dynamic_cast<QuantExprMethod *>(streams[0]->getQuantMethod());
        if(eMethod == NULL)
            Err::errAbort("Must use QuantExprMethod to report summaries");
//        string reportPrefix = o.outDir + PATH_SEPARATOR + o.quantMethod + "." + eMethod->getType();
//        summaryReport = new QuantMethodExprReport(reportPrefix.c_str(),
//                                                  *eMethod, 
//                                                  true,
//                                                  o.doFeatureEffects,
//                                                  o.doResiduals,
//                                                  5);
        summaryReport = new QuantMethodExprReport(info.getNumCols());
        summaryReport->setIsHeaderBuffer(1);
        summaryReport->setDirPath(o.outDir);
        summaryReport->setFileprefix(o.quantMethod + "." + eMethod->getType());
        summaryReport->setFilename("QuantMethodExprReport-DEBUG"); // filename not used -- for debugging
        summaryReport->setFormat(affx::TsvReport::FMT_TSV);
        //
        summaryReport->m_qMethod=eMethod;
        summaryReport->setPrecision(5);
        summaryReport->m_DoSummary=true;
        summaryReport->m_DoFeatureEffects=true;
        summaryReport->m_DoResiduals=true;

        //summaryReport->writeHeader(summaryReport, execGuid.c_str(), qMethodGuid.c_str(),
        //                          timeStr.c_str(), commandLine.c_str(),
        //                          o.sampleNames, execVersion.c_str(), info);

        report->addStdHeaders(report,
                              execGuid,
                              qMethodGuid,
                              timeStr,
                              commandLine,
                              execVersion,
                              info);

        summaryReport->prepare(*eMethod, *dummyMart);
  }

  // Dummies needed by the reporter
  PmOnlyAdjust pmAdjust;
  vector<ChipStream *> iTrans;

  /* register the cel files to open. */
  Verbose::out(1,"Registering CEL files for each channel.");
  for(int chIx=0; chIx<streams.size(); chIx++) {
    vector<string> fileNames;
    for(int celIx = 0; celIx < o.celFiles[chIx].size(); celIx++)
        fileNames.push_back(o.celFiles[chIx][celIx]);
    cReaders[chIx].setFiles(fileNames);
  }

  // Read in cel files
  for(int i=0; i<cReaders.size(); i++) {
      cReaders[i].readFiles();
  }


  /* Loop through all of our probe sets and do an expression quantification for each one. */
  Verbose::progressBegin(1,"Processing probesets",
    40, (int)(layout.getProbeSetCount()/40), layout.getProbeSetCount());
  for(unsigned int psIx = 0; psIx < layout.getProbeSetCount(); psIx++) {
    Verbose::progressStep(1);
    ProbeListPacked pList = layout.getProbeListAtIndex(psIx);
    ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
    ProbeSetGroup group(ps);     // psGroup should delete the memory for ps...
    for(int i=0; i<streams.size(); i++){
        streams[i]->doAnalysis(group, *iMarts[i], true);
    }

    ///@todo much of this should probably be in a MultiQuantLabelZ
    ///      rather than using QuantLabelZ directly. Accessing
    ///      m_Results directly is ugly as well.

    // set probeset name for report
    qBrlmmp->setProbeSetName(ps->name);

    // Figure out the channels for this SNP
    ///@todo default should be global (ie if we have a channel file, then all probesets must be found in it. only default to 0/1 when there is no channel file)
    int alleleAIndex = 0;
    int alleleBIndex = 1;
    if(sampleAlleleMap.find(ps->name) != sampleAlleleMap.end()) {
        alleleAIndex = sampleAlleleMap[ps->name][0];
        alleleBIndex = sampleAlleleMap[ps->name][1];
    } else {
        Verbose::out(1,"Unable to find probeset " + ToStr(ps->name) + " in channel file. Assuming channels 1 and 2.");
    }

    qBrlmmp->setAlleleValues(multiQuant.m_Results[alleleAIndex], multiQuant.m_Results[alleleBIndex]);

    // dump allele summaries
    if(o.summaries) {
        // Allele A
        group.name = Util::cloneString((ToStr(ps->name) + "-A").c_str());
        summaryReport->report(group, *(streams[alleleAIndex]->getQuantMethod()), *(iMarts[alleleAIndex]), iTrans, pmAdjust);
        // Allele B
        group.name = Util::cloneString((ToStr(ps->name) + "-B").c_str());
        summaryReport->report(group, *(streams[alleleBIndex]->getQuantMethod()), *(iMarts[alleleBIndex]), iTrans, pmAdjust);
    }

    // transform the data as specified
    qBrlmmp->transform();
   
    // load up seed genotypes
    qBrlmmp->loadInitialCallsForSnp(ps->name);

    // make the call
    qBrlmmp->computeEstimate();
    
    // report results
    report->report(group, *qMethod, *iMarts[0], iTrans, pmAdjust);

    // clear out results for the current probeset
    multiQuant.m_Results.clear();
  }
  Verbose::progressEnd(1, "Done.");

  // close the report
  report->finish(*qMethod);
  if(summaryReport != NULL) 
    summaryReport->finish(*eMethod);

  // delete the QuantMethod
  delete qBrlmmp;
  
  // Free up memory
  for(int i = 0; i < streams.size(); i++)
    delete streams[i];
  for(int i = 0; i < iMarts.size(); i++)
    delete iMarts[i];
}

void printParams(string algo) {
    cout << endl;
    cout << "+--------------------------------------------------------+" << endl;
    cout << "| Algo Parameters                                        |" << endl;
    cout << "| Not all are relevant when feeding in allele summaries  |" << endl;
    cout << "+--------------------------------------------------------+" << endl;
    cout << endl;

    QuantMethodFactory gFactory(QuantMethodFactory::GenoType);
    vector<vector<SelfDoc> > docs;
    docs.push_back(gFactory.getDocs());

    for(unsigned int i = 0; i < docs.size(); i++) {
      for(unsigned int j = 0; j < docs[i].size(); j++) {
        if(docs[i][j].getDocName() == algo) {
            SelfDoc::printExplanation(docs[i][j], cout);
        }
      }
    }
}

/** Everbody's favorite function... */
int main(int argc, const char *argv[]) {
  const string version ("NON-OFFICIAL-RELEASE");
  const string cvsId ("$Id: multiChannelExample.cpp,v 1.34 2009-10-30 07:20:38 csugne Exp $");
  ofstream logOut;
  string logName;

  AOptions o;
  PgOptions *pgOpts = NULL;
  pgOpts = new PgOptions();
  define_options(pgOpts);
  pgOpts->parseArgv(argv);

  if(pgOpts->getBool("help") == true || argc == 1)
  {
    pgOpts->usage();
    cout << "version:" << endl;
    cout << "   " << version << endl;
    cout << "   " << cvsId << endl;

  }
  else if(pgOpts->get("explain") != "")
    printParams(pgOpts->get("explain"));
  else if(pgOpts->getBool("version"))
    cout << "version: " << version << " " << cvsId << endl;
  else {
    fillInOptions(pgOpts, o);
    logName = Fs::join(o.outDir,"multiChannelExample.log");
    Fs::mustOpenToWrite(logOut, logName);
    LogStream log(3, &logOut);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);
    Verbose::setLevel(o.verbose);
    callGenotypes(o);
  }
  return 0;
}
    
