////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "chipstream/apt-summary-genotype/SummaryGenotypeEngine.h"
//
#include "chipstream/ChipLayout.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/FileGenders.h"
#include "chipstream/GenoSeedTxt.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/QuantMethodGTypeReport.h"
#include "chipstream/SelfDoc.h"
#include "chipstream/SparseMart.h"
#include "chipstream/SpecialSnps.h"
//
#include "util/Fs.h"
#include "util/PgOptions.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//
#ifndef WIN32
#include <unistd.h>
#endif /* WIN32 */

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;


SgIoHelper::SgIoHelper(const std::string &summariesFile,
                   const std::string &summariesFileA, const std::string &summariesFileB,
                   const std::string &file5,
                   const std::string &file5Group) {
  f5SumGroup= NULL;
  f5SumTsv = NULL;
  if(summariesFile != "") {
    string label;
    if(summaries.open(summariesFile) != TSV_OK)
      Err::errAbort("Couldn't open file " + ToStr(summariesFile) + " to read.");
    summaries.cidx2cname(0,0,label);
    if(label != "probeset_id")
      Err::errAbort("probeset_id must be the first column in the summaries file");
    int numCol = summaries.getColumnCount(0);
    for(int cidx=1; cidx < numCol; cidx++) {
      string name;
      summaries.cidx2cname(0,cidx,name);
      // We need *both* sets of column names.
      m_ColNames.push_back(name);
    }
    
  } 
  else if(summariesFileA != "" && summariesFileB != "") {
    string label;
    if(summariesA.open(summariesFileA) != TSV_OK)
      Err::errAbort("Couldn't open file " + ToStr(summariesFileA) + " to read.");
    summariesA.cidx2cname(0,0,label);
    if(label != "probeset_id")
      Err::errAbort("probeset_id must be the first column in the summaries file");
    if(summariesB.open(summariesFileB) != TSV_OK)
      Err::errAbort("Couldn't open file " + ToStr(summariesFileB) + " to read.");
    summariesB.cidx2cname(0,0,label);
    if(label != "probeset_id")
      Err::errAbort("probeset_id must be the first column in the summaries file");
    int numCol = summariesA.getColumnCount(0);
    int numColB = summariesB.getColumnCount(0);
    if(numCol != numColB)
      Err::errAbort("allele A and allele B signal files do not have the same number of columns");
      
    //
    for(int cidx=1; cidx < numCol; cidx++) {
      string name;
      string nameB;
      summariesA.cidx2cname(0,cidx,name);
      summariesB.cidx2cname(0,cidx,nameB);
      if(name != nameB)
        Err::errAbort("Sample names do not match between allele A and B signal files.");
      // We need *both* sets of column names.
      m_ColNames.push_back(name);
    }
      
  } 
  else if(file5 != "" && file5Group != "") {
    int status = f5Summaries.open(file5, affx::FILE5_OPEN_RO);
    if(status != FILE5_OK) {
      Err::errAbort("Can't open: " + file5);
    }
    f5SumGroup = f5Summaries.openGroup(file5Group, affx::FILE5_OPEN);
    f5SumTsv = f5SumGroup->openTsv("summaries", affx::FILE5_OPEN_RO);
    int columnCount = f5SumTsv->getColumnCount(0);
    string name;
    for(int i = 1; i < columnCount; i++) {
      f5SumTsv->getColumnName(0,i,&name);
      m_ColNames.push_back(name);
    }
  }
  else {
    Err::errAbort("Must provide either an interleaved summary file or separate allele A and B summary files.");
  }
}

SgIoHelper::~SgIoHelper() {
  if(f5SumTsv != NULL) {
    f5SumTsv->close();
    Freez(f5SumTsv);
  }
  if(f5SumGroup != NULL) {
    f5SumGroup->close();
    Freez(f5SumGroup);
  }
  f5Summaries.close();
}


void SgIoHelper::getLine(string &label, vector<double> &vals, TsvFile &summaries, vector<string> &colNames) {
  vals.clear();
  for(int i=0; i<colNames.size(); i++) {
    double val;
    if(summaries.get(0,i+1,val)==TSV_OK)
      vals.push_back(val);
    else
        Err::errAbort("Unable to read column " + ToStr(i+1) + " from summaries file.");
  }
  
  if(summaries.get(0,0,label)!=TSV_OK)
    Err::errAbort("Unable ro read label from summaries file!");
}

void SgIoHelper::getFile5Line(string &label, vector<double> &vals, File5_Tsv *summaries, vector<string> &colNames) {
  vals.clear();
  for(int i=0; i<colNames.size(); i++) {
    double val;
    if(summaries->get(0,i+1,&val)==FILE5_OK)
      vals.push_back(val);
    else
      Err::errAbort("Unable to read column " + ToStr(i+1) + " from summaries file.");
  }
  
  if(summaries->get(0,0,&label) != FILE5_OK)
    Err::errAbort("Unable ro read label from summaries file!");
}

bool SgIoHelper::getNextMarkerTsvSingle(vector<string> &colNames, vector<double> &aVals, 
  vector<double> &bVals, string &probesetName) {
  if(summaries.nextLevel(0)==TSV_OK) {
    string label, alleleName;
    
    // parse the "A" allele
    string aProbesetName, aProbesetAllele;
    getLine(label, aVals, summaries, colNames);
    aProbesetName = label.substr(0, label.rfind('-'));
    alleleName = label.substr(label.rfind('-')+1, label.size());
    if(alleleName != "A")
      Err::errAbort("Expecting to see summaries for allele A first for " + aProbesetName);
    
    // parse the "B" allele
    if(!(summaries.nextLevel(0)==TSV_OK)) {
      Err::errAbort("Unable to read Allele B for probeset " + aProbesetName);
    }
    string bProbesetName, bProbesetAllele;
    getLine(label, bVals, summaries, colNames);
    bProbesetName = label.substr(0, label.rfind('-'));
    alleleName = label.substr(label.rfind('-')+1, label.size());
    if(alleleName != "B")
      Err::errAbort("Expecting to see summaries for allele B second for " + aProbesetName);
    
    // check that we have the right pair
    if(aProbesetName != bProbesetName)
      Err::errAbort("Allele summaries do not pair up: " + aProbesetName + " != " + bProbesetName);
    probesetName = aProbesetName;
    
    return true;
  } 
  return false;
}

bool SgIoHelper::getNextMarkerFile5(vector<string> &colNames, vector<double> &aVals, 
                                  vector<double> &bVals, string &probesetName) {
  if(f5SumTsv->nextLevel(0)==FILE5_OK) {
    string label, alleleName;
    
    // parse the "A" allele
    string aProbesetName, aProbesetAllele;
    getFile5Line(label, aVals, f5SumTsv, colNames);
    aProbesetName = label.substr(0, label.rfind('-'));
    alleleName = label.substr(label.rfind('-')+1, label.size());
    if(alleleName != "A")
      Err::errAbort("Expecting to see summaries for allele A first for " + aProbesetName);
      
    // parse the "B" allele
    if(!(f5SumTsv->nextLevel(0)==FILE5_OK)) {
      Err::errAbort("Unable to read Allele B for probeset " + aProbesetName);
    }
    string bProbesetName, bProbesetAllele;
    getFile5Line(label, bVals, f5SumTsv, colNames);
    bProbesetName = label.substr(0, label.rfind('-'));
    alleleName = label.substr(label.rfind('-')+1, label.size());
    if(alleleName != "B")
      Err::errAbort("Expecting to see summaries for allele B second for " + aProbesetName);
      
    // check that we have the right pair
    if(aProbesetName != bProbesetName)
      Err::errAbort("Allele summaries do not pair up: " + aProbesetName + " != " + bProbesetName);
    probesetName = aProbesetName;
      
    return true;
  } 
  return false;
}

bool SgIoHelper::getNextMarkerTwoTsv(vector<string> &colNames, vector<double> &aVals, 
                                   vector<double> &bVals, string &probesetName) {
  if(summariesA.nextLevel(0)!=TSV_OK)
    return false;
  if(summariesB.nextLevel(0)!=TSV_OK)
    return false;
    
  string labelA, labelB;
  getLine(labelA, aVals, summariesA, colNames);
  getLine(labelB, bVals, summariesB, colNames);
  if(labelA != labelB)
    Err::errAbort("Probeset mismatch between allele signal files.");
  probesetName = labelA;
  return true;
}

bool SgIoHelper::getNextMarker(vector<string> &colNames, vector<double> &aVals, 
                   vector<double> &bVals, string &probesetName) {
  if(summaries.is_open()) {
    return getNextMarkerTsvSingle(colNames, aVals, bVals, probesetName);
  } 
  else if(summariesA.is_open() && summariesB.is_open()) {
    return getNextMarkerTwoTsv(colNames, aVals, bVals, probesetName);
  }
  else if(f5Summaries.is_open() && f5SumTsv != NULL) {
    return getNextMarkerFile5(colNames, aVals, bVals, probesetName);
  }
  else {
    Err::errAbort("Couldn't get next marker.");
  }
  return false; // Should never get here, just for compiler
}


SummaryGenotypeEngine::Reg SummaryGenotypeEngine::reg;

SummaryGenotypeEngine * SummaryGenotypeEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == SummaryGenotypeEngine::EngineName())
		return (SummaryGenotypeEngine *)engine;
	return NULL;
}

SummaryGenotypeEngine::SummaryGenotypeEngine() {
    defineOptions();
    defineStates();
}

SummaryGenotypeEngine::~SummaryGenotypeEngine() {
}

void SummaryGenotypeEngine::defineOptions() {
	defineOption("p", "algo-params", PgOpt::STRING_OPT,
		"A string specifying algo parameters. See --explain option for acceptable parameters.",
		"");
	defineOption("", "file5-summaries", PgOpt::STRING_OPT,
		"Interleaved Allele summaries file in file5 format with summaries in /summaries group.",
		"");
	defineOption("s", "summaries", PgOpt::STRING_OPT,
		"Interleaved Allele summaries file.",
		"");
	defineOption("", "summaries-a", PgOpt::STRING_OPT,
		"Allele summarie file for A allele.",
		"");
	defineOption("", "summaries-b", PgOpt::STRING_OPT,
		"Allele summarie file for B allele.",
		"");
	defineOption("", "file5-output", PgOpt::BOOL_OPT,
		"Should output be formatted in file5 (hdf5) output.",
		"false");
	defineOption("", "chrX-snps", PgOpt::STRING_OPT,
		"File containing snps on chrX (non-pseudoautosomal region).",
		"");
	defineOption("", "special-snps",PgOpt::STRING_OPT,
		"File containing all snps of unusual copy (chrX,mito,Y)",
		"");
	defineOption("", "special-sample-snps",PgOpt::STRING_OPT,
		"File containing all sample-specific snps of unusual copy",
		"");
	defineOption("", "write-models", PgOpt::BOOL_OPT,
		"Should we write snp specific models out for analysis? [experimental]",
		"false");
	defineOption("", "read-models-brlmmp", PgOpt::STRING_OPT,
		"File to read precomputed BRLMM-P snp specific models from.",
		"");
	defineOption("", "read-genders", PgOpt::STRING_OPT,
		"Explicitly read genders from a file.",
		"");
	defineOption("", "genotypes", PgOpt::STRING_OPT,
		"File to read seed genotypes from.",
		"");
	defineOption("","select-probes",PgOpt::BOOL_OPT,
		"Should we evaluate model fit for clusters by probe-set/probe[experimental]",
		"false");
	defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
        "Explicitly set the analysis name. This affects output file names (ie prefix) and various meta info.",
		"");
}

void SummaryGenotypeEngine::defineStates() { }


/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void SummaryGenotypeEngine::checkOptionsImp()
{
    defineStates();

	setLibFileOpt("chrX-snps");
	setLibFileOpt("special-snps");
	setLibFileOpt("read-models-brlmmp");

	string outDir = getOpt("out-dir");
	if(!Fs::isWriteableDir(outDir)) {
          if(Fs::mkdirPath(outDir, false) != APT_OK) {
            APT_ERR_ABORT("Can't make or write to directory: " + ToStr(outDir));
          }
	}
}

/**
   This is the "main()" equivalent of the engine. This function has
   two main phases: 1) Creation an initialization of AnalysisStream
   and associated objects. 2) Looping through the list of cdf
   probesets to do the genotyping calls.

   Phase 2 is definitely where most of the time and cpu cycles are
   spent. The list of things that happens for each iteration of phase
   2 is:
   - Read in cel files only keeping data needed for current iteration in RAM.
   - While reading cel files make DM calls.
   - Once cel files are read loop through probesets making BRLMM genotype
   calls.

   Special things that happen in the first iteration of phase 2:
   - Calculate prior if desired.
   - Determine gender using DM calls if desired.
   - Learn chipstream parameters necessary.
   - Attach any reporters to the analysis streams.

   Top of the call stack... This function is getting to be a real
   beast, time for some refactoring...
*/
void SummaryGenotypeEngine::runImp()
{
    // dummy stuff needed to spoof full analysis stream
    ChipLayout layout; 
    std::vector<probeidx_t> analysisOrder;
    std::vector<ChipStream *> iTrans;
    PmOnlyAdjust pmAdjust;
    ProbeSetGroup psGroup;

	string outDir = getOpt("out-dir");
	Util::chompLastIfSep(outDir);
	string summariesFile = getOpt("summaries");
	string summariesFileA = getOpt("summaries-a");
	string summariesFileB = getOpt("summaries-b");
	string quantMethod = "brlmm-p";
	string algoParams = getOpt("algo-params");
	bool selectProbes = getOptBool("select-probes");
	string specialSnpsFile = getOpt("special-snps");
	string specialSampleSnpsFile = getOpt("special-sample-snps");
	string chrXSnps = getOpt("chrX-snps");
	string readModelsBrlmmp=getOpt("read-models-brlmmp");
	bool writeModels = getOptBool("write-models");
	string genderFileIn = getOpt("read-genders");
	string genoTypesFile = getOpt("genotypes");
	string filePrefix = getOpt("set-analysis-name");
	if (filePrefix.empty() == true)
		filePrefix = quantMethod;

    // open up summaries file for input
    vector<string> colNames;
    string summaryFile = getOpt("summaries");
    string summaryAFile = getOpt("summaries-a");
    string summaryBFile = getOpt("summaries-b");
    string summaryA5File = getOpt("file5-summaries");
    SgIoHelper ioHelp(summaryFile, summaryAFile, summaryBFile, summaryA5File, "/");
    colNames = ioHelp.getColNames();

    // Setup dummy sparse mart to spoof the quant method
    SparseMart iMart(analysisOrder, colNames, false);

    // Create Quant Method
    string analysis = quantMethod;
    if(algoParams != "") {
        analysis += "." + algoParams;
    }
    QuantMethodFactory quantMethodFactory(QuantMethodFactory::GenoType);
    Verbose::out(1,"Creating QuantMethod for genotyping from: " + analysis);
    QuantGTypeMethod *qMethod = 
        quantMethodFactory.quantGTypeMethodForString(analysis, layout, QuantMethodFactory::GenoType);
    QuantLabelZ *qBrlmmp = dynamic_cast<QuantLabelZ *>(qMethod);

    // Create Report
    QuantMethodGTypeReport *report = new QuantMethodGTypeReport();
    report->setDirPath(outDir);
    report->setFileprefix(filePrefix);
    report->setFilename("QuantMethodGTypeReport-A5"); // for debugging
    if(getOptBool("file5-output")) {
      report->setCompactFile5Format(20);
    }
    else {
      report->setFormat(affx::TsvReport::FMT_TSV);
      report->setPrecision(5);
    }
    //    report->setFilename("QuantMethodGTypeReport-TXT"); // for debugging
    //    report->setFormat(affx::TsvReport::FMT_A5);

    //

    //
    string qMethodGuid = affxutil::Guid::GenerateNewGuid();
    string execGuid = affxutil::Guid::GenerateNewGuid();
    string timeStr = ""; ///@todo set time string
    string commandLine = ""; ///@todo set command line string
    string execVersion = ""; ///@todo set command line string
    AnalysisInfo info; ///@todo populate analysis info object
    //
    report->addStdHeaders((QuantMethodReport *)report, 
                          execGuid,
                          qMethodGuid,
                          timeStr,
                          commandLine, 
                          execVersion,
                          info);

    report->prepare(*qMethod, iMart);
    
    // setup special SNPs
    SpecialSnpMap specialSnps;
    map<string,bool> haploidSnps;
    if(specialSnpsFile != "") {
        string chip, specialType;
        readSpecialSnpsFile(&specialSnps, &chip, &specialType, specialSnpsFile);
        haploidFromSpecial(specialSnps,haploidSnps);
    } else if(chrXSnps != "") {
        fillInHaploidSnps(chrXSnps, haploidSnps, NULL, true);
        specialFromHaploid(specialSnps, haploidSnps);
    }

    if(!specialSnps.empty())
        qBrlmmp->setSpecialSnps(specialSnps);
    else if(!haploidSnps.empty())
        qBrlmmp->setHaploidSnps(haploidSnps);

    SpecialSampleSnpMap specialSampleSnps;
    SpecialSampleNameVector specialSampleNames;
    if(specialSampleSnpsFile != "") {
        string chip, specialType;
        readSpecialSampleSnpsFile(&specialSampleSnps, &specialSampleNames,
			specialSampleSnpsFile);
    }

    if(!specialSampleSnps.empty())
        qBrlmmp->setSpecialSampleSnps(specialSampleSnps, specialSampleNames);

    // setup genders
    vector<Gender> genders;
    if(genderFileIn != "") {
        FileGenders genderReader(genderFileIn,colNames);
        genders = genderReader.getGenders();
    } else {
        for(int i = 0; i < colNames.size(); i++)
            genders.push_back(UnknownGender);
    }
    qBrlmmp->setGenders(genders);

    // read in SNP models
    if(readModelsBrlmmp!="")
        qBrlmmp->readSnpPriorMap(readModelsBrlmmp.c_str());

    // write SNP models?
    if(writeModels) {
      string outfile = Fs::join(outDir,filePrefix + ".snp-posteriors");
      qBrlmmp->writeSnpPosteriorTsv(outfile,affx::TsvReport::FMT_TSV);
    }


	if (selectProbes)
      {
        string outfile = Fs::join(outDir,filePrefix + ".select-probes.txt");
          qBrlmmp->writeSnpProbeTsv(outfile,affx::TsvReport::FMT_TSV);
      }
	

    // read in seed genotypes
    std::map<const char *,std::vector<affx::GType>, Util::ltstr > knownGenotypes;
    if(!genoTypesFile.empty()) {
      std::vector<std::string> probesets;
      GenoSeedTxt seeds(genoTypesFile.c_str(), layout, colNames, haploidSnps);
      seeds.fillProbesets(probesets);

      for(int i=0; i<probesets.size(); i++) {
        const char *name = Util::cloneString(probesets[i].c_str());
        vector<GType> callVec = seeds.getGenoCalls(name);
        knownGenotypes[name] = callVec;
      }
      qBrlmmp->setKnownGenoTypes(&knownGenotypes);
    }

    // for each SNP -- two rows in summaries file
    Verbose::out(1,"Making calls...");
    vector<double> aVals,bVals;

    string probesetName;
    while(ioHelp.getNextMarker( colNames, aVals, bVals, probesetName)) {

        // set probeset name for report
        qBrlmmp->setProbeSetName(probesetName);

        // fill in allele values
        qBrlmmp->setAlleleValues(aVals,bVals);
    
        // transform the data as specified
        qBrlmmp->transform();
   
        // load up seed genotypes
        qBrlmmp->loadInitialCallsForSnp(probesetName);

        // make the call
        qBrlmmp->computeEstimate();
    
        // report results
        report->report(psGroup, *qMethod, iMart, iTrans, pmAdjust);
    }
    Verbose::out(1,"Done.");

    report->finish(*qMethod);
    // clean up
    delete qBrlmmp;
    delete report;
}
