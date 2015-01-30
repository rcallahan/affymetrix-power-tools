////////////////////////////////////////////////////////////////
//
// Copyright (C) 2013 Affymetrix, Inc.
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

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <file/TsvFile/TsvFile.h>
#include <util/Err.h>
#include <util/Util.h>
#include <util/Fs.h>
#include <util/md5sum.h>
#include <util/Guid.h>
#include <translation/MetabolizerPhenotypingEngine.h>

using namespace affymetrix;
using namespace std;

#define MET_PROGRAM_VERSION "1.0"
#define MET_PROGRAM_NAME "apt-metabolizer"

// Combines two alleles
#define COMBINE_ALLELES(a1, a2) (a1 + "/" + a2)

// Combine activities
#define COMBINE_ACTIVITIES(a1, a2) (a1 + "/" + a2)

// Unknown phenotype call
#define UNKNOWN_PHENOTYPE "unknown"

MetabolizerPhenotypingEngine::Reg MetabolizerPhenotypingEngine::reg;

MetabolizerPhenotypingEngine * MetabolizerPhenotypingEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == MetabolizerPhenotypingEngine::EngineName())
		return (MetabolizerPhenotypingEngine *)engine;
	return NULL;
}

MetabolizerPhenotypingEngine::MetabolizerPhenotypingEngine()
{
	defineOptions();
}

MetabolizerPhenotypingEngine::~MetabolizerPhenotypingEngine()
{
}

void MetabolizerPhenotypingEngine::Init()
{
	metabolizerResults.clear();
	metabolizerPhenotypeMap.clear();
	phenotypeCallDescs.clear();
	userInformation.clear();
}

void MetabolizerPhenotypingEngine::ReadMetabolizerBinFile()
{
	std::string gene;
	std::string allele1;
	std::string allele2;
	std::string phenotype;
	std::string activity1;
	std::string activity2;
	std::string allelePair;
	PhenotypeCallType pcall;

	// Open the file.
	affx::TsvFile tsv;
	std::string metabolizerBinFile = getOpt("metabolizer-file");
	if (tsv.open(metabolizerBinFile) != affx::TSV_OK)
	{
		std::string err = "Unable to open the metabolizer bin file.";
		Err::errAbort(err);
	}

	// Get the phenotype description and information
	phenotypeCallDescs.clear();
	tsv.getHeaderMatchingKeySubstr("PhenotypeCallDesc", phenotypeCallDescs);
	userInformation.clear();
	tsv.getHeaderMatchingKeySubstr("Info", userInformation);

	// Read the body of the file for the mapping data.
	tsv.bind(0, "gene", &gene, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "allele_1", &allele1, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "allele_2", &allele2, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "phenotype", &phenotype, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "activity_1", &activity1, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "activity_2", &activity2, affx::TSV_BIND_REQUIRED);
	while (tsv.nextLine() == affx::TSV_OK)
	{
		if (metabolizerPhenotypeMap.find(gene) == metabolizerPhenotypeMap.end())
			metabolizerPhenotypeMap[gene].clear();

		allelePair = COMBINE_ALLELES(allele1, allele2);
		pcall = PhenotypeCallType(phenotype, activity1, activity2);
		if (metabolizerPhenotypeMap[gene].find(allelePair) != metabolizerPhenotypeMap[gene].end())
		{
			if (pcall.phenotype != metabolizerPhenotypeMap[gene][allelePair].phenotype)
				Err::errAbort("Duplicated allele pair with conflicting information: Gene:" + gene + ", Alleles:" + allele1 + "/" + allele2);
		}
		metabolizerPhenotypeMap[gene][allelePair] = pcall;

		allelePair = COMBINE_ALLELES(allele2, allele1);
		pcall = PhenotypeCallType(phenotype, activity2, activity1);
		if (metabolizerPhenotypeMap[gene].find(allelePair) != metabolizerPhenotypeMap[gene].end())
		{
			if (pcall.phenotype != metabolizerPhenotypeMap[gene][allelePair].phenotype)
				Err::errAbort("Duplicated allele pair with conflicting information: Gene:" + gene + ", Alleles:" + allele2 + "/" + allele1);
		}
		metabolizerPhenotypeMap[gene][allelePair] = pcall;
	}
	tsv.close();
}

static int FindStartingAttributeColumnIndex(affx::TsvFile &tsv)
{
	int n = tsv.getColumnCount(0);
	int idx = n;
	for (int i=0; i<n; i++)
	{
		std::string name = tsv.getColumnName(0, i);
		if (name == "Override Comment")
		{
			idx = i+1;
			break;
		}
	}
	return idx;
}

static void WriteHeader(MetabolizerPhenotypingEngine *eng, affx::TsvFile &tsv, int idx, std::ofstream &str)
{
	// Output the header
	str << "# For research use only. Not for diagnostic purposes." << std::endl;
	str << "#%report-guid=" << affxutil::Guid::GenerateNewGuid() << std::endl;
	str << "#%Program=" << MET_PROGRAM_NAME << std::endl;
	str << "#%Version=" << MET_PROGRAM_VERSION << std::endl;
	str << "#%Date=" << Util::getTimeStamp() << std::endl;
	str << "#%MetabolizerFile=" << eng->getOpt("metabolizer-file") << std::endl;
	str << "#%AlleleFile=" << eng->getOpt("allele-file") << std::endl;
	for (int i=0; i<(int)eng->PhenotypeCallDescs().size(); i++)
		str << "#%PhenotypeCallDesc=" << eng->PhenotypeCallDescs()[i] << std::endl;
	for (int i=0; i<(int)eng->UserInformation().size(); i++)
		str << "#%Info=" << eng->UserInformation()[i] << std::endl;
	str << "Index" << "\t"
		<< "CHP File" << "\t"
		<< "Gene" << "\t"
		<< "Phenotype Call" << "\t"
		<< "Gene Activity" << "\t"
		<< "Known Call" << "\t"
		<< "Unknown Call" << "\t"
		<< "Interpretation Code";
	int n = tsv.getColumnCount(0);
	for (int i=idx; i<n; i++)
	{
		std::string name = tsv.getColumnName(0, i);
		str << "\t" << name;
	}
	str << std::endl;
}

static void OpenAlleleFile(MetabolizerPhenotypingEngine *eng, affx::TsvFile &tsv, MetabolizerResultType &result, int &sampleIdx, int &ncols, std::vector<std::string> &sampleCols)
{
	std::string translationTableFile = eng->getOpt("allele-file");
	tsv.m_optQuoteChar1 = 0;
	tsv.m_optQuoteChar2 = 0;
	if (tsv.open(translationTableFile) != affx::TSV_OK)
	{
		std::string err = "Unable to open the metabolizer bin file.";
		Err::errAbort(err);
	}
	tsv.bind(0, "Index", &result.index, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "CHP Filename", &result.chpFile, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "Gene", &result.gene, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "Known Call", &result.knownCall, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "Unknown Call", &result.unknownCall, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "Interpretation Code", &result.code, affx::TSV_BIND_REQUIRED);
	sampleIdx = FindStartingAttributeColumnIndex(tsv);
	ncols = tsv.getColumnCount(0);
	sampleCols.resize(ncols-sampleIdx);
	for (int i=sampleIdx; i<ncols; i++)
		tsv.bind(0, i, &sampleCols[i-sampleIdx], affx::TSV_BIND_REQUIRED);
}

static void OpenResultFile(MetabolizerPhenotypingEngine *eng, std::ofstream &str)
{
	std::string outputFile = eng->getOpt("output-file");
	if (outputFile.empty() == false)
	{
		str.open(outputFile.c_str(), std::ios::out);
		if (!str)
		{
			std::string err = "Unable to write to the output file: " + outputFile;
			Err::errAbort(err);
		}
	}
}

static void ModifyIndex(std::string &index)
{
	size_t idx = index.rfind('-');
	if (idx != std::string::npos)
		index = index.substr(0, idx);
}

static void WriteResult(std::ofstream &str, MetabolizerResultType &result, int sampleIdx, int ncols, std::vector<std::string> &sampleCols)
{
	str << result.index << "\t"
		<< result.chpFile << "\t"
		<< result.gene << "\t"
		<< result.finalPhenotype << "\t"
		<< result.finalActivity << "\t"
		<< result.knownCall << "\t"
		<< result.unknownCall << "\t"
		<< result.code;
	for (int i=sampleIdx; i<ncols; i++)
		str << "\t" << sampleCols[i-sampleIdx];
	str << std::endl;
}

void MetabolizerPhenotypingEngine::ProcessAlleleTranslationFile()
{
	// Open the file and bind to the required columns.
	affx::TsvFile tsv;
	MetabolizerResultType result;
	int sampleIdx;
	int ncols;
	std::vector<std::string> sampleCols;
	OpenAlleleFile(this, tsv, result, sampleIdx, ncols, sampleCols);

	// Open the results file for writing
	std::ofstream str;
	OpenResultFile(this, str);
	if (str)
		WriteHeader(this, tsv, sampleIdx, str);

	// Process each line of the translation file.
	std::map<std::string, bool> geneFound;
	std::string chp;
	while (tsv.nextLine() == affx::TSV_OK)
	{
		result.finalPhenotype.clear();
		result.phenotypes.clear();
		result.finalActivity.clear();
		result.activities.clear();

		// Skip any genes not found in the translation map.
		if (metabolizerPhenotypeMap.find(result.gene) == metabolizerPhenotypeMap.end())
			continue;

		// Check if the gene of the current CHP file was found. If already processed then skip the line since
		// the allele calls are the same for the repeated lines.
		if (chp != result.chpFile)
			geneFound.clear();
		chp = result.chpFile;
		if (geneFound.find(result.gene) != geneFound.end())
			continue;
		geneFound[result.gene] = true;

		// Modify the index (strip off allele call)
		ModifyIndex(result.index);

		// Parse each diplotype and process it
		std::vector<std::string> calls;
		Util::chopString(result.knownCall, ',', calls);
		for (std::vector<std::string>::iterator it=calls.begin(); it!=calls.end(); it++)
		{
			std::string &call = *it;
			std::string activities = "";
			std::string phenotype = UNKNOWN_PHENOTYPE;
			if (metabolizerPhenotypeMap[result.gene].find(call) != metabolizerPhenotypeMap[result.gene].end())
			{
				phenotype = metabolizerPhenotypeMap[result.gene][call].phenotype;
				activities = COMBINE_ACTIVITIES(metabolizerPhenotypeMap[result.gene][call].activity1, metabolizerPhenotypeMap[result.gene][call].activity2);
			}
			result.phenotypes.push_back(phenotype);
			result.activities.push_back(activities);
		}

		// Compute the final call.
		if (result.phenotypes.size() == 0)
		{
			result.finalPhenotype = UNKNOWN_PHENOTYPE;
			result.finalActivity = "";
		}
		else if (result.phenotypes.size() == 1)
		{
			result.finalPhenotype = result.phenotypes[0];
			result.finalActivity = result.activities[0];
		}
		else
		{
			bool allSame = true;
			std::string first;
			std::string all;
			std::string firstAct;
			std::string allAct;
			for (int i=0; i<(int)result.phenotypes.size(); i++)
			{
				std::string p = result.phenotypes[i];
				std::string a = result.activities[i];
				if (first.empty() == true)
				{
					first = p;
					firstAct = a;
				}
				else if (first != p)
				{
					allSame = false;
				}

				if (all.empty() == true)
				{
					all = p;
					allAct = a;
				}
				else
				{
					all += "," + p;
					allAct += "," + a;
				}
			}
			result.finalPhenotype = (allSame == true ? first : all);
			result.finalActivity = (allSame == true ? firstAct : allAct);
		}
		metabolizerResults.push_back(result);

		// Output the results
		if (str)
			WriteResult(str, result, sampleIdx, ncols, sampleCols);
	}
	if (str)
		str.close();
	tsv.close();
}

void MetabolizerPhenotypingEngine::runImp()
{
	Init();
	ReadMetabolizerBinFile();
	ProcessAlleleTranslationFile();
	CreateMd5File();
}

void MetabolizerPhenotypingEngine::CreateMd5File()
{
	std::string alleleFile = getOpt("allele-file");
	std::string outputMetabolizerFile = getOpt("output-file");
	std::string md5file = getOpt("md5-file");
	std::string salt = getOpt("salt");
	std::string ttSum;
	std::string mbSum;
	affx::md5sum md5;

	// Check if the file has been set.
	if (md5file.empty() == true)
		return;

	// Store the MD5's in a file.
	std::ofstream md5out(md5file.c_str(), ios::out);
	md5out << "File\tMD5" << std::endl;
	if (alleleFile.empty() == false)
	{
		md5.ofFile(alleleFile, ttSum, salt);
		md5out << Fs::basename(alleleFile) << "\t" << ttSum << std::endl;
	}
	if (outputMetabolizerFile.empty() == false)
	{
		md5.ofFile(outputMetabolizerFile, mbSum, salt);
		md5out << Fs::basename(outputMetabolizerFile) << "\t" << mbSum << std::endl;
	}
	md5out.close();
}

void MetabolizerPhenotypingEngine::defineOptions()
{
	defineOptionSection("Inputs");
	defineOption("", "allele-file", PgOpt::STRING_OPT,
					"The file containing the allele calls", "");
	defineOption("", "metabolizer-file", PgOpt::STRING_OPT,
					 "The metabolizer bin file for translating allele calls to metabolizer phenotypes", "");
	defineOption("", "salt", PgOpt::STRING_OPT,
					 "The salt value to add to the MD5s.", "salt");
	defineOptionSection("Outputs");
	defineOption("", "output-file", PgOpt::STRING_OPT,
					 "The output file for the metabolizer bin phenotypes", "");
	defineOption("", "md5-file", PgOpt::STRING_OPT,
					 "The output file to store the MD5 values.", "");
}

void MetabolizerPhenotypingEngine::checkOptionsImp()
{
	setLibFileOpt("metabolizer-file");

	std::string metabolizerBinFile = getOpt("metabolizer-file");
	std::string translationTableFile = getOpt("allele-file");

	if (metabolizerBinFile.empty() == true)
		Err::errAbort("You must specify the metabolizer bin matrix file.");

	if (Fs::exists(metabolizerBinFile) == false)
		Err::errAbort("The specified metabolizer bin matrix file does not exist.");

	if (translationTableFile.empty() == true)
		Err::errAbort("You must specify the translation table file.");

	if (Fs::exists(translationTableFile) == false)
		Err::errAbort("The specified translation table file does not exist.");
}
