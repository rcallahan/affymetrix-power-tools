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
#include "chipstream/apt-impute/ImputeEngine.h"
//
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// stolen from impute.cc
#include <stdlib.h>
#include <stdio.h>
#include <iomanip> 
#include <map> 
#include <iterator> 
#include <algorithm> 
#include <numeric> 
#include <functional> 
#include <cmath> 
#include <math.h>
#include "time.h"

// [ToDO]undefining apt min
#ifdef min
	#undef min
#endif
// 
#include "chipstream/apt-impute/Oxford-Impute-version2/data.h"

ImputeEngine::Reg ImputeEngine::reg;

ImputeEngine * ImputeEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == ImputeEngine::EngineName())
		return (ImputeEngine *)engine;
	return NULL;
}

// constructor
ImputeEngine::ImputeEngine() {
  defineOptions();
}

ImputeEngine::~ImputeEngine() {
  clear();
}

void ImputeEngine::clear()
{
}

void ImputeEngine::defineOptions() 
{
  defineOptionSection("Input Options");
  defOptMult("", "haplotype-file", PgOpt::STRING_OPT,
                     "File or know haplotypes, which one row per SNP and one column per haplotype.",
                     "");

  defOptMult("l", "legend-file", PgOpt::STRING_OPT,
                     "Legend file with information about the SNPs in the <-h> (known haplotypes) file.",
                     "");

  defineOption("m", "Recombination-map", PgOpt::STRING_OPT,
                     "Fine-scale recombination map.",
                     "");

  defineOption("g", "Genotype-file", PgOpt::STRING_OPT,
                     "File containing genotypes for a study sample in which we want to impute untyped SNPs.",
                     "");

  defineOption("", "strand_g", PgOpt::STRING_OPT,
					"File for aligning the strand orientation of the <-g> (genotypes) file with the alleles in the reference panel.",
					"");

  defineOption("", "g_ref", PgOpt::STRING_OPT,
                     "File containing unphased genotypes to use as a reference panel for imputation.  Should follow the same format as the <-g> file.",
                     "");

  defineOption("", "strand_g_ref", PgOpt::STRING_OPT,
                     "File specifying strand orientation of genotypes relative to some external reference, such as the Human Genome reference sequence, in the same format as the <-strand_g> file.  This argument is optional, although it is crucial to make sure that the <-g> and <-g-ref> files are aligned with each other (either by making sure they are aligned initially, or aligning each of them to an external reference using the <-strand_g> and <-strand_g_ref> files.",
                     "");

  defineOption("", "output", PgOpt::STRING_OPT,
					"Name of main output file that will contain the imputed genotypes.",
					"");

  defineOption("", "o_gz", PgOpt::BOOL_OPT,
					"Specifies that the output file should be gzipped.",
					"");

  defineOption("", "int", PgOpt::STRING_OPT,
					"Lower and upper boundaries (in base pair position) of region in which to carry out imputation; useful for splitting chromosomes into manageable chunks for analysis.",
                    "");

  defineOption("", "Ne", PgOpt::INT_OPT,
                     "Effective size of the population in which inference is being performed.",
                     "11418");

  defineOption("", "buffer", PgOpt::DOUBLE_OPT,
					"Length of buffer region (in kb) to include on either side of imputation interval to avoid edge effects.",
					"0.0");

  defineOption("", "include_snps", PgOpt::STRING_OPT,
					"A file containing a list of reference-panel-only SNPs to include in the output file; other SNPs of this type will not be imputed.",
					"");

  defineOption("", "iter", PgOpt::INT_OPT,
					"Total number of MCMC iterations to perform, including burn-in.  Default is 30.",
					"-1");

  defineOption("k", "", PgOpt::INT_OPT,
					"Number of 'conditioning states' to use for phasing updates during an MCMC run. Default is 40.",
					"-1");

  defineOption("", "seed", PgOpt::INT_OPT,
					"Random seed.",
					"-1");
}

void ImputeEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ImputeEngine::checkOptionsImp()
{
    defineStates();
	string file;
	vector<string> hap_files;
	vector<string> legned_files;

	// reference panel can be either haplotype or reference_genotype files
	hap_files = getOptVector("haplotype-file");
	legned_files = getOptVector("l");

	if (hap_files.size() == 0 && getOpt("g_ref") == "")
		Err::errAbort("Must specify at least one haplotye file or reference genotype file.");

	if (hap_files.size() != 0)
	{
		if (legned_files.size() == 0)
			Err::errAbort("Haploid file is used.  Must specify at least one lengend file.");

		if (hap_files.size() != legned_files.size())
		{
			Err::errAbort("Number of Haplotype files need to be equal to number of legend files.");
		}
	}

	if (getOpt("g_ref") != "")
	{
		if (getOpt("strand_g_ref") != "")
			Err::errAbort("Reference genotype file is used.  Must specify referenece strand file.");
	}

	if (getOpt("g") == "")
		Err::errAbort("Must specify a genotype file.");

	if (getOpt("m") == "")
		Err::errAbort("Must specify a recombination map file.");

	if (getOpt("int") == "")
		Err::errAbort("Must specify lower and upper boundary of a region.");

	// [ToDo], check if int has lower and upper boundary
	// and make sure the first is lower than the second
	

	if (getOpt("output") == "")
		Err::errAbort("Must specify name of main output file.");

	if (getOpt("strand_g") == "")
		Err::errAbort("Must specify file for aligning the strand orientation.");

	if (getOptInt("Ne") == -1 || getOptInt("Ne") < 0 || getOptInt("Ne") > 30000)
		Err::errAbort("Ne needs to be set in range (0, 30000).");

}

// Compute the various metrics over all chips.
void ImputeEngine::runImp() {
/*
	SnpSummaryReporter reporter;
	int size = getOptInt("buffer-size");

	bool rc = reporter.CreateSnpSummaryReport(getOpt("summary-out-file"), getOptVector("chps"), getOpt("output-format"), size, getOpt("call-file"));	
    if(!rc)
        Err::errAbort("Unable to extract SNP summaries.");
*/

	// copied from impute.cc
	
  long int seed = -1;
  int i, j;
  string in_str;
  bool option_used;
  parameters p;
  data dat;

  // [ToDo]
/*
  // store command line input for printing in impute() function
  for (i = 0; i < argc; ++i) {
    in_str = argv[i];
    p.cmd_call += " " + in_str;
  }
*/

  //////////////////////
  // INPUT PARAMETERS //
  //////////////////////
	vector<string> files;
	files = getOptVector("haplotype-file");
	p.haps_files.clear();
	for (unsigned int i = 0; i < files.size(); i++)
	{
		 p.haps_files.push_back(files[i]);
	}

	files = getOptVector("l");
	p.legend_files.clear();
	for (unsigned int i = 0; i < files.size(); i++)
	{
		 p.legend_files.push_back(files[i]);
	}

    p.gens_file = getOpt("g");
	p.map_file = getOpt("m");

	p.strand_file_g = getOpt("strand_g");

	// using diploid panel
	p.gens_file_ref = getOpt("g_ref"); // file containing reference panel genotypes
	p.strand_file_g_ref = getOpt("strand_g_ref");

	p.out_file = getOpt("output");
	p.o_gz = (getOptBool("o_gz") == true);

	string boundaries = getOpt("int");
	int spacePos = boundaries.find(" ");
	int boundariesLength = boundaries.length();
	string lowerBoundary = boundaries.substr(0, spacePos);

	while (boundaries[spacePos] == ' ' )
	{
		spacePos++;
		if (spacePos > boundariesLength)
			break;
	}

	string upperBoundary = boundaries.substr(spacePos, boundaries.length() - spacePos);

	p.lower = atof(lowerBoundary.c_str());
	p.upper = atof(upperBoundary.c_str());


	p.Ne = getOptDouble("Ne");

	if (getOptDouble("buffer") != 0.0)
		p.buffer = getOptDouble("buffer"); 

	if (getOpt("include_snps") != "")
	{
		p.include_snps_file = getOpt("include_snps");
	}

	if (getOptInt("iter") != -1)
	{
		p.n_mcmc_iter = getOptInt("iter");
	}

	if (getOptInt("k") != -1)
	{
		p.n_cond_dip = getOptInt("k");
	}

  // if random seed has not been supplied manually, choose one at random using the system clock
	if (getOptInt("seed") == -1)
	{
		srand(time(NULL));
		seed = rand();
	}
	else
	{
		seed = getOptInt("seed");
	}

	MT urng((double) seed / RAND_MAX);
	Random::Set(urng);
	p.seed = seed;

	dat.params = p;

	dat.impute();

	cout << endl << "Have a nice day!" << endl;
//	return 0;
}




