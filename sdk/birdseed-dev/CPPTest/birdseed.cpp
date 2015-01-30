////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/// @file birdseed.cpp

//
#include "birdseed-dev/ClustersReader.h"
#include "birdseed-dev/FitSNPGaussiansPriors3.h"
#include "birdseed-dev/GenotypeCaller.h"
#include "birdseed-dev/IntensitiesParser.h"
#include "birdseed-dev/PriorsReader.h"
//
#include "broadutil/BroadException.h"
#include "broadutil/BroadUtil.h"
#include "broadutil/CelUtil.h"
#include "chipstream/SpecialSnps.h"
#include "util/AptVersionInfo.h"
#include "util/PgOptions.h"
//
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>
//

using namespace std;
using namespace birdseed::dev;

static const char *BIRDSEEDVERSION = "1.2";

static void define_birdseed_options(PgOptions* opts)
{
  opts->setUsage("Usage: birdseed <datafile> <outputcallsfile> <outputconfidencesfile>\n");
  
  /// Input prior files.
  opts->defineOption("", "priors-binary", PgOpt::STRING_OPT,
                     "Read priors from this binary file.  Either this option or --priors-text must be specified.",
                     "");
  opts->defineOption("", "priors-tsv", PgOpt::STRING_OPT,
                     "A Tsv file containing priors (one per SNP).",
                     "");
  opts->defineOption("", "priors-text", PgOpt::STRING_OPT,
                     "A text file containing priors (one per SNP). "
                     "This can also be used when writing a binary priors file to specify the priors to put in the binary file. "
                     "In that case all priors in the text file are included in the binary file.",
                     "");
  
  /// Writing prior files.
  opts->defineOption("", "write-priors-binary", PgOpt::STRING_OPT,
                     "Create a binary priors file. "
                     "If using this option, --priors-directory or --priors-text option must be used to specify input priors.",
                     "");
  opts->defineOption("", "write-priors-tsv", PgOpt::STRING_OPT,
                     "Create a tsv priors file. "
                     "If using this option, one of the --priors- options must be used to specify input priors.",
                     "");
  opts->defineOption("", "write-priors-text", PgOpt::STRING_OPT,
                     "Create a text priors file. "
                     "If using this option, one of the --priors- options must be used to specify input priors.",
                     "");
  
  opts->defineOption("", "clusters", PgOpt::STRING_OPT,
                     "Instead of calculating clusters, read clusters for each SNP from the given file, and use these to call genotypes. "
                     "This option is mutually exclusive with --priors-text and --priors-binary",
                     "");
  
  opts->defineOption("", "chrX-snps", PgOpt::STRING_OPT,
                     "(deprecated) "
                     "File containing snps on chrX (non-pseudoautosomal region). "
                     "This file must have a header line \"all_chrx_no_par\". "
                     "This opt and --special-snps are mutually exclusive. "
                     "'--special-snps' is preferred to this option. "
                     "For a chip with mitochondrial SNPs, --special-snps should be used.",
                     "");
  
  opts->defineOption("", "special-snps", PgOpt::STRING_OPT,
                     "File containing snps of unusual copy number (X, Y, mito). "
                     "This opt and --chrX-snps are mutually exclusive.",
                     "");
  
  opts->defineOption("", "gender-file", PgOpt::STRING_OPT,
                     "File containing a line for each sample in the data file. "
                     "0=female, 1=male, 2=unknown. This file must have a header line 'gender'.",
                     "");
  
  opts->defineOption("", "chip-type", PgOpt::STRING_OPT,
                     "The chip type to be written into a binary priors file, and that is required when reading a binary priors file.",
                     BinaryPriorsReader::kDefaultChipType);
  
  opts->defineOption("c", "correction-factor", PgOpt::DOUBLE_OPT,
                     "Use the supplied value rather than the regular mechanism for calculating the prior correction factor. "
                     "This is useful in debugging.",
                     "-1.0");

  opts->defineOption("", "snp-specific-correction-factor", PgOpt::BOOL_OPT,
                     "Calculate a snp-specific prior correction factor, instead of averaging all the intensities. ",
                     "true");
  
  opts->defineOption("", "average-correction-factor", PgOpt::BOOL_OPT,
                     "Calculate a prior correction factor by averaging all the intensities.",
                     "false");
  
  opts->defineOption("", "write-clusters", PgOpt::STRING_OPT,
                     "Write clusters found to this file.",
                     "");
  
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages "
                     "0 - quiet, 1 - write clusters to stdout, 2 - more messages, 3 - even more messages.",
                     "0");

  opts->defineOption("V", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");

  opts->defineOption("", "disallow_unlikely_clusters", PgOpt::BOOL_OPT,
                     "Don't allow 2-cluster model to match to AA && BB priors, and allow 1-cluster model to match AB prior.",
                     "false");
}

//
BasicPriorsReader* open_priorsSource (PgOptions& opts, std::string *path)
{
  if (opts.get("priors-binary")!="") {
    *path = opts.get("priors-binary");
    return new BinaryPriorsReader(*path);
  }
  else if (opts.get("priors-tsv")!="") {
    *path = opts.get("priors-tsv");
    return new TsvPriorsReader(*path);
  }
  else if (opts.get("priors-text")!="") {
    *path = opts.get("priors-text");
    return new TextPriorsReader(*path);
  }
  else {
    throw BroadException("No priors file. Use: --priors-FMT: ",__FILE__, __LINE__);
  }
}

int main(int argc, const char **argv)
{
  try {
    //
    PgOptions opts;
    define_birdseed_options(&opts);

    PgOpt **algorithmicOpts = getAlgorithmicOpts();
    size_t numOpts = 0;
    for (PgOpt **p = algorithmicOpts; *p != NULL; ++p) {
      ++numOpts;
      opts.addPgOpt(*p);
    }
 
    //
    opts.parseArgv(argv);
    // opts.dump(); // debug
        
    if(opts.getBool("help") == true)
      {
        opts.usage();
        cout << "version:" << endl;
        cout << "   " << BIRDSEEDVERSION << endl;
        return 0;
      }
    else if(opts.getBool("version")) {
      cout << "version: " << BIRDSEEDVERSION << endl;
      return 0;
    }

    if (opts.getBool("disallow_unlikely_clusters")) {
      opts.mustFindOpt("allow_unlikely_clusters")->setValue("false");
    }

    setAlgorithmicParameters(opts);
        
    int verbosity = opts.getInt("verbose");
    setFitSNPGaussiansVerbosity(verbosity);

    double correctionFactor = opts.getDouble("correction-factor");

    //
    // Convert the file to a different format.
    //
    if (opts.get("write-priors-binary")!="" ||
        opts.get("write-priors-tsv")!="" ||
        opts.get("write-priors-text")!="") {
          
      //
      if ((opts.get("priors-binary")=="") &&
          (opts.get("priors-tsv")=="") &&
          (opts.get("priors-text")=="")) {
        fprintf(stderr, "ERROR: When using converting files the '%s', '%s' or '%s' option must also be used.\n",
                "priors-binary",
                "priors-tsv",
                "priors-text");
        opts.usage();
        return 1;
      }

      BasicPriorsReader* priorsSource;
      std::string priorsSource_path;

      priorsSource=open_priorsSource(opts, &priorsSource_path);
            
      cout << "Reading priors from " << priorsSource_path;

      if (opts.get("write-priors-binary")!="") {
        cout << " and writing to binary file " << opts.get("write-priors-binary") << endl;
        priorsSource->writeBinaryPriorsFile(opts.get("write-priors-binary"),opts.get("chip-type"));
      }
      else if (opts.get("write-priors-tsv")!="") {
        cout << " and writing to tsv file " << opts.get("write-priors-tsv") << endl;
        priorsSource->writeTsvPriorsFile(opts.get("write-priors-tsv"),opts.get("chip-type"));
      }
      else if (opts.get("write-priors-text")!="") {
        cout << " and writing to text file " << opts.get("write-priors-text") << endl;
        priorsSource->writeTextPriorsFile(opts.get("write-priors-text"),opts.get("chip-type"));
      }
      else {
        throw BroadException("No priors file. Use: --priors-FMT: ",__FILE__, __LINE__,"");
      }

      return 0;
    }
        
    if (opts.getArgCount() != 3) {
      fprintf(stderr, "ERROR: Incorrect number of arguments\n");
      opts.usage();
      return 1;
    }

    ofstream callsStream(opts.getArg(1).c_str());
    ofstream confidencesStream(opts.getArg(2).c_str());
        
    std::string dataPath = opts.getArg(0);
    bool calculateSNPSpecificCorrectionFactor = opts.getBool("snp-specific-correction-factor") &&
                                                !opts.getBool("average-correction-factor");
    IntensitiesParser intensitiesParser(dataPath, (correctionFactor < 0) && !calculateSNPSpecificCorrectionFactor);
    if (verbosity >= 2) {
      if (calculateSNPSpecificCorrectionFactor) {
        cout << "Using SNP-specific correction factors." << endl;
      } else {
        cout.precision(10);
        cout << "Intensity correction factor: " << intensitiesParser.getCorrectionFactor() << endl;
      }
    }
    if (correctionFactor < 0  && !calculateSNPSpecificCorrectionFactor) {
      correctionFactor = intensitiesParser.getCorrectionFactor();
    }

    if ((opts.get("special-snps")!="") &&
        (opts.get("chrX-file")!="")) {
        fprintf(stderr, "ERROR: Cannot specify both %s and %s options\n", 
                "special-snps", // SPECIAL_SNPS_FILE_OPT.longName, 
                "chrX-snps" // CHRX_FILE_OPT.longName
                );
        opts.usage();
        return 1;
    }

    SpecialSnpMap specialSnps;
    string specialSnpsChipType;
    string specialSnpVersion;
    
    if (opts.get("special-snps")!="") {
      readSpecialSnpsFile(&specialSnps, &specialSnpsChipType, &specialSnpVersion, opts.get("special-snps"));
    } 
    else if (opts.get("chrX-snps")!="") {
      fprintf(stderr, "WARNING: Using deprecated %s option\n", "chrX-snps");
      makeSpecialSnpsFromChrXFile(&specialSnps, opts.get("chrX-snps"));
    }

    vector<Gender> genders;
    if (opts.get("gender-file")!="") {
      vector<int> genderInts;
      readSingleColumnFile(&genderInts, opts.get("gender-file"), "gender");
      for (vector<int>::const_iterator it = genderInts.begin();
           it != genderInts.end(); ++it) {
        genders.push_back(Gender(*it));
      }
    }
        
    callsStream << intensitiesParser.getHeader();
    confidencesStream << intensitiesParser.getHeader();


    BasicPriorsReader* priorsSource;
    std::string priorsSource_name;
    bool usingPriors = true;
    auto_ptr<PriorsReader> priorsReader;
    auto_ptr<ClustersReader> clustersReader;

    if ((opts.get("priors-binary")!="") ||
        (opts.get("priors-tsv")!="") ||
        (opts.get("priors-text")!="")) {
      priorsSource=open_priorsSource(opts, &priorsSource_name);
      priorsReader.reset(new PriorsReader(specialSnps, priorsSource));
    }
    else {
      if (opts.get("clusters")=="") {
        fprintf(stderr, "ERROR: One of --%s, --%s or --%s options must be specified.\n",
                "priors-binary", "priors-text", "clusters");
        opts.usage();
        return 1;
      }
      usingPriors = false;
      clustersReader.reset(new ClustersReader(specialSnps, new TextClustersReader(opts.get("clusters"))));
    }

    auto_ptr<ofstream> clusterOstrm;
    if (opts.get("write-clusters")!="") {
      clusterOstrm.reset(new ofstream(opts.get("write-clusters").c_str()));
    }
    
    while (intensitiesParser.advanceSNP()) {
      auto_ptr<GenotypeCaller> caller;
      if (usingPriors) {
        caller.reset(new GenderAwareGenotypeCaller
                     (intensitiesParser.getCurrentIntensities(),
                      genders,
                      priorsReader.get(),
                      intensitiesParser.getCurrentSNPName().c_str(),
                      correctionFactor,
                      verbosity,
                      clusterOstrm.get()));
      }
      else {
        caller.reset(new GenderAwareForcedClusterGenotypeCaller
                     (intensitiesParser.getCurrentIntensities(),
                      genders,
                      clustersReader.get(),
                      intensitiesParser.getCurrentSNPName().c_str(),
                      correctionFactor,
                      verbosity,
                      clusterOstrm.get()));
      }
            
      callsStream << intensitiesParser.getCurrentSNPName();
      confidencesStream << intensitiesParser.getCurrentSNPName();
      for (size_t i = 0; i < intensitiesParser.getCurrentIntensities().numRows(); ++i) {
        CallAndConfidence call = caller->nextCall();
        // Convert NOCALL to -1
        callsStream << "\t" << (call.call == CallAndConfidence::NOCALL? -1: call.call);
        confidencesStream << "\t" << call.confidence;
      }
      callsStream << endl;
      confidencesStream << endl;
    }
    callsStream.close();
    confidencesStream.close();
    if (clusterOstrm.get() != NULL) {
        clusterOstrm->close();
    }
  } catch (BroadException ex) {
    ex.report(stderr);
    return 1;
  }
  //
  //delete[] combinedOpts;

  return 0;
}

/******************************************************************/
/**************************[END OF birdseed.cpp]*******************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
