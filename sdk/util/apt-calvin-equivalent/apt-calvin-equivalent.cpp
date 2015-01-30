////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

/// @file  apt-calvin-equivalent.cpp
/// @brief Simple program to check that two Calvin files are equivalent.

//
#include "calvin_files/utils/src/Calvin.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AptVersionInfo.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace std;


void define_aptcalvinEquiv_options(PgOptions* opts)
{
  opts->setUsage("apt-calvin-equivalent - check to see if two Calvin files contain equivalent values.\n\n"
                 "This can be annoying to do with a regular diff as differences in \n"
                 "precision and compiler rounding errors can lead to representational \n"
                 "differences for two numbers that are the same (i.e. 1.0 == 0.99999981).\n\n"
                 "usage:\n"
                 "   apt-calvin-equivalent -epsilon 0.001 -fraction 0.0001 -c 0.9999 file1.cychp file2.cychp");
                 
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("k", "check-header", PgOpt::BOOL_OPT,
                     "Check header parameters and values match.",
                     "true");
  opts->defineOption("i", "ignore-params-file", PgOpt::STRING_OPT,
                     "Path to TSV file with list of header parameters, datasets, and \n" 
					 "dataset columns to ignore. (One column with header line IgnoreParameters)",
                     "");
  opts->defineOption("m", "epsilon-map-file", PgOpt::STRING_OPT,
                     "Path to TSV file with mapEpsilon parameter threshold values.\n" 
					 "(Two columns with header line MapParameter MapEpsilon)",
                     "");
  opts->defineOption("l", "message-limit", PgOpt::INT_OPT,
                     "Maximum number messages to report.",
                     "1000");
  opts->defineOption("e", "epsilon", PgOpt::DOUBLE_OPT,
                     "Epsilon tolerance to differences.",
                     ".001");
  opts->defineOption("f", "fraction", PgOpt::DOUBLE_OPT,
                     "Fractional tolerance to differences.",
                     "0.0");
  opts->defineOption("c", "correlation", PgOpt::DOUBLE_OPT,
                     "Correlation tolerance to differences.",
                     "1.0");
  opts->defineOption("v", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("n", "report-nan-num-diff", PgOpt::BOOL_OPT,
                     "Report differences between NaN and non-NaN values.",
                     "false");
}

/**
 * Check data in all sections of binary Calvin file to determine eqivalence.
 *
 * @param opts - Options supplied to program.
 * @return status 0 if equivalent, non-zero if not not equivalent.
 */
int calvinEquiv(PgOptions *opts) {
  float fEpsilon = (float) opts->getDouble("epsilon");
  float fFraction = (float) opts->getDouble("fraction");
  double dCorrelationCutoff = opts->getDouble("correlation");
  bool bCheckHeader = opts->getBool("check-header");
  int iMessageLimit = opts->getInt("message-limit");
  bool bFlagNaNNumDiff = opts->getBool("report-nan-num-diff");
  std::string setIgnoreFile = opts->get("ignore-params-file");
  std::string mapEpsilonFile = opts->get("epsilon-map-file");
  if(opts->getArgCount() != 2) {
    Err::errAbort("Must have exactly two file name arguments.");
  }
  Verbose::setLevel(2);
  // load list if parameters to ignore from text file
  std::set<std::string> setIgnore;
  if(setIgnoreFile!="") {
    std::vector<std::string> vectorIgnore;
    affx::TsvFile::extractColToVec(setIgnoreFile,"IgnoreParameters",&vectorIgnore);
    // copy vector to set required by interface
    vector <std::string>::iterator vectorIgnore_Iter = vectorIgnore.begin();
    for( vectorIgnore_Iter=vectorIgnore.begin() ; vectorIgnore_Iter!=vectorIgnore.end() ; vectorIgnore_Iter++ ) {
      setIgnore.insert(*vectorIgnore_Iter);
	}
    Verbose::out(1, "Read " + ToStr(setIgnore.size()) + " ignored parameters from: " + Fs::basename(setIgnoreFile));
  }
  std::map<std::string, float> mapEpsilon;
  // load mapping of parameter names to epsilon values from text file
  if(mapEpsilonFile!="") {
    std::vector<std::string> vectorMapParameter;
    std::vector<std::string> vectorMapEpsilon;
    affx::TsvFile::extractColToVec(mapEpsilonFile,"MapParameter",&vectorMapParameter);
    affx::TsvFile::extractColToVec(mapEpsilonFile,"MapEpsilon",&vectorMapEpsilon);
    // copy vector to set required by interface
	vector <std::string>::iterator vectorMapEpsilon_Iter = vectorMapEpsilon.begin();
    for( vector <std::string>::iterator vectorMapParameter_Iter = vectorMapParameter.begin() ; 
		 vectorMapParameter_Iter != vectorMapParameter.end() ; 
		 vectorMapParameter_Iter++ ) {
      bool success;
      double dMapEpsilon = Convert::toDoubleCheck(*vectorMapEpsilon_Iter, &success);
	  if(success) {
        mapEpsilon[*vectorMapParameter_Iter] = (float)dMapEpsilon;
	  }
	  else {
        Err::errAbort("Could not convert MapEpsilon value '" + *vectorMapEpsilon_Iter + "' to numeric value for parameter '" + *vectorMapParameter_Iter + "'.");	    
	  }
	  vectorMapEpsilon_Iter++;
	}
    Verbose::out(1, "Read " + ToStr(mapEpsilon.size()) + " MapEpsilon parameter value pairs from: " + Fs::basename(mapEpsilonFile));
  }
  // compare Calvin files
  const bool ok = Calvin::equivalent(opts->getArg(0), opts->getArg(1), setIgnore, setIgnore, mapEpsilon, fEpsilon, dCorrelationCutoff,bCheckHeader, iMessageLimit, fFraction, bFlagNaNNumDiff);
  // A return of 0 is good, 1 is bad.
  return ! ok;
}

/** Everybody's favorite function. */
int main(int argc, const char *argv[]) {
  try {
  const string version = AptVersionInfo::versionToReport();

  int same = 0;

  PgOptions *opts = NULL;
  opts = new PgOptions();
  define_aptcalvinEquiv_options(opts);
  opts->parseArgv(argv);

  // Do we need help?
  if(opts->getBool("help") || argc == 1) {
    opts->usage();
    cout << "version: " << version << endl;
  }
  else if(opts->getBool("version")) {
    cout << "version: " << version << endl;
  }
  else {
    same = calvinEquiv(opts);
  }
  delete opts;
  return same;
  } 
  catch(...) {
    Verbose::out(1,"Unexpected Error: uncaught exception.");
    return 1;
  }
}
