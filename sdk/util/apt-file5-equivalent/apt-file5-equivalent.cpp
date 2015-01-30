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

/// @file  apt-file5-equivalent.cpp
/// @brief Simple program to check that two HDF5 files are equivalent.

//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "util/AptVersionInfo.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <list>

using namespace std;


void define_aptfile5Equiv_options(PgOptions* opts)
{
  opts->setUsage("apt-file5-equivalent - check to see if two HDF5 files contain equivalent values.\n\n"
                 "This can be annoying to do with a regular diff as differences in \n"
                 "precision and compiler rounding errors can lead to representational \n"
                 "differences for two numbers that are the same (i.e. 1.0 == 0.99999981).\n\n"
                 "usage:\n"
                 "   apt-file5-equivalent -epsilon 0.001 -c 0.9999 file1.cychp file2.cychp");
                 
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("i", "ignore-columns-file", PgOpt::STRING_OPT,
                     "Path to TSV file with list of columns to ignore.\n" 
					 "(One column with header line IgnoreColumns)",
                     "");
  opts->defineOption("d", "datasets-file", PgOpt::STRING_OPT,
                     "Path to TSV file with list of datasets to compare.\n" 
					 "(Two columns with header line GroupName TsvName)",
                     "");
  opts->defineOption("e", "epsilon", PgOpt::DOUBLE_OPT,
                     "Epsilon tolerance to differences.",
                     ".001");
  opts->defineOption("c", "correlation", PgOpt::DOUBLE_OPT,
                     "Correlation tolerance to differences.",
                     "1.0");
  opts->defineOption("v", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("n", "report-nan-num-diff", PgOpt::BOOL_OPT,
                     "Report differences between NaN and non-NaN values.",
                     "false");
  opts->defineOption("s", "sign-allow-negation", PgOpt::BOOL_OPT,
                     "Allow sign difference between pairs of values compared.\n",
                     "false");
}

/**
 * Check data in all sections of binary HDF5 file to determine eqivalence.
 *
 * @param opts - Options supplied to program.
 * @return status 0 if equivalent, non-zero if not not equivalent.
 */
int file5Equiv(PgOptions *opts) {
  std::string setIgnoreFile = opts->get("ignore-columns-file");
  std::string datasetsFile = opts->get("datasets-file");
  double dEpsilon = opts->getDouble("epsilon");
  double dCorrelationCutoff = opts->getDouble("correlation");
  bool bFlagNaNNumDiff = opts->getBool("report-nan-num-diff");
  bool bAllowNegation = opts->getBool("sign-allow-negation");
  if(opts->getArgCount() != 2) {
    Err::errAbort("Must have exactly two file name arguments.");
  }
  std::string strFileName1 = opts->getArg(0);
  std::string strFileName2 = opts->getArg(1);
  Verbose::setLevel(2);
  // load list of columns to ignore from text file, example for list "Cyto2.Parameters.Parameter"
  std::set<std::string> setIgnore;
  if(setIgnoreFile!="") {
    std::vector<std::string> vectorIgnore;
    affx::TsvFile::extractColToVec(setIgnoreFile,"IgnoreColumns",&vectorIgnore);
    // copy vector to set required by interface
    vector <std::string>::iterator vectorIgnore_Iter = vectorIgnore.begin();
    for( vectorIgnore_Iter=vectorIgnore.begin() ; vectorIgnore_Iter!=vectorIgnore.end() ; vectorIgnore_Iter++ ) {
      setIgnore.insert(*vectorIgnore_Iter);
	}
    Verbose::out(1, "Read " + ToStr(setIgnore.size()) + " ignored parameters from: " + Fs::basename(setIgnoreFile));
  }
  std::vector<std::string> vectorGroupName;
  std::vector<std::string> vectorTsvName;
  // load list GroupName-TsvName for datasets to compare
  if(datasetsFile!="") {
    affx::TsvFile::extractColToVec(datasetsFile,"GroupName",&vectorGroupName);
    affx::TsvFile::extractColToVec(datasetsFile,"TsvName",&vectorTsvName);
    Verbose::out(1, "Read " + ToStr(vectorGroupName.size()) + " GroupName-TsvName value pairs from: " + Fs::basename(datasetsFile));
  }
  // compare HDF5 files
  bool bPassed = true;
  vector <std::string>::iterator vectorGroupName_Iter = vectorGroupName.begin();
  vector <std::string>::iterator vectorTsvName_Iter = vectorTsvName.begin();
  for( vectorGroupName_Iter=vectorGroupName.begin() ; vectorGroupName_Iter!=vectorGroupName.end() ; vectorGroupName_Iter++ ) {
    if (!affx::File5_File::equivalent(strFileName1, strFileName2, *vectorGroupName_Iter, *vectorTsvName_Iter, setIgnore, dEpsilon, dCorrelationCutoff, bAllowNegation, bFlagNaNNumDiff) ) {
      bPassed = false;
      Verbose::out(1, "File5 " + *vectorGroupName_Iter + "." + *vectorTsvName_Iter + " Comparison failed\n" \
		  "for files "+strFileName1 + "\n" + \
		  "      and "+strFileName2);
    }
    vectorTsvName_Iter++;
  }

  // A return of 0 is good, 1 is bad.
  return ! bPassed;
}

/** Everybody's favorite function. */
int main(int argc, const char *argv[]) {
  try {
  const string version = AptVersionInfo::versionToReport();

  int same = 0;

  PgOptions *opts = NULL;
  opts = new PgOptions();
  define_aptfile5Equiv_options(opts);
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
    same = file5Equiv(opts);
  }
  delete opts;
  return same;
  } 
  catch(...) {
    Verbose::out(1,"Unexpected Error: uncaught exception.");
    return 1;
  }
}
