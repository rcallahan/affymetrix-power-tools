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

/// @file  apt-matrix-diff.cpp
/// @brief Simple program to check that two matrices are the same.

//
#include "util/AptVersionInfo.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/FsPath.h"
#include "util/MatrixCheck.h"
#include "util/MixedFileCheck.h"
#include "util/PgOptions.h"
#include "util/RowFile.h"
#include "util/TableFile.h"
#include "util/TextFileCheck.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace std;


void define_aptmatrixdiff_options(PgOptions* opts)
{
  opts->setUsage("apt-matrix-diff - check to see if two matrices contain the same values.\n\n"
                 "This can be annoying to do with a regular diff as differences in \n"
                 "precision and outputting functions can lead to representational \n"
                 "differences for two numbers that are the same (i.e. 1.0 == 1).\n\n"
                 "usage:\n"
                 "   apt-matrix-diff -col-skip 1 -line-skip 1 file1.tab file2.tab");
                 
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("m", "match-rows", PgOpt::BOOL_OPT,
                     "Try to match rows based on row names (useful for matrices in different orders).",
                     "false");
  opts->defineOption("a", "allowed-mismatch", PgOpt::INT_OPT,
                     "Maximum number of non-matching numeric values accepted.",
                     "0");
  opts->defineOption("p", "print-mismatch", PgOpt::BOOL_OPT,
                     "Print out entries that don't match.",
                     "false");
  opts->defineOption("", "print-mismatch-max", PgOpt::INT_OPT,
                     "Maximum number of entries to that don't match to report.",
                     "-1");
  opts->defineOption("t", "text", PgOpt::BOOL_OPT,
                     "Compare files as straight text, without conversion to numeric "
                     "data.  This behaves in a fashion very similar to diff, except that "
                     "line-skip number of header lines and all line ending differences "
                     "are ignored.",
                     "false");
  opts->defineOption("", "mixed", PgOpt::BOOL_OPT,
                     "Compare files as a mixture of numeric and non-numeric data. If a "
                     "field can be converted to a floating point number, compare numerically; "
                     "if not, compare as text.",
                     "false");
  opts->defineOption("l", "line-skip", PgOpt::INT_OPT,
                     "Number of lines (header?) to skip.",
                     "0");
  opts->defineOption("c", "col-skip", PgOpt::INT_OPT,
                     "Number of columns (row names?) to skip.",
                     "0");
  opts->defineOption("e", "epsilon", PgOpt::DOUBLE_OPT,
                     "Epsilon tolerance to differences.",
                     ".001");
  opts->defineOption("f", "fraction", PgOpt::DOUBLE_OPT,
                     "Fractional tolerance to differences.",
                     "0.0");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
}

/**
 * Check each entry in two matrices to see if they are the same.  If
 * doing 'match-rows' we will attempt to find the matching row by the
 * row name.
 *
 * @param opts - Options supplied to program.
 * @return number of entries where difference is greater than epsilon
 */
int matrixDiff(PgOptions *opts) {
  int colSkip = opts->getInt("col-skip");
  int rowSkip = opts->getInt("line-skip");
  double epsilon = opts->getDouble("epsilon");
  double fraction = opts->getDouble("fraction");
  int mmAllowed = opts->getInt("allowed-mismatch");
  bool verbose = opts->getBool("print-mismatch");
  int mmMaxShow = opts->getInt("print-mismatch-max");
  const bool text = opts->getBool("text");
  const bool mixed = opts->getBool("mixed");
  bool matchRows = opts->getBool("match-rows");
  if(opts->getArgCount() != 2) {
    Err::errAbort("Must have exactly two file name arguments.");
  }
  if(text && mixed) {
    Err::errAbort("Options 'text' and 'mixed' cannot both be true.");
  }
  Verbose::setLevel(2);
  if (text) {
    TextFileCheck textCheck (opts->getArg(0), opts->getArg(1), rowSkip);
    string errorMsg;
    const bool textOk = textCheck.check (errorMsg);
    if (! textOk && verbose) {
      Verbose::out (1, errorMsg);
    }
    // A return of 0 is good, 1 is bad.
    return ! textOk;
  }
  else if (mixed) {
    MixedFileCheck mixedCheck (opts->getArg(0), opts->getArg(1), epsilon, rowSkip, mmAllowed, fraction);
    if( mmMaxShow >= 0 ) {
      mixedCheck.setMaxError(mmMaxShow);
    }
    string errorMsg;
    const bool mixedOk = mixedCheck.check (errorMsg);
	// errorMsg, if any, is already (unconditionally) reported by check method
	return ! mixedOk;
  }
  else {
    MatrixCheck matrixCheck (opts->getArg(0), opts->getArg(1), epsilon, rowSkip, colSkip, matchRows, mmAllowed, fraction);
	matrixCheck.setPrintMismatch(verbose);
    if( mmMaxShow >= 0 ) {
      matrixCheck.setPrintMismatchMax(mmMaxShow);
    }
    string errorMsg;
    const bool matrixOk = matrixCheck.check (errorMsg);
    if (! matrixOk && verbose) {
      Verbose::out (1, errorMsg);
    }
    return ! matrixOk;
  }
}

/** Everybody's favorite function. */
int main(int argc, const char *argv[]) {
  try {
  const string version = AptVersionInfo::versionToReport();

  int same = 0;

  PgOptions *opts = NULL;
  opts = new PgOptions();
  define_aptmatrixdiff_options(opts);
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
    same = matrixDiff(opts);
  }
  delete opts;
  return same;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
}
