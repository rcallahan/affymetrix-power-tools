////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   apt-snp-compare.cpp
 * @author Chuck Sugnet
 * @date   Fri Feb 10 10:34:35 2006
 * 
 * @brief  Program for comparing the results of a snp prediction method to a known set.
 */

//
#include "file/TsvFile/SnpTable.h"
#include "util/AptVersionInfo.h"
#include "util/PgOptions.h"
#include "util/Util.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//

using namespace std;


class StatCounts {

public:
  int predicted;
  int correct;
  int wrong;

  StatCounts() {
    predicted = 0;
    correct = 0;
    wrong = 0;
  }
};

class GenoTypeCounts {

public:
  StatCounts homoCalls;
  StatCounts hetCalls;
  StatCounts totalCalls;
  int total;
  int noCall;
  int notDetermined;
  GenoTypeCounts() {
    total = 0;
    noCall = 0;
    notDetermined = 0;
  }
};

void define_aptsnpcompare_options(PgOptions* opts)
{
  opts->setUsage("apt-snp-compare compares the results of a snp prediction method to a known results set.\n"
                 "It takes as arguments a matrix file of predicted genotypes (the query) and a matrix\n"
                 "file of predicted genotypes (the target) and computes some summary statistics describing\n"
                 "how well the query matches the the target.\n"
                 "Matrix files should be tab separated, with the first column being snp identifiers\n"
                 "and the first row a header of 'probeset_id' followed by by chip names.\n"
                 "Columns and rows do not have to be matched; the program will determine the genotypes\n"
                 "to compare, using row and column identifiers.\n"
                 "usage:\n"
                 "   apt-snp-compare snp-predictions.mtx snp-known.mtx");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Print help message.",
                     "false");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
}

void scoreSnp(SnpTable &qTable, SnpTable &tTable, GenoTypeCounts &counts, int rowIx, int colIx) {
  static int printWrongCount = 10;
  const string &rowName = qTable.getRowName(rowIx);
  const string &colName = qTable.getColName(colIx);
  int tRowIx = tTable.getRowIndex(rowName);
  int tColIx = tTable.getColIndex(colName);
  if(tRowIx == -1) 
    return;
    //    Err::errAbort("Can't find known SNP for identifier: " + rowName);
  if(tColIx == -1)  
    Err::errAbort("Can't find chip for identifier: " + colName);
  int qType = qTable.getGenotypeForSnp(rowIx, colIx);
  int tType = tTable.getGenotypeForSnp(tRowIx, tColIx);
  StatCounts *counter = NULL;
  counts.total++;
  /* If we don't actually know the genotype of this snp. */
  if(tType == SnpTable::NN) {
    counts.notDetermined++;
    return;
  }
 
  /* Determine if we are counting this as heterozygous or homozygous. 
     Good to break out as the het calls are often harder. */
  if(tType == SnpTable::AB) 
    counter = &counts.hetCalls;
  else if(tType == SnpTable::AA || tType == SnpTable::BB) 
    counter = &counts.homoCalls;
  else 
    Err::errAbort("Don't recognize genotype of: " + ToStr(tType));
  
  /* Score according to correctness. */
  if(qType == SnpTable::NN) {
    counts.noCall++;
    return;
  }
  counter->predicted++;
  counts.totalCalls.predicted++;
  if(qType == tType) {
    counter->correct++;
    counts.totalCalls.correct++;
  }
  else {
    if(printWrongCount > 0) {
      printWrongCount--;
      const string &tRowName = tTable.getRowName(tRowIx);
      const string &tColName = tTable.getColName(tColIx);
      Verbose::out(1, "Mistake example at predicted: " + rowName + " (" + ToStr(rowIx) + ") " + colName + " (" + ToStr(colIx) + 
                   ") vs known: " +
                   tRowName + " (" + ToStr(tRowIx) + ") " + tColName + " (" + ToStr(tColIx) + ") '" + ToStr(qType) + "' vs '" + ToStr(tType) + "'");
    }
    counter->wrong++;
    counts.totalCalls.wrong++;
  }
}

double accuracy(StatCounts &counts) {
  return (double) counts.correct / counts.predicted;
}

void reportResults(GenoTypeCounts &counts) {
  
  double callRate = (double) (counts.total - counts.noCall) / counts.total * 100;
  double callDeterminedRate = (double) (counts.total - counts.noCall - counts.notDetermined) / (counts.total - counts.notDetermined) * 100;
  Verbose::out(1,"Call rate: " + ToStr(callRate) + " For SNPs with known type: " + ToStr(callDeterminedRate));
  Verbose::out(1,"Overall Accuracy: " + ToStr(accuracy(counts.totalCalls) * 100));
  Verbose::out(1,"Het Accuracy: " + ToStr(accuracy(counts.hetCalls) * 100));
  Verbose::out(1,"Hom Accuracy: " + ToStr(accuracy(counts.homoCalls) * 100));
}

    
void compareGenotypes(const std::string& qFile, const std::string& tFile) {
  GenoTypeCounts counts;

  /* Get the data into memory. */
  SnpTable qTable, tTable;
  Verbose::out(1, "Reading file: " + qFile);
  qTable.open(qFile);
  Verbose::out(1, "Reading file: " + tFile);
  tTable.open(tFile);

  Verbose::out(1, "Comparing tables.");
  for(int rowIx = 0; rowIx < qTable.getNumRows(); rowIx++) {
    for(int colIx = 0; colIx < qTable.getNumCols(); colIx++) {
      scoreSnp(qTable, tTable, counts, rowIx, colIx);
    }
  }

  reportResults(counts);
}

int main(int argc, char *argv[]) {
  try {
    const string version = AptVersionInfo::versionToReport();

    PgOptions* opts=new PgOptions;
    define_aptsnpcompare_options(opts);
    opts->parseArgv(argv);
    
    if(opts->getBool("help") || argc == 1) {
        opts->usage();
        cout << "version: " << version << endl;
    }
    else if(opts->getBool("version"))
        cout << "version: " << version << endl;
    else
        compareGenotypes(opts->getArg(0), opts->getArg(1));
    
    delete opts;
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 0;
}
