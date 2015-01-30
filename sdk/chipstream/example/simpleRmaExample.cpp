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
 * @file   simpleRmaExample.cpp
 * @author Chuck Sugnet
 * @date   Thu Jan 12 11:06:30 2006
 * 
 * @brief A simple example of RMA as implemented in the ChipStream
 * framework to give a flavor of how things work. To keep example
 * simple and readable very little error checking is done as we
 * progress. Please do not use this example for real life data
 * analysis.
 */

#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "chipstream/DiskIntensityMart.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//

/** 
 * Perform Robust Multichip Average (RMA) analysis on the cel files
 * specified. This simple example doesn't do very much in the way of
 * error checking and shouldn't be included as error checking.
 * 
 * @param cdfFile - Layout of probesets on chip.
 * @param outDir - Directory to write results to.
 * @param numCelFiles - How many cel files are in the array.
 * @param celFiles - Array of cel file names to be used.
 */
void doRma(const string cdfFile, const string outDir, vector<string> &celFiles) {
  ChipLayout layout;               // Specifies probesets, locations of features on chip.
  AnalysisStream *aStream = NULL;  // Composite of chipstream, pmadjuster and quantification method to do anlysis.
  AnalysisStreamFactory asFactory; // Factory used to construct an analysis stream according to string representation.
  CelReader cReader;               // Object for getting data from cel files.

  /* Get the layout for the chip. */
  Verbose::out(1,"Opening cdf file: " + ToStr(cdfFile));
  layout.openCdfAll(cdfFile);

  vector<int> desiredOrder;
  for(unsigned int psIx = 0; psIx < layout.getProbeSetCount(); psIx++) {
      ProbeListPacked pl = layout.getProbeListAtIndex(psIx);
      for(int pIx = 0; pIx < pl.probe_cnt(); pIx++) {
          desiredOrder.push_back(pl.get_probeId(pIx));
      }
  }
  DiskIntensityMart iMart(
                desiredOrder,
                celFiles,
                50 * 1048576, 
                outDir,
                "multiChannelExample.tmp", 
                true);

  /* Construct the analysis we are going to do. In this case an
     rma background subtraction followed by a quantile normalization
     using only perfect match probes and summarize probe sets with a 
     median polish. Also known as RMA. Commas separate modules, periods
     deliminate object specific parameters. */
  aStream = asFactory.constructAnalysisStream("rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish",
                                              layout);
  /* Uncomment next line to have "simplerPlierExample instead. */
  //   aStream = asFactory.constructAnalysisStream("quant-norm.pm-mm,plier", layout);

  //@todo need to fix due to refactor for A5
  //aStream->addTextReporter(outDir, 5, true, false, false); 
  /* As cel files are read the chipstream wants to see them and learn parameters. */
  cReader.registerStream(aStream->getChipStreamHead());
  /* As cel files are read the data cache needs to store intensities. */
  cReader.registerIntensityMart(&iMart);

  /* Specify the cel files to open. */
  cReader.setFiles(celFiles);
  Verbose::out(1,"Reading cel files.");
  /* Open, read, and send cel file data one by one to chipstream and imart. */
  cReader.readFiles();

  /* Loop through all of our probe sets and do RMA for each one. */
  Verbose::out(1,"Doing analysis.");
  AnalysisInfo info = aStream->getInfo();
  // @todo: fixme -jhg
  aStream->addStdHeaders("",
                         "",
                         "simpleRmaExample", 
                         // iMart.getRowNames(),
                         "simpleRmaExample", 
                         info);
  for(unsigned int psIx = 0; psIx < layout.getProbeSetCount(); psIx++) {
    ProbeListPacked pList = layout.getProbeListAtIndex(psIx);
    ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
    ProbeSetGroup group(ps);     // psGroup should delete the memory for ps...
    aStream->doAnalysis(group, iMart, true);
  }
  delete aStream;
}

/** Everybody's favorite function. */
int main(int argc, const char *argv[]) {
  if(argc < 4) {
    cerr << 
      "simpleRmaExample - Program to demonstrate AnalysisStream\n"
      "framework. Does RMA on a series of cel files. Does very little error\n"
      "checking to keep example code easy to read. Please do not use this\n"
      "program for real life analysis.\n"
      "\n"
      "usage:\n"
      "   simpleRmaExample cdfFile.cdf outputDir celFile1 celFile2 ... celFileN\n" 
         << endl;
    exit(1);
  }
  string cdfFile = argv[1];
  string outDir = argv[2];
  vector<string> celFiles;
  for(int i=3; i < argc; i++)
      celFiles.push_back(argv[i]);
  doRma(cdfFile, outDir, celFiles);
}
      
    
