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
 * @file   ChipStreamTest.cpp
 * @author csugne
 * @date   Fri Nov  4 11:09:48 PST 2005
 * 
 * @brief  Testing the ChipStream functions.
 * 
 */

//
#include "chipstream/ChipStream.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/MedNormTran.h"
#include "chipstream/RmaBgTran.h"
#include "chipstream/SketchQuantNormTran.h"
//
#include "normalization/normalization.h"
#include "util/Convert.h"
#include "util/RowFile.h"
#include "util/TableFile.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//

using namespace std;

/**
 * @class ChipStreamTest
 * @brief cppunit class for testing conversion functions.
 */
class ChipStreamTest : public CppUnit::TestFixture
{
public:
  CPPUNIT_TEST_SUITE( ChipStreamTest );
  CPPUNIT_TEST( testRmaBg );
  CPPUNIT_TEST( testSketchQuantNormTran );
  CPPUNIT_TEST_SUITE_END();


  // blank test.
  void testRmaBg();
  void testSketchQuantNormTran();
  bool testChipStream(ChipStream *stream, const char *fileIn, const char *goldFile);
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( ChipStreamTest );

bool ChipStreamTest::testChipStream(ChipStream *stream, const char *fileIn, const char *goldFile) {
  TableFile toNorm('\t','#',false,false), goldNorm('\t','#',false,false);
  int rowIx, colIx, i;
  vector<vector <float> > data;
  toNorm.open(fileIn);
  goldNorm.open(goldFile);

  std::vector<int> order(toNorm.numCols());
  for (i = 0; i < order.size(); i++) {
    order[i] = i;
  }
  std::vector<std::string> names(toNorm.numRows());
  for (i = 0; i < names.size(); i++) {
    names[i] = ToStr(i);
  }
  DiskIntensityMart* diskMart = new DiskIntensityMart(order, names, names.size() * order.size(), string("."));

  // Read in data.
  for(rowIx = 0; rowIx < toNorm.numRows(); rowIx++) {
    data.push_back( vector<float>() );
    for(colIx = 0; colIx < toNorm.numCols(); colIx++) {
      data[rowIx].push_back(Convert::toFloat(toNorm.getData(rowIx,colIx).c_str()));
    }
    diskMart->setProbeIntensity(rowIx, data[rowIx]);
  }
  
  stream->newDataSet(diskMart);
  stream->endDataSet();
  for(rowIx = 0; rowIx < toNorm.numRows(); rowIx++) {
    for(colIx = 0; colIx < toNorm.numCols(); colIx++) {
      float transformed = stream->getTransformedIntensity(colIx, rowIx);
      float gold = Convert::toFloat(goldNorm.getData(rowIx, colIx).c_str());
      if( !Convert::doubleCloseEnough(gold,transformed,3) )
        return false;
    }
  }
  return true;
} 

/**
   Use R to make some test data sets...
   library(affy)
   x <- c(1:100)
   y <- c(1:100)
   z <- c(1:100)
   y <- 2 *y
   z <- 3 * z
   x <- c(1:100)
   y <- c(1:100)
   z <- c(1:100)
   x[x >=10 & x <= 20] <- 17
   y[y >=15 & y <= 20] <- 19
   z[z >= 20 & z <= 29] <- 23
   y <- y * 2
   z <- z * 3
   y[100] <- 1000
   mat <- rbind(x,y,z)
   nmat <- normalize.quantiles(t(mat))
   write.table(t(nmat), file="qnorm-mat.txt", col.names=F, sep="\t", row.names=F, quote=F)
   y.norm <- normalize.constant(mat[2,], med, FUN=median)
   x.norm <- normalize.constant(mat[1,], med, FUN=median)
   z.norm <- normalize.constant(mat[3,], med, FUN=median)
   mat.medNorm <- rbind(x.norm,y.norm,z.norm)
   write.table(mat.medNorm, file="median-norm.txt", sep="\t", row.names=F, col.names=F, quote=F)
   med <- apply(mat, 1, median)
   med <- median(med)
   
   ave <- apply(mat, 1, mean)
   ave <- mean(ave)
   y.norm <- normalize.constant(mat[2,], ave, FUN=mean)
   x.norm <- normalize.constant(mat[1,], ave, FUN=mean)
   z.norm <- normalize.constant(mat[3,], ave, FUN=mean)
   mat.meanNorm <- rbind(x.norm,y.norm,z.norm)
   write.table(mat.meanNorm, file="mean-norm.txt", sep="\t", row.names=F, col.names=F, quote=F)
*/
void ChipStreamTest::testSketchQuantNormTran() {

  /* Set up chipstreams. */
  SketchQuantNormTran biocNorm(100, true, false, false, 0.0, false); 
  SketchQuantNormTran biocSketchNorm(53,true, false, false, 0.0, false); 
  SketchQuantNormTran affyNorm(100, false, false, false, 0.0, false);
  SketchQuantNormTran affySketchNorm(53, false, false, false, 0.0, false);
  MedNormTran medTran(0, false, true, false);
  MedNormTran medTranTarget(101, false, false, false);
  MedNormTran meanTran(0, true, true, false);
  MedNormTran meanTranTarget(103.65, true, false, false);
  vector<vector<double> > qData;
  RowFile::matrixFromFile("input/norm-data.txt", qData);
  QuantileNormalization<vector<vector<double>::iterator >::iterator, double> qn;
  vector<vector <double>::iterator> cellDataItVec;

  for(int i = 0; i < qData.size(); i++) {
    cellDataItVec.push_back(qData[i].begin());
  }

  qn(cellDataItVec.begin(), cellDataItVec.end(), *(cellDataItVec.begin()) + qData[0].size());
  /* Are the names correct? */
  CPPUNIT_ASSERT( biocNorm.getType() == SKETCHQUANTNORMSTR );
  CPPUNIT_ASSERT( medTran.getType() == MEDNORMSTR );

  /* Match files from incoming data to expected results. */
  CPPUNIT_ASSERT( testChipStream(&biocNorm, "input/norm-data.txt", "expected/qnorm-bioc-mat.txt") );
  CPPUNIT_ASSERT( testChipStream(&biocSketchNorm, "input/norm-data.txt", "expected/bioc-sketch.txt"));
  //  CPPUNIT_ASSERT( testChipStream(&affySketchNorm, "input/norm-data.txt", "expected/affy-sketch.txt"));
  CPPUNIT_ASSERT( testChipStream(&affyNorm, "input/norm-data.txt", "expected/affy.txt"));
  CPPUNIT_ASSERT( testChipStream(&medTran, "input/norm-data.txt", "expected/median-norm.txt"));
  CPPUNIT_ASSERT( testChipStream(&medTranTarget, "input/norm-data.txt", "expected/median-norm.txt"));
  CPPUNIT_ASSERT( testChipStream(&meanTran, "input/norm-data.txt", "expected/mean-norm.txt"));
  CPPUNIT_ASSERT( testChipStream(&meanTranTarget, "input/norm-data.txt", "expected/mean-norm.txt"));
//   TableFile toNorm('\t','#',false,false), biocGoldNorm('\t','#',false,false);
//   ofstream affySketchOut, affyOut, biocSketchOut, medNormOut;
//   vector<vector<float> >data;
//   int rowIx, colIx;
//   Util::mustOpenToWrite(affySketchOut,"expected/affy-sketch.2.txt");
//   float transformed;
//   Util::mustOpenToWrite(affyOut, "expected/affy.2.txt");
//   Util::mustOpenToWrite(biocSketchOut, "expected/bioc-sketch.2.txt");
//   Util::mustOpenToWrite(medNormOut, "expected/med-norm.2.txt");
//   toNorm.open("input/norm-data.txt");
//   biocGoldNorm.open("expected/qnorm-bioc-mat-1000.txt");
//   // Read in data.
//   for(rowIx = 0; rowIx < toNorm.numRows(); rowIx++) {
//     data.push_back( vector<float>() );
//     for(colIx = 0; colIx < toNorm.numCols(); colIx++) {
//       data[rowIx].push_back(Convert::toFloat(toNorm.getData(rowIx,colIx).c_str()));
//     }
//     biocNorm.newChip(data[rowIx]);
//     biocSketchNorm.newChip(data[rowIx]);
//     affyNorm.newChip(data[rowIx]);
  //   affySketchNorm.newChip(data[rowIx]);
   //     medTran.newChip(data[rowIx]);
  //   }
   //   biocNorm.endChips();
   //   biocSketchNorm.endChips();
   //   affyNorm.endChips();
  //   affySketchNorm.endChips();
   //   medTran.endChips();
  //   for(rowIx = 0; rowIx < toNorm.numRows(); rowIx++) {
  //     for(colIx = 0; colIx < toNorm.numCols(); colIx++) {
       //       float transformed = biocNorm.transform(colIx, rowIx, data[rowIx][colIx]);
       //       float gold = Convert::toFloat(biocGoldNorm.getData(rowIx, colIx).c_str());
       //       CPPUNIT_ASSERT( Convert::doubleCloseEnough(gold,transformed) );
       //       transformed = biocSketchNorm.transform(colIx, rowIx, data[rowIx][colIx]);
       //       biocSketchOut << transformed << "\t";
  //       transformed = affySketchNorm.transform(colIx, rowIx, data[rowIx][colIx]);
  //       affySketchOut << transformed << "\t";
       //       transformed = affyNorm.transform(colIx, rowIx, data[rowIx][colIx]);
       //       affyOut << transformed << "\t";
       //       transformed = medTran.transform(colIx, rowIx, data[rowIx][colIx]);
       //       medNormOut << transformed << "\t";
       //     }
       //     biocSketchOut << endl;
  //       affySketchOut << endl;
       //     affyOut << endl;
       //     medNormOut << endl;
  //     }
  //   }
}
   
void ChipStreamTest::testRmaBg() {
  RowFile input, gold;
  vector<float> dataVec, goldVec;
  RmaBgTran rmaBg(512);
  vector<bool> pmProbes;
  vector<string> inputWord, goldWord;
  input.open("input/rma-bg.txt");
  gold.open("expected/rma-bg-gold.txt");
  input.nextRow(inputWord);
  gold.nextRow(goldWord);
  CPPUNIT_ASSERT(goldWord.size() == inputWord.size());
  // Data and gold value from bioconductor "truth"
  for(int i = 0; i < inputWord.size(); i++) {
    dataVec.push_back(Convert::toFloat(inputWord[i].c_str()));
    goldVec.push_back(Convert::toFloat(goldWord[i].c_str()));
    pmProbes.push_back(true);
  }
  rmaBg.setPmProbes(pmProbes);
  rmaBg.initializeData(dataVec);
  rmaBg.endDataSet();
  for(int i = 0; i < dataVec.size(); i++) {
    CPPUNIT_ASSERT( Convert::doubleCloseEnough(goldVec[i], rmaBg.transform(i, 0, dataVec[i])) );
  }
  CPPUNIT_ASSERT( Util::sameString(rmaBg.getType().c_str(),RMABGSTR) );
}
