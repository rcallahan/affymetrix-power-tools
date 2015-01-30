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

//
#include "normalization/normalization.h"
//
#include "util/Convert.h"
#include "util/RowFile.h"
#include "util/Verbose.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
#include <algorithm>
#include <numeric>
#include <vector>
//
#include "normalization/CPPTest/NormalizationTest.h"


CPPUNIT_TEST_SUITE_REGISTRATION( NormalizationTest );

//////////

typedef double InterpT;
typedef double IntensityT;

////////////////////

template <typename T1>
void
util_vector_pushrange(vector<T1>& vec1,T1 s,T1 e,T1 step)
{
  for (T1 i=s;i<e;i=i+step) {
    vec1.push_back(i);
  }
}

template <typename T1>
void
util_vector_pushptrcnt(vector<T1>& vec1,int* ptr,int cnt)
{
  for (T1 i=0;i<cnt;i++) {
    vec1.push_back(*(ptr++));
  }
}

template <typename T1>
void
util_vector_print(vector<T1>& vec1)
{
  printf("VEC==[");
  for (int i=0;i<vec1.size();i++) {
    //printf("%s%6.2f",((i==0)?"":","),((float)vec1[i]));
    printf("%s%2d",((i==0)?"":","),((int)vec1[i]));
  }
  printf("];\n");
}


template <typename T1>
void
util_vecvec_resize(vector<vector<T1> >& vecvec,int size_x,int size_y)
{
  vecvec.resize(size_y);
  for (int y=0;y<vecvec.size();y++) {
    vecvec[y].resize(size_x);
  }
}

template <typename T1>
void
util_vecvec_fill(vector<vector<T1> >& vecvec,T1 val)
{
  for (int y=0;y<vecvec.size();y++) {
    fill(vecvec[y].begin(),vecvec[y].end(),val);
  }
}
  
template <typename T1>
void
util_vecvec_print(vector<vector<T1> >& vecvec)
{
  printf("[\n");
  for (int y=0;y<vecvec.size();y++) {
    printf(" [");
    for (int x=0;x<vecvec[y].size();x++) {
      printf("%6.2f,",(float)vecvec[y][x]);
    }
    printf("],\n");
  }
  printf("]\n");
 }
    

//////////

void 
NormalizationTest::setUp()
{
  // none needed
}

void 
NormalizationTest::tearDown()
{
  // none needed
}

//////////

// TestData ranges
// not zero as we want to test values less than zero too
#define TD_MIN -10 
#define TD_MAX 100

void
NormalizationTest::test_ExtractSketch()
{
  vector<int> data;
  vector<int> sketch;

  // printf("NormalizationTest::test_ExtractSketch()...\n");

  // interleaved
  util_vector_pushrange(data,TD_MIN  ,TD_MAX,2);
  util_vector_pushrange(data,TD_MIN+1,TD_MAX,2);
  // 
  //util_vector_print(data); // dump

  //
  for (int sz=2;sz<21;sz++) {
    // this is required
    sketch.resize(sz);
    // fill with bad data
    fill(sketch.begin(),sketch.end(),-1);
    //
    ExtractSketch<vector<int>::iterator> es(sz);
    es(data.begin(),data.end(),sketch.begin());
    // debugging
    // util_vector_print(sketch);
    // the endpoints should be max and min.
    CPPUNIT_ASSERT(sketch.size()==sz);
    CPPUNIT_ASSERT(sketch[0]==TD_MIN);
    CPPUNIT_ASSERT(sketch[sketch.size()-1]==TD_MAX-1);
  }
}



#define SMIDGE_SIZE 7

void
NormalizationTest::test_InterpolateSketch()
{
  vector<InterpT> dat_sketch;
  vector<InterpT> avg_sketch;
  vector<InterpT> scale_vec;

  for (int scale=1;scale<=10;scale++) {
    // identity vectors
    dat_sketch.clear();
    util_vector_pushrange(dat_sketch,0.0,100.0,1.0);
    avg_sketch.clear();
    util_vector_pushrange(avg_sketch,0.0,100.0,1.0);
    // scale
    scale_vec.resize(avg_sketch.size());
    fill(scale_vec.begin(),scale_vec.end(),scale);
    transform(scale_vec.begin(), scale_vec.end(), dat_sketch.begin(), avg_sketch.begin(),multiplies<InterpT>());
  
    // make our function object
    InterpolateSketch<std::vector<InterpT>::iterator,InterpT> 
      interpolate(dat_sketch.begin(),dat_sketch.end(),avg_sketch.begin(),avg_sketch.end());

    //
    for (int idx=0;idx<dat_sketch.size();idx++) {
      double val_dat,val_out,val_ref,val_err;
      for (int smidge=0;smidge<SMIDGE_SIZE;smidge++) {
        val_dat=dat_sketch[idx]+(smidge/SMIDGE_SIZE);
        val_out=interpolate(val_dat);
        val_ref=val_dat*scale;
        val_err=fabs(val_out-val_ref);
        //printf("inter: %8.2f => %8.2f  (%8.2f %8.6f)\n",val_dat,val_out,val_ref,val_err);
        CPPUNIT_ASSERT(val_err<0.00001);
        //
      }
    }
  }
}




void 
NormalizationTest::test_MedianNormalization () {

  MedianNormalization<vector<vector<IntensityT>::iterator>::iterator> mn;

  vector<vector<IntensityT> > data_vecvec;
  util_vecvec_resize(data_vecvec,10,10);
  util_vecvec_fill(data_vecvec,5.0);
  // debugging
  // util_vecvec_print(data_vecvec);

  vector <vector<IntensityT>::iterator> cellDataItVec;
  for (int i=0;i<data_vecvec.size();i++) {
    cellDataItVec.push_back(data_vecvec[i].begin());
  }

  //
  mn(cellDataItVec.begin(),cellDataItVec.end(), *(cellDataItVec.begin())+(data_vecvec[0].size()),10.0);
  // util_vecvec_print(data_vecvec);
  CPPUNIT_ASSERT(data_vecvec[0][0]==10.0);
  CPPUNIT_ASSERT(data_vecvec[4][9]==10.0);

  //
  //std::iota(data_vecvec[0].begin(),data_vecvec[0].end(),2);
  for (int i=0;i<data_vecvec[0].size();i++) {
    data_vecvec[0][i]=i;
  }
  mn(cellDataItVec.begin(),cellDataItVec.end(), *(cellDataItVec.begin())+(data_vecvec[0].size()),10.0);
  // util_vecvec_print(data_vecvec);
  CPPUNIT_ASSERT(data_vecvec[0][0]== 0.0);
  CPPUNIT_ASSERT(data_vecvec[0][9]==20.0);
  CPPUNIT_ASSERT(data_vecvec[4][9]==10.0);

}



void 
NormalizationTest::test_Normalization () {
  // a list of tests 
  // test_ExtractSketch();
  
}

// Test to see if same result as normalize.quantiles() in bioconductor
void 
NormalizationTest::test_BiocCompat() {
  vector<vector<double> > inputData;
  vector<vector<double> > expectedData;
  double maxDiff = .001;

  QuantileNormalization<vector<vector<double>::iterator >::iterator, double> qn(true);  
  vector<vector <double>::iterator> celDataItVec;

  // read in raw data and result from bioconductor.
  RowFile::matrixFromFile("input/norm-data.txt", inputData);
  RowFile::matrixFromFile("expected/qnorm-bioc-mat.txt", expectedData);
  
  for(int i = 0; i < inputData.size(); i++) {
    celDataItVec.push_back(inputData[i].begin());
  }
  qn(celDataItVec.begin(), celDataItVec.end(), *(celDataItVec.begin()) + inputData[0].size());
  for(int i = 0; i < inputData.size(); i++) {
    for(int j = 0; j < inputData[i].size(); j++) {
      CPPUNIT_ASSERT( fabs(inputData[i][j] - expectedData[i][j]) < maxDiff );
    }
  }
}

bool
NormalizationTest::test_InterpolateVsCoreNorm(vector<vector<float> > inputCoreData,
                                              vector<vector<float> > inputInterpData,
                                              bool doBioc, float maxDiff,
                                              int sketchSize) {
  vector< InterpolateSketch< vector<float>::iterator, double> > m_Interpolators;
  vector< float > m_AverageSketch;
  std::vector< std::vector<float>::iterator > m_Sketches;
  QuantileNormalization<vector<vector<float>::iterator >::iterator, double> qn(doBioc);  
  vector<vector <float>::iterator> celDataItVec;
  // Do in core normalization.
  for(int i = 0; i < inputCoreData.size(); i++) {
    celDataItVec.push_back(inputCoreData[i].begin());
  }
  qn(celDataItVec.begin(), celDataItVec.end(), *(celDataItVec.begin()) + inputCoreData[0].size());

  // Do the sketch normalization with sketch size == data size.
  // make sketch size == data size.
  if(sketchSize == 0)
    sketchSize = inputInterpData[0].size();
    
  ExtractSketch<vector<float>::iterator > extractor(sketchSize);
  // Exract sketches.
  for(int i = 0; i < inputInterpData.size(); i++) {
    std::vector<float> *pVec = new vector<float>(sketchSize);
    m_Sketches.push_back(pVec->begin());
    extractor(inputInterpData[i].begin(), inputInterpData[i].end(), pVec->begin());
  }
  // Average sketches together.
  AverageSketch< vector< vector<float>::iterator >::iterator, vector<float>::iterator, float> as;
  m_AverageSketch.resize(sketchSize);
  as(m_Sketches.begin(), m_Sketches.end(), (*m_Sketches.begin()) + sketchSize, m_AverageSketch.begin());
  // Make interpolator for each chip.
  for(int i = 0; i < m_Sketches.size(); i++) {
    m_Interpolators.push_back(InterpolateSketch<vector<float>::iterator, double>(m_Sketches[i], m_Sketches[i] + sketchSize,
                                                                                 m_AverageSketch.begin(), m_AverageSketch.end(), doBioc));
  }
  // Compare two different normalizations.
  for(int i = 0; i < inputInterpData.size(); i++) {
    for(int j = 0; j < inputInterpData[i].size(); j++) {
      float coreValue = inputCoreData[i][j];
      float interpValue = m_Interpolators[i](inputInterpData[i][j]);
      float value = coreValue - interpValue;
      if(fabs(value) >= maxDiff)
        return false;
    }
  }
  return true;
}

void 
NormalizationTest::fillInWRandData(unsigned int seed, int maxVal, int numChips, int numProbes, 
                     vector<vector<float> > &toFill) {
  int chipIx = 0, probeIx = 0;
  srand(seed);
  toFill.clear();
  for(chipIx = 0; chipIx < numChips; chipIx++) {
    toFill.push_back(vector<float>());
    for(probeIx = 0; probeIx < numProbes; probeIx++) {
      int val = (int) rand() % maxVal;
      toFill[chipIx].push_back(val);
    }
  }
}

void 
NormalizationTest::test_DocExample() {
  vector<vector<float> > data;
  vector<vector <float>::iterator> dataIterVec;
  RowFile::matrixFromFile("input/doc-example.txt", data);
  QuantileNormalization<vector<vector<float>::iterator >::iterator, double> qn(true);  
  for(int i = 0; i < data.size(); i++) {
    dataIterVec.push_back(data[i].begin());
  }
  qn(dataIterVec.begin(), dataIterVec.end(), *(dataIterVec.begin()) + data[0].size());
  for(int rowIx = 0; rowIx < data.size(); rowIx++) {
    for(int colIx = 0; colIx < data[rowIx].size(); colIx++) {
      Verbose::out(2, ToStr(data[rowIx][colIx]) + "\t", false);
    }
    Verbose::out(2, "");
  }
}


void 
NormalizationTest::test_InterpolateVsCore() {
  bool allOk = true;
  vector<vector<float> > inputCoreData;
  vector<vector<float> > inputInterpData;
  // read in raw data and result from bioconductor.
  RowFile::matrixFromFile("input/norm-data.txt", inputCoreData);
  RowFile::matrixFromFile("input/norm-data.txt", inputInterpData);
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, true)) {
    Verbose::out(1,"Interpolate normalization different than core normalization for bioconductor mode.");
    allOk = false;
  }
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, false)) {
    Verbose::out(1,"Interpolate normalization different than core normalization for partial sums (regular) mode.");
    allOk = false;
  }
  inputCoreData.clear();
  inputInterpData.clear();
  fillInWRandData(100, 100, 10, 1000, inputCoreData);
  inputInterpData = inputCoreData;
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, true)) {
    Verbose::out(1,"Interpolate normalization different than core normalization for bioconductor mode with rand data.");
    allOk = false;
  }
  inputCoreData.clear();
  inputInterpData.clear();
  fillInWRandData(100, 100, 10, 1000, inputCoreData);
  inputInterpData = inputCoreData;
  Verbose::out(2, "Doing partial sums normalization with " + ToStr(inputInterpData.size()) + 
               " chips and " + ToStr(inputInterpData[0].size()) + " probes.");
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, false)) {
    Verbose::out(1,"Interpolate normalization different than core normalization for partial sums (regular) mode with rand data.");
    allOk = false;
  }
  inputCoreData.clear();
  inputInterpData.clear();
  fillInWRandData(100, 100, 10, 1000, inputCoreData);
  inputInterpData = inputCoreData;
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, true, 5, 500)) {
    Verbose::out(1,"Interpolate normalization significantly different than core normalization for bioconductor mode with rand data sketchSize = 50%.");
    allOk = false;
  }

  inputCoreData.clear();
  inputInterpData.clear();
  fillInWRandData(100, 100, 10, 1000, inputCoreData);
  inputInterpData = inputCoreData;
  if(!test_InterpolateVsCoreNorm(inputCoreData, inputInterpData, false, .5, 500)) {
    Verbose::out(1,"Interpolate normalization significantly different than core normalization for partial sums (regular) mode with rand data sketchSize = 50%.");
    allOk = false;
  }
  CPPUNIT_ASSERT(allOk);
}
