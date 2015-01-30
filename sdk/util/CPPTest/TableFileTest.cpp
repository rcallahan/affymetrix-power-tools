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

//
#include "util/AffxString.h"
#include "util/CPPTest/Setup.h"
#include "util/Fs.h"
#include "util/Util.h"
#include <util/TableFile.h>
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
//
using namespace std;
/**
 * @class TableFileTest
 * @brief cppunit class for testing conversion functions.
 */
class TableFileTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( TableFileTest );
  CPPUNIT_TEST( testOpen );
  CPPUNIT_TEST( testWrite );
  CPPUNIT_TEST( testColIndex );
  CPPUNIT_TEST( testRowIndex );
  //vliber
  CPPUNIT_TEST( testOpenData );
  CPPUNIT_TEST( testWriteData );
  CPPUNIT_TEST(columnFromFileTest);
  CPPUNIT_TEST(columnFromRowFileTest);
  CPPUNIT_TEST(writeVectorTest);
  CPPUNIT_TEST(getHeaderValueTest);
  CPPUNIT_TEST_SUITE_END();

public:
  /** Test reading rows. */
  void testOpen();
  void testWrite();
  void testColIndex();
  void testRowIndex();
  void setUp();
  //vliber
  void testOpenData();
  void testWriteData();
  void columnFromFileTest();
  void columnFromRowFileTest();
  void writeVectorTest();
  void getHeaderValueTest();
  TableFile tf, tf1, tf2;
};

void TableFileTest::setUp() {
  tf.setUseRowNames();
  tf.open(INPUT+"/TableFileTest.txt");
  //if( !Util::directoryWritable(OUTPUT+"/") )
  //		Util::createDir(OUTPUT+"/");
  //}
  if ( !Fs::dirExists(OUTPUT) ) {
    Fs::mkdirPath(OUTPUT);
  }
}

void TableFileTest::testColIndex() {
  CPPUNIT_ASSERT( tf.colIndex((char *)"colName1") == 0);
  CPPUNIT_ASSERT( tf.colIndex((char *)"colName2") == 1);
  CPPUNIT_ASSERT( tf.colIndex((char *)"colName3") == 2);
  //vliber
  CPPUNIT_ASSERT( tf.colIndex((char *)"colName4") == TableFile::npos);
}

void TableFileTest::testRowIndex() {
  CPPUNIT_ASSERT( tf.rowIndex((char *)"row1") == 0);
  CPPUNIT_ASSERT( tf.rowIndex((char *)"row2") == 1);
  CPPUNIT_ASSERT( tf.rowIndex((char *)"blah") == TableFile::npos);
}

void TableFileTest::testWrite()
{
  CPPUNIT_ASSERT(tf.write("output/TableFileTest.txt"));
}
void TableFileTest::testOpen() {
  cout<<endl;
  Verbose::out(1, "***TableFile testcases***");
  Verbose::out(1, "TableFileTest::testOpen");
  CPPUNIT_ASSERT( tf.numRows() == 2 );
  CPPUNIT_ASSERT( tf.numCols() == 3 );
  CPPUNIT_ASSERT( tf.getData(0,0) == string("rowData1") );
  CPPUNIT_ASSERT( tf.getData(0,1) == string("1") );
  CPPUNIT_ASSERT( tf.getData(1,2) == string("2.0") );
  CPPUNIT_ASSERT( tf.getColName(1) == string("colName2") );
  CPPUNIT_ASSERT( tf.getRowName(1) == string("row2") );
  //****vliber****
  //empty input file
  //FATAL ERROR: Nothing after header in file: input/TableFileTest4.txt
  TableFile tf4 ('\t','#',true,true);
  NEGATIVE_TEST(tf4.open(INPUT+"/TableFileTest4.txt"),Except);
}

// test open();getData();getColName();getRowName();numRows();numCols(); getUseRowNames(); setUseRowNames(); vliber 
void TableFileTest::testOpenData() {

  Verbose::out(1, "TableFileTest::testOpenData");
  NEGATIVE_TEST(tf.getData(2,3), std::out_of_range); //todo vliber use index validation before output data
  //string s=tf.getColName(4);
  NEGATIVE_TEST(tf.getColName(4), std::out_of_range);//todo vliber use index validation before output data
  NEGATIVE_TEST(tf.getRowName(3), std::out_of_range);//todo vliber use index validation before output data
  

  vector<std::string> requiredCols1;
  string s1 ("colName1");
  string s2 ("colName2");
  string s3 ("colName3");
  requiredCols1.push_back(s1);
  requiredCols1.push_back(s2);
  requiredCols1.push_back(s3);
  //UseRowNames=true
  tf1.setUseRowNames();
  CPPUNIT_ASSERT(tf1.getUseRowNames()==true);
  tf1.open(INPUT+"/TableFileTest.txt",requiredCols1);
  CPPUNIT_ASSERT( tf1.numRows() == 2 );
  CPPUNIT_ASSERT( tf1.numCols() == 3 );
  CPPUNIT_ASSERT( tf1.getData(0,0) == string("rowData1") );
  CPPUNIT_ASSERT( tf1.getData(0,1) == string("1") );
  CPPUNIT_ASSERT( tf1.getData(1,2) == string("2.0") );
  CPPUNIT_ASSERT( tf1.getColName(1) == string("colName2") );
  CPPUNIT_ASSERT( tf1.getRowName(1) == string("row2") );
  // //UseRowNames=false
  tf2.open(INPUT+"/TableFileTest1.txt",requiredCols1);
  CPPUNIT_ASSERT(tf2.getUseRowNames()==false);
  CPPUNIT_ASSERT( tf2.numRows() == 2 );
  CPPUNIT_ASSERT( tf2.numCols() == 4 );
  CPPUNIT_ASSERT( tf2.getData(0,1) == string("rowData1") );
  CPPUNIT_ASSERT( tf2.getData(0,2) == string("1") );
  CPPUNIT_ASSERT( tf2.getData(1,3) == string("2.0") );
  CPPUNIT_ASSERT( tf2.getColName(1) == string("colName2") );
  //check UseRowNames
  CPPUNIT_ASSERT(tf.getUseRowNames()==true);
  tf.setUseRowNames(false);
  CPPUNIT_ASSERT(tf.getUseRowNames()==false);
  tf.setUseRowNames(true);
  CPPUNIT_ASSERT(tf.getUseRowNames()==true);
  
  //test row and col names
  TableFile tf3 ('\t','#',true,true);
  CPPUNIT_ASSERT(tf.getUseRowNames()==true);
  requiredCols1.clear();
  string s4 ("colName1");
  string s5 ("colName4");
  requiredCols1.push_back(s4);
  requiredCols1.push_back(s5);
  //goal is to receive the messages:
  //FATAL ERROR: Didn't find required column name: 'colName4' in file: input/TableFileTest2.txt
  NEGATIVE_TEST(tf3.open(INPUT+"/TableFileTest2.txt",requiredCols1),Except);
  
  //FATAL ERROR: Expecting 3 words but got 4 at line 5
  TableFile tf4 ('\t','#',true,false);
  NEGATIVE_TEST(tf4.open(INPUT+"/TableFileTest2.txt",requiredCols1),Except);
  
  //FATAL ERROR: Duplicate name: row1 in row names.
  TableFile tf5 ('\t','#',true,true);
  requiredCols1.clear();
  string s6 ("colName1");
  string s7 ("colName3");
  requiredCols1.push_back(s6);
  requiredCols1.push_back(s7);
  NEGATIVE_TEST(tf5.open(INPUT+"/TableFileTest5.txt",requiredCols1),Except);
}


// test that files are identical vliber
void TableFileTest::testWriteData()
{
  Verbose::out(1, "TableFileTest::testWriteData");
  CPPUNIT_ASSERT(tf.write(OUTPUT+"/TableFileTest.txt"));
  AffxString s1=INPUT+"/TableFileTest.txt";
  ifstream ifs1(s1.c_str());
  AffxString s2=OUTPUT+"/TableFileTest.txt";
  ifstream ifs2(s2.c_str());
  string line1;
  string line2;
  while (!ifs2.eof() || !ifs1.eof())
  {
	  getline(ifs2,line2);
	  getline(ifs1,line1);
	  #ifdef WIN32
	  CPPUNIT_ASSERT(line1==line2);
      #endif
  }
  ifs1.close();
  ifs2.close();
}


// test columnFromFile() vliber
void TableFileTest::columnFromFileTest()
{
  Verbose::out(1, "TableFileTest::columnFromFileTest");
  vector<std::string> colValue;
  //verify data
  TableFile::columnFromFile(INPUT+"/TableFileTest2.txt",colValue,"colName2",0,false);
  CPPUNIT_ASSERT(colValue[0]== string("4.1"));
  CPPUNIT_ASSERT(colValue[1]==string("4.2"));
  CPPUNIT_ASSERT(colValue[2]==string("4.2"));
  colValue.clear();
 // duplication in data rows 2 & 3 
 //Warning: column name: colName2 occurs multiple times in: ./input/TableFileTest2.txt using first column.
 //FATAL ERROR: Entry: '4.2' has already been seen in file: input/TableFileTest2.txt in column 1.
  NEGATIVE_TEST(TableFile::columnFromFile(INPUT+"/TableFileTest2.txt",colValue,"colName2",0,true),Except);
  colValue.clear();
  //no data after header
  //FATAL ERROR: Nothing after header in file: input/TableFileTest3.txt -OK
  NEGATIVE_TEST(TableFile::columnFromFile(INPUT+"/TableFileTest3.txt",colValue,"colName2",1,false),Except);
}

// test columnFromRowFile vliber
void TableFileTest::columnFromRowFileTest()
{
  Verbose::out(1, "TableFileTest::columnFromRowFileTest");
  RowFile rf,rf1,rf2;
  rf.open(INPUT+"/TableFileTest2.txt");
  vector<std::string> colValue;
  //coulumn out of range
  //FATAL ERROR: Trying to read column: 5 from row with only: 4 columns at line: 4 in file: input/TableFileTest2.txt
  NEGATIVE_TEST(TableFile::columnFromRowFile(rf,colValue,5,false),Except);
  CPPUNIT_ASSERT(colValue.size()==0);
  rf.close();
  //positive
  rf.open(INPUT+"/TableFileTest2.txt");
  CPPUNIT_ASSERT(TableFile::columnFromRowFile(rf,colValue,2,false));
  CPPUNIT_ASSERT(colValue.size()==4);
  CPPUNIT_ASSERT(colValue[0]== string("colName3"));
  CPPUNIT_ASSERT(colValue[1]== string("1"));
  CPPUNIT_ASSERT(colValue[2]== string("2"));
  CPPUNIT_ASSERT(colValue[3]== string("3"));
  colValue.clear();
  rf.close();
  // ignore duplication data in column
  rf1.open(INPUT+"/TableFileTest2.txt");
  POSITIVE_TEST(TableFile::columnFromRowFile(rf1,colValue,1,false));
  CPPUNIT_ASSERT(colValue.size()==4);
  colValue.clear();
  rf1.close();
  // duplication data in column
  rf2.open(INPUT+"/TableFileTest2.txt");
  //FATAL ERROR: Entry: '4.2' has already been seen in file: input/TableFileTest2.txt in column 1.
  NEGATIVE_TEST(TableFile::columnFromRowFile(rf2,colValue,1,true),Except);
  CPPUNIT_ASSERT(colValue.size()==3);
  rf2.close();
}


// test writeVector vliber
void TableFileTest::writeVectorTest()
{
  Verbose::out(1, "TableFileTest::writeVectorTest");
  //create vector
  RowFile rf;
  rf.open(INPUT+"/TableFileTest2.txt");
  vector<std::string> colValue;
  CPPUNIT_ASSERT(TableFile::columnFromRowFile(rf,colValue,2,false));
  CPPUNIT_ASSERT(colValue.size()==4);
  //write vector to the file
  filebuf fb;
  fb.open ("./output/test2.txt",ios::out);
  ostream out (&fb) ;
  TableFile::writeVector(out,colValue,'%');
  rf.close();
  //test data in file against vector data with delimiter
  ifstream ifs("./output/test2.txt");
  string line1;	  
  getline(ifs,line1);
  CPPUNIT_ASSERT(line1=="colName3%1%2%3");
  ifs.close();
}


// test getHeaderValue vliber
void TableFileTest::getHeaderValueTest()
{
  Verbose::out(1, "TableFileTest::getHeaderValueTest");
  POSITIVE_TEST(tf.getHeaderValue("this"));
  //FATAL ERROR: Don't recognize header key: this1
  NEGATIVE_TEST(tf.getHeaderValue("this1"),Except);
  vector<std::string> colValue;
  colValue=tf.getHeaderValue("this");
  CPPUNIT_ASSERT(colValue.size()==1);
  CPPUNIT_ASSERT(Util::sameString(colValue[0].c_str(),"header"));
}


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( TableFileTest );
