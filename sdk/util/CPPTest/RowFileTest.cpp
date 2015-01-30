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

#ifndef ROWFILETEST_H
#define ROWFILETEST_H

#include "util/AffxString.h"
#include "util/RowFile.h"
#include "util/Util.h"
//
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
//
#include "util/CPPTest/Setup.h"

#include "util/Fs.h"
using namespace std;
//
/**
 * @class RowFileTest
 * @brief cppunit class for testing conversion functions.
 * last change by vliber on 05/19/09
 */
class RowFileTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( RowFileTest );
  CPPUNIT_TEST( testReadHeader );
  CPPUNIT_TEST( testGetSet );
  CPPUNIT_TEST( testLineEnding );
  CPPUNIT_TEST( testReaders );
  CPPUNIT_TEST( testNextRow );
  //vliber
  CPPUNIT_TEST(nextLineTest);
  CPPUNIT_TEST(nextLineReuseTest);
  CPPUNIT_TEST(nextRealLineTest);
  CPPUNIT_TEST(nextCStringRowExpectTest);
  CPPUNIT_TEST(nextCStringRowTest);
  CPPUNIT_TEST(writeHeaderTest);
  CPPUNIT_TEST(openTest);
  CPPUNIT_TEST(matrixFromFileTest);
  CPPUNIT_TEST_SUITE_END();

public:
  void testReadHeader();
  // Test our getters and setters
  void testGetSet();
  // Test checking the type of files
  void testLineEnding();
  // Test some of the readers.
  void testReaders();
  // Test reading the next row.
  void testNextRow();
  //vliber
  void nextLineTest();
  void nextLineReuseTest();
  void nextRealLineTest();
  void nextCStringRowExpectTest();
  void nextCStringRowTest();
  void writeHeaderTest();
  void openTest();
  void matrixFromFileTest();

  void setUp(){
    if ( !Fs::dirExists(OUTPUT) ) {
      Fs::mkdirPath(OUTPUT, false);
    }
  };
};

#endif

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RowFileTest );

void RowFileTest::testReadHeader()
{
  cout<<endl;
  Verbose::out(1, "***RowFile testcases***");
  Verbose::out(1, "RowFileTest::testReadHeader");
  RowFile rf;
  vector<string> lines;
  map<string,vector<string> > header;
  rf.open(INPUT+"/RowFileTestHeader.txt");
  RowFile::readHeader(rf, header, lines);
  CPPUNIT_ASSERT( header["this"][0] == "that" );
  CPPUNIT_ASSERT( header["this"][1] == "that2" );
  CPPUNIT_ASSERT( header["testing"][0] == "good" );
  CPPUNIT_ASSERT( header["love"][0] == "life" );
  CPPUNIT_ASSERT( header.find("Hello") == header.end() );
  //vliber
  CPPUNIT_ASSERT( header["this"][2] == "" );
}
void RowFileTest::testGetSet()
{
  Verbose::out(1, "RowFileTest::testGetSet");
  RowFile rf;
  rf.setDelim(',');
  CPPUNIT_ASSERT( ',' == rf.getDelim() );
  rf.setComment('*');
  CPPUNIT_ASSERT( '*' == rf.getComment() );
}

void RowFileTest::testLineEnding()
{
  Verbose::out(1, "RowFileTest::testLineEnding");
  RowFile rf, dos, mac, unknown;
  rf.open(INPUT+"/RowFileTest.1.txt");
  dos.open(INPUT+"/msdosFormat.txt");
  mac.open(INPUT+"/macFormat.txt");
  unknown.open(INPUT+"/noLineEnding.txt");
  CPPUNIT_ASSERT( rf.getFileType() == RowFile::UNIX);
  //Under CYGWIN installed to use DOS line endings, this will look like UNIX
  //CPPUNIT_ASSERT( dos.getFileType() == RowFile::DOS);
  CPPUNIT_ASSERT( mac.getFileType() == RowFile::MAC);
  CPPUNIT_ASSERT( unknown.getFileType() == RowFile::UNKNOWN);
}

void RowFileTest::testReaders()
{
  Verbose::out(1, "RowFileTest::testReaders");
  RowFile rf;
  const string *s = NULL;
  vector<const char *> cstrings;
  vector<string> words;
  rf.open(INPUT+"/RowFileTest.1.txt");
  rf.nextCStringRow(cstrings);
  CPPUNIT_ASSERT( Util::sameString("Hello", cstrings[0]) );
  CPPUNIT_ASSERT( Util::sameString("am", cstrings[3]) );
  s = rf.nextLine();
  CPPUNIT_ASSERT ( *s == "This	row		has	a	null" );
  s = rf.nextRealLine();
  CPPUNIT_ASSERT ( *s == "here	is	the	third	line	of	input" );
}

/** Check and see if we get the next row in our test file.
*/
void RowFileTest::testNextRow()
{
  Verbose::out(1, "RowFileTest::testNextRow");
  RowFile rf;
  vector<string> words;
  rf.open(INPUT+"/RowFileTest.1.txt");
  rf.nextRow(words);
  CPPUNIT_ASSERT( words.size() == 4 );
  CPPUNIT_ASSERT( words[1] == "Here" );

  rf.nextRow(words);
  CPPUNIT_ASSERT( words.size() == 6 );
  CPPUNIT_ASSERT( words[1] == "row" );
  CPPUNIT_ASSERT( words[2] == "" );

  rf.nextRow(words);
  CPPUNIT_ASSERT( words.size() == 7 );
  CPPUNIT_ASSERT( words[1] == "is" );

  CPPUNIT_ASSERT(rf.nextRow(words) == false);
  rf.close();
}
//test nextLine() vliber
void RowFileTest::nextLineTest()
{
   Verbose::out(1, "RowFileTest::nextLineTest");
   RowFile rf;
   rf.open(INPUT+"/RowFileTest.1.txt");
   const string *line = NULL;
   line = rf.nextLine();   
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"Hello	Here	I	am"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==1);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"This	row		has	a	null"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==2);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"# This is a comment line"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==3);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),""));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==4);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"here	is	the	third	line	of	input"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==5);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),""));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==6);
   //eof
   line = rf.nextLine();
   CPPUNIT_ASSERT(line==NULL);
}
void RowFileTest::nextLineReuseTest()
{   
   Verbose::out(1, "RowFileTest::nextLineReuseTest");
   RowFile rf;
   rf.open(INPUT+"/RowFileTest.1.txt");
   const string *line = NULL;
   line = rf.nextLine();   
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"Hello	Here	I	am"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==1);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"This	row		has	a	null"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==2);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"# This is a comment line"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==3);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),""));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==4);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"here	is	the	third	line	of	input"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==5);
   
   //reuse
   rf.reuseLine();
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"here	is	the	third	line	of	input"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==5);
   line = rf.nextLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),""));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==6);
   //eof
   line = rf.nextLine();
   CPPUNIT_ASSERT(line==NULL);
}
//test nextRealLine() vliber
void RowFileTest::nextRealLineTest()
{
   Verbose::out(1, "RowFileTest::nextRealLineTest");
   RowFile rf;
   rf.open(INPUT+"/RowFileTest.1.txt");
   const string *line = NULL;
   line = rf.nextRealLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"Hello	Here	I	am"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==1);
   line = rf.nextRealLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"This	row		has	a	null"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==2);
   line = rf.nextRealLine();
   CPPUNIT_ASSERT(Util::sameString((*line).c_str(),"here	is	the	third	line	of	input"));
   CPPUNIT_ASSERT(rf.getCurrentLineNumber()==5);
   //eof
   CPPUNIT_ASSERT(rf.nextRealLine()==NULL);
 }

//test nextCStringRow() vliber
void RowFileTest::nextCStringRowExpectTest()
{
    Verbose::out(1, "RowFileTest::nextCStringRowExpectTest");
    RowFile rf1, rf2;
	vector<const char *> words;
    rf1.open(INPUT+"/RowFileTest.1.txt");
    CPPUNIT_ASSERT( rf1.nextCStringRowExpect(words,4));
	CPPUNIT_ASSERT(words.size()==4);
	words.clear();
	//negative
	//Goal is to receive a message: FATAL ERROR: Got 6 words instead of 5 expected at line: 2 
	NEGATIVE_TEST( rf1.nextCStringRowExpect(words,5),Except);
	CPPUNIT_ASSERT( rf1.nextCStringRowExpect(words,7));
	CPPUNIT_ASSERT(words.size()==7);
	    
	//additional empty cell at the end of the first row should be ignore
	words.clear();
	rf2.open(INPUT+"/RowFileTest.3.txt");
	CPPUNIT_ASSERT(rf2.nextCStringRowExpect(words,4));
	CPPUNIT_ASSERT(words.size()==4);
}	

//test nextCStringRow() vliber
void RowFileTest::nextCStringRowTest()
{
	Verbose::out(1, "RowFileTest::nextCStringRowTest");
    RowFile rf;
    vector<const char *> words;
	rf.setDelim(',');
    rf.open(INPUT+"/RowFileTest.2.txt");
	CPPUNIT_ASSERT(rf.nextCStringRow(words));
	CPPUNIT_ASSERT(words.size()==7);
	//eof
	CPPUNIT_ASSERT(rf.nextCStringRow(words)==false);
}

//test writeHeader(); readHeader() vliber
void RowFileTest::writeHeaderTest()
{
	Verbose::out(1, "RowFileTest::writeHeaderTest");
	filebuf fb;
	AffxString s=OUTPUT+"/header1.txt";
	const char * path1=s.c_str();
    fb.open (path1,ios::out);
    ostream output1 (&fb) ;
	RowFile rf1,rf2;
	vector<string> lines1, lines2;
	map<string,vector<string> > header1;
	lines1.push_back("#%this=one");
    lines1.push_back("#%this=two");
	lines1.push_back("#%this=3");
	lines1.push_back("#%this=");
	lines1.push_back("#%this=5");
	rf1.writeHeader(output1,lines1);
	//test	    
	rf2.open(OUTPUT+"/header1.txt");
    RowFile::readHeader(rf2, header1, lines2);
    CPPUNIT_ASSERT( header1["this"][0] == string("one"));
	CPPUNIT_ASSERT( header1["this"][1] == string("two"));
    CPPUNIT_ASSERT( header1["this"][2] == string("3"));
    CPPUNIT_ASSERT( header1["this"][3] == string("" ));
    CPPUNIT_ASSERT( header1["this"][4] == string("5"));
	fb.close();

	//negative
	AffxString s1=OUTPUT+"/header2.txt";
	const char * path2=s1.c_str();
	fb.open (path2,ios::out);
    ostream output2 (&fb) ;
	RowFile rf3,rf4;	
	vector<string> lines3,lines4;
	map<string,vector<string> > header2;
	lines3.push_back("#%this=one");
    lines3.push_back("#%this=two");
	lines3.push_back("#%this\\3");
	lines3.push_back("#%this=");
	lines3.push_back("#%this=5");
	rf3.writeHeader(output2,lines3);
	//test
	rf4.open(OUTPUT+"/header2.txt");
	//Goal is to receive a message: FATAL ERROR: Couldn't find delimiter: '=' in line: #%this\3
	NEGATIVE_TEST(RowFile::readHeader(rf4, header2, lines4),Except);
    CPPUNIT_ASSERT( header1["this"][0] == string("one"));
	CPPUNIT_ASSERT( header1["this"][1] == string("two"));
    CPPUNIT_ASSERT( header1["this"][3] == string("" ));
    CPPUNIT_ASSERT( header1["this"][4] == string("5"));
	fb.close();
}

//test open(); vliber
void RowFileTest::openTest()
{
	Verbose::out(1, "RowFileTest::openTest");
	RowFile rf1,rf2;
	POSITIVE_TEST(rf1.open(string(OUTPUT+"/header2.txt")));
	//negative
	//Goal is to receive a message: FATAL ERROR: Can't open file ./input1/header3.txt to read. -OK
	NEGATIVE_TEST(rf2.open(string("./input1/header3.txt")),Except);
}

//test matrixFromFile() vliber
void RowFileTest::matrixFromFileTest()
{
	Verbose::out(1, "RowFileTest::matrixFromFileTest");
	vector<vector<double> > inputData;
	RowFile::matrixFromFile(INPUT+"/matrix1.txt", inputData);
	vector< vector<double> >::iterator iter_ii;
	vector<double>::iterator iter_jj;
	for(iter_ii=inputData.begin(); iter_ii!=inputData.end(); iter_ii++)
	{      
		for(iter_jj=(*iter_ii).begin(); iter_jj!=(*iter_ii).end(); iter_jj++)
		{  
			CPPUNIT_ASSERT( inputData[0][0] == 1.0 );
			CPPUNIT_ASSERT( inputData[0][1] == 1.1 );
			CPPUNIT_ASSERT( inputData[0][2] == 1.2 );
			CPPUNIT_ASSERT( inputData[0][3] == 1.3 );
			CPPUNIT_ASSERT( inputData[1][0] == 2.0 );
			CPPUNIT_ASSERT( inputData[1][1] == 2.1 );
			CPPUNIT_ASSERT( inputData[1][2] == 2.2 );
			CPPUNIT_ASSERT( inputData[1][3] == 2.3 );
			CPPUNIT_ASSERT( inputData[2][0] == 3.0 );
			CPPUNIT_ASSERT( inputData[2][1] == 3.1 );
			CPPUNIT_ASSERT( inputData[2][2] == 3.2 );
			CPPUNIT_ASSERT( inputData[2][3] == 3.3 );
		} 
	}
	inputData.clear();
	RowFile::matrixFromFile(INPUT+"/matrix1.txt", inputData,1);
	for(iter_ii=inputData.begin(); iter_ii!=inputData.end(); iter_ii++)
	{      
		for(iter_jj=(*iter_ii).begin(); iter_jj!=(*iter_ii).end(); iter_jj++)
		{ 
			CPPUNIT_ASSERT( inputData[0][0] == 2.0 );
			CPPUNIT_ASSERT( inputData[0][1] == 2.1 );
			CPPUNIT_ASSERT( inputData[0][2] == 2.2 );
			CPPUNIT_ASSERT( inputData[0][3] == 2.3 );
			CPPUNIT_ASSERT( inputData[1][0] == 3.0 );
			CPPUNIT_ASSERT( inputData[1][1] == 3.1 );
			CPPUNIT_ASSERT( inputData[1][2] == 3.2 );
			CPPUNIT_ASSERT( inputData[1][3] == 3.3 );
		} 
	}
	inputData.clear();
	RowFile::matrixFromFile(INPUT+"/matrix1.txt", inputData,1,1);
	for(iter_ii=inputData.begin(); iter_ii!=inputData.end(); iter_ii++)
	{      
		for(iter_jj=(*iter_ii).begin(); iter_jj!=(*iter_ii).end(); iter_jj++)
		{ 
			CPPUNIT_ASSERT( inputData[0][0] == 2.1 );
			CPPUNIT_ASSERT( inputData[0][1] == 2.2 );
			CPPUNIT_ASSERT( inputData[0][2] == 2.3 );
			CPPUNIT_ASSERT( inputData[1][0] == 3.1 );
			CPPUNIT_ASSERT( inputData[1][1] == 3.2 );
			CPPUNIT_ASSERT( inputData[1][2] == 3.3 );
			//cout << *iter_jj << endl;
		} 
	}
	//negative
	//if column number to skip is > then number of columns
	inputData.clear();
	//Goal is to receive a message: FATAL ERROR: RowFile::matrixFromFile() - Number of skipCols >= number of cols.
    NEGATIVE_TEST(RowFile::matrixFromFile(INPUT+"/matrix1.txt", inputData,1,4),Except);
	//if row number to skip is > then number of rows
	inputData.clear();
    POSITIVE_TEST(RowFile::matrixFromFile(INPUT+"/matrix1.txt", inputData,3,1));
	CPPUNIT_ASSERT(inputData.empty());
}

