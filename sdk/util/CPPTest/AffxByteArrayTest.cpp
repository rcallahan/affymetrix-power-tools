////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "util/AffxArray.h"

#include "util/AffxString.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
#include "util/CPPTest/Setup.h"

using namespace std;

/**
 * @class AffxByteArrayTest
 * @brief cppunit class for testing AffxByteArray functions.
 * last change by vliber on 05/26/09
 */

class AffxByteArrayTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(AffxByteArrayTest);
    CPPUNIT_TEST(constructorTest);
	CPPUNIT_TEST(operatorTest);
	CPPUNIT_TEST(growingArrayFunctionTest);
	CPPUNIT_TEST(insertRemoveTest);
	CPPUNIT_TEST(getLineTest);
	CPPUNIT_TEST(trimTest);
	CPPUNIT_TEST(reverseComplementTest);
	CPPUNIT_TEST(firstWordTest);
	CPPUNIT_TEST(equalsTest);
	CPPUNIT_TEST(parameterTest);
	CPPUNIT_TEST(indexOfTest);
	CPPUNIT_TEST(compareToTest);
	CPPUNIT_TEST(getMaxRunTest);
	CPPUNIT_TEST(getCountInWindowTest);
	CPPUNIT_TEST_SUITE_END();

   public:
    void constructorTest();
	void operatorTest();
	void growingArrayFunctionTest();
	void insertRemoveTest();
	void getLineTest();
	void trimTest();
	void reverseComplementTest();
	void firstWordTest();
	void equalsTest();
	void parameterTest();
	void indexOfTest();
	void compareToTest();
	void getMaxRunTest();
	void getCountInWindowTest();
    
};
CPPUNIT_TEST_SUITE_REGISTRATION(AffxByteArrayTest );

void AffxByteArrayTest::constructorTest()
{
	cout<<endl;
	Verbose::out(1, "***AffxByteArrayTest testcases***");
    Verbose::out(1, "AffxByteArrayTest::constructorTest");
	//AffxByteArray();
	AffxByteArray ar1;
	CPPUNIT_ASSERT(ar1.length()==0);
	CPPUNIT_ASSERT(ar1.length()==ar1.getLength());
	CPPUNIT_ASSERT(ar1.getData()==0);
	CPPUNIT_ASSERT(ar1.getAllocLength()==0);
	CPPUNIT_ASSERT(ar1.isLocked()==false);
	
	//AffxByteArray(const AffxString& ba);
	AffxString s="acgtacgt";
	AffxByteArray ar2(s);
	CPPUNIT_ASSERT(ar2.getLength()==8);
	CPPUNIT_ASSERT(ar2.getAllocLength()==8);
	CPPUNIT_ASSERT(ar2.getUpperBound()==7);
	CPPUNIT_ASSERT(ar2.getAt(0)=='a');
	CPPUNIT_ASSERT(ar2.getAt(7)=='t');
    CPPUNIT_ASSERT(ar2.elementAt(0)==ar2.getAt(0));
	CPPUNIT_ASSERT(ar2.elementAt(8)==ar2.getAt(8));
	CPPUNIT_ASSERT(ar2.toString()=="acgtacgt");
	ar2.setAt(0,'T');
	ar2.setAt(7,'A');
	CPPUNIT_ASSERT(ar2.toString()=="TcgtacgA");
	
	//AffxByteArray(const AffxByteArray& ba);
	AffxByteArray ar3(ar2);
	CPPUNIT_ASSERT(ar3.getLength()==8);
	CPPUNIT_ASSERT(ar3.getAllocLength()==8);
	CPPUNIT_ASSERT(ar3.getUpperBound()==7);
	CPPUNIT_ASSERT(ar3.getAt(0)=='T');
	CPPUNIT_ASSERT(ar3.getAt(7)=='A');
    CPPUNIT_ASSERT(ar3.elementAt(0)==ar2.getAt(0));
	CPPUNIT_ASSERT(ar3.elementAt(7)==ar2.getAt(7));
	CPPUNIT_ASSERT(ar3.toString()=="TcgtacgA");
	ar3.setAt(0,'a');
	ar3.setAt(7,'t');
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgt");	
}
void AffxByteArrayTest::operatorTest()
{
	Verbose::out(1, "AffxByteArrayTest::operatorTest");
	AffxString s="acgtacgt";
	AffxByteArray ar1(s);
	CPPUNIT_ASSERT(ar1.toString()=="acgtacgt");
	

	//const AffxByteArray& AffxByteArray::operator=(AffxByteArray& ba)
	AffxByteArray ar2=ar1;
	CPPUNIT_ASSERT(ar2.getLength()==8);
	CPPUNIT_ASSERT(ar2.getAllocLength()==8);
	CPPUNIT_ASSERT(ar2.getUpperBound()==7);
	CPPUNIT_ASSERT(ar2.getAt(0)=='a');
	CPPUNIT_ASSERT(ar2.getAt(7)=='t');
    CPPUNIT_ASSERT(ar2.elementAt(0)==ar2.getAt(0));
	CPPUNIT_ASSERT(ar2.elementAt(7)==ar2.getAt(7));
	CPPUNIT_ASSERT(ar2.toString()=="acgtacgt");
	
	//bool operator==(AffxByteArray& ba)
	CPPUNIT_ASSERT(ar2.getLength()==ar1.getLength());
	for(int i=0;i<ar2.getLength();i++)
      CPPUNIT_ASSERT(ar2[i]==ar1[i]);
	CPPUNIT_ASSERT(ar1==ar2); 

	//char& operator[](int nIndex);
	ar2[0]='T';
	CPPUNIT_ASSERT(ar2.toString()=="Tcgtacgt");
	CPPUNIT_ASSERT(ar1!=ar2);

	//char operator[](int nIndex) const;
	CPPUNIT_ASSERT(ar2[7]=='t');
}

void AffxByteArrayTest::growingArrayFunctionTest()
{
	Verbose::out(1, "AffxByteArrayTest::growingArrayFunctionTest");
	AffxString s="acgtacgt";
	AffxByteArray ar1(s);
	//void setAtGrow(int nIndex, char newElement);
	ar1.setAtGrow(8,'A');
	CPPUNIT_ASSERT(ar1.elementAt(8)=='A');
	CPPUNIT_ASSERT(ar1.toString()=="acgtacgtA");
	ar1.setAtGrow(1,'A');
	CPPUNIT_ASSERT(ar1.toString()=="aAgtacgtA");
	ar1.setAtGrow(9,'T');
	CPPUNIT_ASSERT(ar1.toString()=="aAgtacgtAT");
	ar1.setAtGrow(11,'T');//todo vliber out of bounds
	CPPUNIT_ASSERT(ar1.toString()=="aAgtacgtAT");
    
	//int add(char newElement);
	AffxByteArray ar2(s);
	ar2.add('G');
	CPPUNIT_ASSERT(ar2.toString()=="acgtacgtG");
	ar2.add('G');
	CPPUNIT_ASSERT(ar2.toString()=="acgtacgtGG");
    

	//int append(const AffxString& src);
	AffxByteArray ar3(s);
	ar3.append("ACGT");
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGT");
	ar3.append(ar2);
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGTacgtacgtGG");
	
	//int append(int i)
    ar3.append(-100);
    CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGTacgtacgtGG-100");
	//int append(unsigned int i)
    ar3.append(400);
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGTacgtacgtGG-100400");
	//int AffxByteArray::append(char by)
    ar3.append('G');
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGTacgtacgtGG-100400G");
	//int append(double d)
	ar3.append(666600.98765);
	//cout<<ar3.toString()<<endl;
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtACGTacgtacgtGG-100400G666600.98765");


	//AffxByteArray& append(char bytes[], int offset, int len);
	AffxByteArray ar4(s);
	char in[]="GGTTCCAA";
	CPPUNIT_ASSERT(ar4.append(in,3,3).toString()=="acgtacgtTCC");
	CPPUNIT_ASSERT(ar4.append(in,5,3).toString()=="acgtacgtTCCCAA");
	CPPUNIT_ASSERT(ar4.append(in,5,2).toString()=="acgtacgtTCCCAACA");
	CPPUNIT_ASSERT(ar4.append(in,1,3).toString()=="acgtacgtTCCCAACAGTT");
	CPPUNIT_ASSERT(ar4.append(in,0,7).toString()=="acgtacgtTCCCAACAGTTGGTTCCA");
	CPPUNIT_ASSERT(ar4.append(in,0,7).toString()=="acgtacgtTCCCAACAGTTGGTTCCAGGTTCCA");
	

	//int AffxByteArray::append(const AffxByteArray& src)
    AffxByteArray ar5(s);
    ar5.append(ar1);
    CPPUNIT_ASSERT(ar5.toString()=="acgtacgtaAgtacgtAT");
    

	//void AffxByteArray::copy(const AffxByteArray& src)
	AffxByteArray ar6;
	ar6.copy(ar5);
	CPPUNIT_ASSERT(ar6.toString()=="acgtacgtaAgtacgtAT");
	

	//void assign(const AffxString& str) {setSize(0); append(str);}
    AffxByteArray ar7;
	ar7.assign("");
	CPPUNIT_ASSERT(ar7.toString()=="");
	ar7.assign("aaaaAAAA");
	CPPUNIT_ASSERT(ar7.toString()=="aaaaAAAA");
   

	//void copy(int iDestStart, AffxByteArray& src, int iSrcStart, int iLength);
     AffxByteArray ar8;
	 ar8.copy(0,ar7,0,8);
	 CPPUNIT_ASSERT(ar8.toString()=="aaaaAAAA");
	 ar8.copy(2,ar7,4,4);
	 CPPUNIT_ASSERT(ar8.toString()=="aaAAAAAA");
}
void AffxByteArrayTest::insertRemoveTest()
{
	Verbose::out(1, "AffxByteArrayTest::insertRemoveTest");
	AffxString s="acgtacgt";
	AffxByteArray ar1(s);
	//void insertAt(int nIndex, char newElement, int nCount = 1);
	ar1.insertAt(0,'C');
	CPPUNIT_ASSERT(ar1.toString()=="Cacgtacgt");
	ar1.insertAt(0,'G',4);
	CPPUNIT_ASSERT(ar1.toString()=="GGGGCacgtacgt");
	ar1.insertAt(3,'A',3);
	CPPUNIT_ASSERT(ar1.toString()=="GGGAAAGCacgtacgt");
	//void removeAt(int nIndex, int nCount = 1);
    ar1.removeAt(0);
	CPPUNIT_ASSERT(ar1.toString()=="GGAAAGCacgtacgt");
	ar1.removeAt(14);
	CPPUNIT_ASSERT(ar1.toString()=="GGAAAGCacgtacg");
	ar1.removeAt(0,2);
	CPPUNIT_ASSERT(ar1.toString()=="AAAGCacgtacg");
	ar1.removeAt(9,2);
	CPPUNIT_ASSERT(ar1.toString()=="AAAGCacgtg");
	ar1.removeAt(3,4);
	CPPUNIT_ASSERT(ar1.toString()=="AAAgtg");
	

	//void insertAt(int nStartIndex, AffxByteArray* pNewArray);
	AffxByteArray ar2(s);
	AffxByteArray ar3;
	ar3.insertAt(0,&ar2);
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgt");
	ar3.insertAt(8,&ar2);
	CPPUNIT_ASSERT(ar3.toString()=="acgtacgtacgtacgt");
	ar3.insertAt(3,&ar2);
	CPPUNIT_ASSERT(ar3.toString()=="acgacgtacgttacgtacgtacgt");
	ar3.insertAt(30,&ar2);//todo vliber out of bounds
	CPPUNIT_ASSERT(ar3.toString()=="acgacgtacgttacgtacgtacgt");
}


void AffxByteArrayTest::getLineTest()
{
	Verbose::out(1, "AffxByteArrayTest::getLineTest");
	AffxByteArray ba,line;
	CPPUNIT_ASSERT(ba.readFile(INPUT+"/test_1lqGold1.txt"));
	//first line
    CPPUNIT_ASSERT(ba.getLine(line));
	CPPUNIT_ASSERT(line.toString()=="X\tY\tSEQUENCE\tDESTYPE\tFEATURE\tQUALIFIER\tEXPOS\tPLEN\tPOS\tCBASE\tPBASE\tTBASE\tIPBASE\tUNIT\tBLOCK\tATOM");
	CPPUNIT_ASSERT(line.getColumnCount()==16);
	AffxByteArray column;
	line.getColumn(0,column);
	CPPUNIT_ASSERT(column.toString()=="");
	line.getColumn(1,column);
	CPPUNIT_ASSERT(column.toString()=="X");
	line.getColumn(16,column);
	CPPUNIT_ASSERT(column.toString()=="ATOM");
	//second line
	CPPUNIT_ASSERT(ba.getLine(line));
	CPPUNIT_ASSERT(line.getColumnCount()==16);
	line.getColumn(1,column);
	CPPUNIT_ASSERT(column.toString()=="0");
	CPPUNIT_ASSERT(column.parseInt()==0);
	line.getColumn(4,column);
	CPPUNIT_ASSERT(column.toString()=="-3.025");
	CPPUNIT_ASSERT(column.parseDouble()==-3.025);
	line.getColumn(10,column);
	CPPUNIT_ASSERT(column.toString()=="C");
	CPPUNIT_ASSERT(column.parseChar()=='C');
	line.getColumn(14,column);
	CPPUNIT_ASSERT(column.toString().toLowerCase()=="true");
	CPPUNIT_ASSERT(column.parsebool()==true);
	line.getColumn(16,column);
	CPPUNIT_ASSERT(column.toString()=="0");
	
	//bool nextLine(AffxByteArray& ba);
	//third line
    AffxByteArray nLine, col;
    CPPUNIT_ASSERT(ba.nextLine(nLine)==true);
	CPPUNIT_ASSERT(nLine.getColumnCount()==16);
	nLine.getColumn(1,col);
	CPPUNIT_ASSERT(col.toString()=="1");
	nLine.getColumn(3,col);
	CPPUNIT_ASSERT(col.toString()=="CAGCAGTTCTACGATGGCAAGTCCT");
	nLine.getColumn(6,col);
	CPPUNIT_ASSERT(col.toString()=="default_at");
	nLine.getColumn(7,col);
	CPPUNIT_ASSERT(col.toString()=="-1");
	//fourth line
	CPPUNIT_ASSERT(ba.nextLine(nLine)==true);
	CPPUNIT_ASSERT(nLine.getColumnCount()==16);
	nLine.getColumn(1,col);
	CPPUNIT_ASSERT(col.toString()=="2");
	nLine.getColumn(3,col);
	CPPUNIT_ASSERT(col.toString()=="AGGACTTGCCATCGTAGAACTGCTG");
	
	//AffxString getWord(int iWordIndex);
    CPPUNIT_ASSERT(col.toString()==nLine.getWord(3));
	CPPUNIT_ASSERT(nLine.getWord(3)=="AGGACTTGCCATCGTAGAACTGCTG");
	CPPUNIT_ASSERT(nLine.getWord(0)=="");
	CPPUNIT_ASSERT(nLine.getWord(17)=="");//out of bounds
	
	//AffxByteArray getWord(int iWordIndex, AffxByteArray& ba);
	AffxByteArray test;
    nLine.getWord(3,test);
    CPPUNIT_ASSERT(test.toString()=="AGGACTTGCCATCGTAGAACTGCTG");
	nLine.getWord(0,test);
	CPPUNIT_ASSERT(test.toString()=="");
	nLine.getWord(17,test);
	CPPUNIT_ASSERT(test.toString()=="");
   

	//AffxByteArray nextColumn(AffxByteArray& ba);
    AffxByteArray col1,col2,col3,col4;
    nLine.nextColumn(col1);
	CPPUNIT_ASSERT(col1.toString()=="2");
	nLine.nextColumn(col2);
	CPPUNIT_ASSERT(col2.toString()=="0");
	nLine.nextColumn(col3);
	CPPUNIT_ASSERT(col3.toString()=="AGGACTTGCCATCGTAGAACTGCTG");
	nLine.nextColumn(col4);
	CPPUNIT_ASSERT(col4.toString()=="3");
}
void AffxByteArrayTest::trimTest()
{
	Verbose::out(1, "AffxByteArrayTest::trimTest");
	AffxString s=" acgtacgt ";
	AffxByteArray ar1(s);
	//AffxByteArray& trim();
	CPPUNIT_ASSERT(ar1.trim().toString()=="acgtacgt");
	
	AffxByteArray& trimTabs();
	s="\tACGTacgt\t";
	AffxByteArray ar2(s);
	CPPUNIT_ASSERT(ar2.trimTabs().toString()=="ACGTacgt");
	
	AffxByteArray ar3(s);
	//AffxByteArray& toLowerCase();
    CPPUNIT_ASSERT(ar3.trimTabs().toLowerCase().toString()=="acgtacgt");
	//AffxByteArray& toUpperCase();
	CPPUNIT_ASSERT(ar3.toUpperCase().toString()=="ACGTACGT");
	CPPUNIT_ASSERT(ar3.toLowerCase().toString()=="acgtacgt");

	//bool startsWith(const AffxString& strCompare);
    CPPUNIT_ASSERT(ar3.startsWith("acgt"));
	CPPUNIT_ASSERT(!ar3.startsWith("ACGT"));

	//AffxString substring(int iIndex, int iEndIndex = -1);
	CPPUNIT_ASSERT(ar3.substring(3)=="tacgt");
	CPPUNIT_ASSERT(ar3.substring(0,5)=="acgta");
    CPPUNIT_ASSERT(ar3.substring(4,2)==""); //end<start
	CPPUNIT_ASSERT(ar3.substring(20,2)==""); //out of bounds
	ar3[8]='\0';
	CPPUNIT_ASSERT(ar3.substring(4,20)=="acgt");
   
	//void AffxByteArray::trimInternal()
	s="\tacgt\tGTCA\tatgc\t";
	AffxByteArray ar4(s);
	ar4.trimInternal();
	CPPUNIT_ASSERT(ar4.toString()=="acgtGTCAatgc");
	s="\tacgt GTCA\tatgc\t";
	AffxByteArray ar5(s);
	ar5.trimInternal();
	CPPUNIT_ASSERT(ar4.toString()=="acgtGTCAatgc");
}


void AffxByteArrayTest::reverseComplementTest()
{
	Verbose::out(1, "AffxByteArrayTest::reverseComplementTest");
	AffxString s="acgtaTCG";//even
	AffxByteArray ar(s);	
	CPPUNIT_ASSERT(ar.reverseComplement().toString()=="cgatacgt");
	s="";
	AffxByteArray ar1(s);	
	CPPUNIT_ASSERT(ar1.reverseComplement().toString()=="");
	s="afgvaTCG";
	AffxByteArray ar2(s);
	CPPUNIT_ASSERT(ar2.reverseComplement().toString()=="cgatbcft");
	s="aggcAaTCG";//odd
	AffxByteArray ar3(s);
	//CPPUNIT_ASSERT(ar.reverseComplement().toString()=="cgattgcct");
}


void AffxByteArrayTest::firstWordTest()
{
	Verbose::out(1, "AffxByteArrayTest::firstWordTest");
	AffxString s="acgt aTCG\tgtca";//even
	AffxByteArray ar(s);	
	CPPUNIT_ASSERT(ar.firstWord(0)=="acgt");
	CPPUNIT_ASSERT(ar.firstWord(1)=="cgt");
	CPPUNIT_ASSERT(ar.firstWord(3)=="t");
	CPPUNIT_ASSERT(ar.firstWord(4)=="aTCG");
	CPPUNIT_ASSERT(ar.firstWord(6)=="TCG");
	CPPUNIT_ASSERT(ar.firstWord(9)=="gtca");
	CPPUNIT_ASSERT(ar.firstWord(13)=="a");
	CPPUNIT_ASSERT(ar.firstWord(14)=="");
}

void AffxByteArrayTest::equalsTest()
{
	Verbose::out(1, "AffxByteArrayTest::equalsTest");
	AffxString s="acgt aTCG\tgtca";
	AffxByteArray ar(s);
	CPPUNIT_ASSERT(!ar.equals("acgt aaCG\tgtca"));
	CPPUNIT_ASSERT(ar.equalsIgnoreCase("acgt aTCG\tgtca"));
	s="";
	AffxByteArray ar1(s);
	CPPUNIT_ASSERT(ar1.equals(""));
}


void AffxByteArrayTest::parameterTest()
{
	Verbose::out(1, "AffxByteArrayTest::parameterTest");
	AffxByteArray ba,line;
	CPPUNIT_ASSERT(ba.readFile("./input/patameterFile.txt"));
	CPPUNIT_ASSERT(ba.parameterCount()==40);
	ba.nextLine(line);
	CPPUNIT_ASSERT(line.parameterCount()==4);
    CPPUNIT_ASSERT(line.toString()=="; start input properties");
	ba.nextLine(line);
	CPPUNIT_ASSERT(line.parameterCount()==0);
	CPPUNIT_ASSERT(line.getParameter(1).toString()=="");
	ba.nextLine(line);
	CPPUNIT_ASSERT(line.parameterCount()==2);
	CPPUNIT_ASSERT(line.getParameter(1).toString()=="LogFile:");
    CPPUNIT_ASSERT(line.getParameter(2).toString()=="O:/vliber/DesignRequestChecker/4.1/data/period/<SW:ARRAYNAME>/ExpressionDesignRequestChecker_<SW:ARRAYNAME>_<SW:TIMESTAMP>.log");
	ba.nextLine(line);
	ba.nextLine(line);
	ba.nextLine(line);
	ba.nextLine(line);
	CPPUNIT_ASSERT(line.parameterCount()==2);
	CPPUNIT_ASSERT(line.getParameter(1).toString()=="ArrayFormat:");
	CPPUNIT_ASSERT(line.getParameter(2).toString()=="049");
}
void AffxByteArrayTest::indexOfTest()
{
	Verbose::out(1, "AffxByteArrayTest::indexOfTest");
	AffxString s="acgt aTCG\tgtca";
	AffxByteArray ar(s);
	CPPUNIT_ASSERT(ar.indexOf("u")==-1);
	CPPUNIT_ASSERT(ar.indexOf("a")==0);
	CPPUNIT_ASSERT(ar.indexOf(" a")==4);
	CPPUNIT_ASSERT(ar.indexOf("\tg")==9);
	CPPUNIT_ASSERT(ar.indexOf("T")==6);
	CPPUNIT_ASSERT(ar.indexOf("t")==3);
}


void AffxByteArrayTest::compareToTest()
{
	Verbose::out(1, "AffxByteArrayTest::compareToTest");
	AffxString s="acgtACGT";
	AffxByteArray ar(s);
	//int compareTo(AffxString that);
	CPPUNIT_ASSERT(ar.compareTo("acgt")>0);
	CPPUNIT_ASSERT(ar.compareTo("acgtacgtacgt")<0);
	CPPUNIT_ASSERT(ar.compareTo("acgtACGT")==0);
	CPPUNIT_ASSERT(ar.compareTo("AcgtACGT")>0);
	CPPUNIT_ASSERT(ar.compareTo("bcgtACGT")<0);

	//int compareTo(AffxByteArray that);
	AffxString s1="acgt";
    AffxByteArray ar1(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar1)>0);
	s1="acgtacgtacgt";
    AffxByteArray ar2(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar2)<0);
	s1="acgtACGT";
    AffxByteArray ar3(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar3)==0);
	s1="AcgtACGT";
    AffxByteArray ar4(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar4)>0);
	s1="bcgtACGT";
    AffxByteArray ar5(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar5)<0);

	//int compareTo(AffxByteArray obj, int iCompareCode);
	s1="acgt";
    AffxByteArray ar6(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar6,0)>0);
	s1="acgtacgtacgt";
    AffxByteArray ar7(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar7,0)<0);
	s1="acgtACGT";
    AffxByteArray ar8(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar8,0)==0);
	s1="AcgtACGT";
    AffxByteArray ar9(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar9,0)>0);
	s1="bcgtACGT";
    AffxByteArray ar10(s1);
	CPPUNIT_ASSERT(ar.compareTo(ar10,0)<0);
	//
	s1="124";
    AffxByteArray ari(s1);
	AffxString s2="124";
    AffxByteArray ari1(s2);
	CPPUNIT_ASSERT(ari.compareTo(ari1,1)==0);
	s2="125";
    AffxByteArray ari2(s2);
	CPPUNIT_ASSERT(ari.compareTo(ari2,1)==-1);
	s2="105";
    AffxByteArray ari3(s2);
	CPPUNIT_ASSERT(ari.compareTo(ari3,1)==1);
}


void AffxByteArrayTest::getMaxRunTest()
{
	Verbose::out(1, "AffxByteArrayTest::getMaxRunTest");
	AffxString s="aaaaa";
	AffxByteArray ar1(s);
	CPPUNIT_ASSERT(ar1.getMaxRun('c')==0);
	CPPUNIT_ASSERT(ar1.getMaxRun('a')==5);
	s="aaca";
	AffxByteArray ar2(s);
	CPPUNIT_ASSERT(ar2.getMaxRun('a')==2);
	s="aacaaa";
	AffxByteArray ar3(s);
	CPPUNIT_ASSERT(ar3.getMaxRun('a')==3);
	s="aaccccaaa";
	AffxByteArray ar4(s);
	CPPUNIT_ASSERT(ar4.getMaxRun('c')==4);
	s="acgtaaccggttaaacccgggtttaaaa";
	AffxByteArray ar5(s);
	CPPUNIT_ASSERT(ar5.getMaxRun('a')==4);
}


void AffxByteArrayTest::getCountInWindowTest()
{
	Verbose::out(1, "AffxByteArrayTest::getCountInWindowTest");
	AffxString s="aaaaa";
	AffxByteArray ar1(s);
	CPPUNIT_ASSERT(ar1.getCountInWindow('c',5)==0);
	CPPUNIT_ASSERT(ar1.getCountInWindow('a',4)==4);
	s="aaca";
	AffxByteArray ar2(s);
	CPPUNIT_ASSERT(ar2.getCountInWindow('a',4)==3);
	s="aacaaa";
	AffxByteArray ar3(s);
	CPPUNIT_ASSERT(ar3.getCountInWindow('a',4)==3);
	s="aaccccaaa";
	AffxByteArray ar4(s);
	CPPUNIT_ASSERT(ar4.getCountInWindow('c',6)==4);
	s="acgtaaccggttaaacccgggtttaaaa";
	AffxByteArray ar5(s);
	CPPUNIT_ASSERT(ar5.getCountInWindow('a',4)==4);
	s="acgtggctttaaggctt";
	AffxByteArray ar6(s);
	CPPUNIT_ASSERT(ar6.getCountInWindow('t',4)==3);
	CPPUNIT_ASSERT(ar6.getCountInWindow('t',10)==5);
	CPPUNIT_ASSERT(ar6.getCountInWindow('t',14)==6);
}
