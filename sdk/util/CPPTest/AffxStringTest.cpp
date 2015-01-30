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
 * @class AffxStringTest
 * @brief cppunit class for testing AffxString functions.
 * last change by vliber on 02/25/08
 */

class AffxStringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(AffxStringTest);
  CPPUNIT_TEST(isTabTest);
  CPPUNIT_TEST(isWhiteSpaceTest);
  CPPUNIT_TEST(isAlphabeticTest);
  CPPUNIT_TEST(isAlphanumericTest);
  CPPUNIT_TEST(isBlankIndexTest);
  CPPUNIT_TEST(isBlankTest);
  CPPUNIT_TEST(isBlankOrEmptyTest);
  CPPUNIT_TEST(isPunctuationTest);
  CPPUNIT_TEST(isUppercaseTest);
  CPPUNIT_TEST(isLowercaseTest);
  CPPUNIT_TEST(isNumericTest);
  CPPUNIT_TEST(isHexNumericTest);
  CPPUNIT_TEST(toUpperTest);
  CPPUNIT_TEST(toLowerTest);
  CPPUNIT_TEST(getLengthTest);
  CPPUNIT_TEST(midTest);
  CPPUNIT_TEST(stripPrecedingBlanksTest);
  CPPUNIT_TEST(stripTrailingBlanksTest);
  CPPUNIT_TEST(padBlanksTest);
  CPPUNIT_TEST(stripAndPadTest);
  CPPUNIT_TEST(stripTest);
  CPPUNIT_TEST(trimTest);
  CPPUNIT_TEST(intToStringTest);
  CPPUNIT_TEST(doubleToStringTest);
  CPPUNIT_TEST(setAtTest);
  CPPUNIT_TEST(getAtTest);
  CPPUNIT_TEST(charAtTest);
  CPPUNIT_TEST(toLowerCaseTest);
  CPPUNIT_TEST(toUpperCaseTest);
  CPPUNIT_TEST(substringTest);
  CPPUNIT_TEST(replaceTest);
  CPPUNIT_TEST(equalsTest);
  CPPUNIT_TEST(reverseComplementTest);
  CPPUNIT_TEST(reverseTest);
  CPPUNIT_TEST(complementTest);
  CPPUNIT_TEST(convertAmbiguitiesTest);
  CPPUNIT_TEST(startsWithTest);
  CPPUNIT_TEST(compareToTest);
  CPPUNIT_TEST(endsWithTest);
  CPPUNIT_TEST(indexOfTest);
  CPPUNIT_TEST(nextIndexOfTest);
  CPPUNIT_TEST(lastIndexOfTest);
  CPPUNIT_TEST(operatorPlusTest);
  CPPUNIT_TEST(removeSurroundingQuotesTest);
  CPPUNIT_TEST(isBaseAtTest);
  CPPUNIT_TEST(isSnpAtTest);
  CPPUNIT_TEST(isIUPACCodeAtTest);
  CPPUNIT_TEST(getSnpStringAtTest);
  CPPUNIT_TEST(hashTest);
  CPPUNIT_TEST(constructorTest);
  CPPUNIT_TEST_SUITE_END();

public:
  
  void isTabTest();
  void isWhiteSpaceTest();
  void isAlphabeticTest();
  void isAlphanumericTest();
  void isBlankIndexTest();
  void isBlankTest();
  void isBlankOrEmptyTest();
  void isPunctuationTest();
  void isUppercaseTest();
  void isLowercaseTest();
  void isNumericTest();
  void isHexNumericTest();
  void toUpperTest();
  void toLowerTest(); 
  void getLengthTest();
  void midTest();
  void stripPrecedingBlanksTest();
  void stripTrailingBlanksTest();
  void padBlanksTest();
  void stripAndPadTest();
  void stripTest();
  void trimTest();
  void intToStringTest();
  void doubleToStringTest();
  void setAtTest();
  void getAtTest();
  void charAtTest();
  void toLowerCaseTest();
  void toUpperCaseTest();
  void substringTest();
  void replaceTest();
  void equalsTest();
  void reverseComplementTest();
  void reverseTest();
  void complementTest();
  void convertAmbiguitiesTest();
  void startsWithTest();
  void compareToTest();
  void endsWithTest();
  void indexOfTest();
  void nextIndexOfTest();
  void lastIndexOfTest();
  void operatorPlusTest();
  void removeSurroundingQuotesTest();
  void isBaseAtTest();
  void isSnpAtTest();
  void isIUPACCodeAtTest();
  void getSnpStringAtTest();
  void hashTest();
  void constructorTest();
};


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(AffxStringTest );



void AffxStringTest::isTabTest()
{
	cout<<endl;
	cout<<endl;
	Verbose::out(1, "***AffxStringTest testcases***");
    Verbose::out(1, "AffxStringTest::isTabTest");
	AffxString s="	ac&gt	a cgt	";
	CPPUNIT_ASSERT(s.isTab(0));
	CPPUNIT_ASSERT(!s.isTab(1));
	CPPUNIT_ASSERT(!s.isTab(3));
	CPPUNIT_ASSERT(s.isTab(6));
	CPPUNIT_ASSERT(!s.isTab(8));
	CPPUNIT_ASSERT(s.isTab(12));
	NEGATIVE_TEST(s.isTab(13),Except);
}

void AffxStringTest::isWhiteSpaceTest()
{
	Verbose::out(1, "AffxStringTest::isWhiteSpaceTest");
	AffxString s="	ac&gt acgt	\n";
	CPPUNIT_ASSERT(s.isWhiteSpace(0)!=0);
	CPPUNIT_ASSERT(s.isWhiteSpace(1)==0);
	CPPUNIT_ASSERT(s.isWhiteSpace(6)!=0);
	CPPUNIT_ASSERT(s.isWhiteSpace(11)!=0);
	CPPUNIT_ASSERT(s.isWhiteSpace(12)!=0);
	NEGATIVE_TEST(s.isWhiteSpace(13),Except);//out of length
}

void AffxStringTest::isAlphabeticTest()
{
	Verbose::out(1, "AffxStringTest::isAlphabeticTest");
	AffxString s="	az&Gt aAgZ\n";
	CPPUNIT_ASSERT(s.isAlphabetic(0)==0);
	CPPUNIT_ASSERT(s.isAlphabetic(1)!=0);
	CPPUNIT_ASSERT(s.isAlphabetic(2)!=0);
	CPPUNIT_ASSERT(s.isAlphabetic(6)==0);
	CPPUNIT_ASSERT(s.isAlphabetic(8)!=0);
	CPPUNIT_ASSERT(s.isAlphabetic(10)!=0);
	CPPUNIT_ASSERT(s.isAlphabetic(11)==0);
	NEGATIVE_TEST(s.isAlphabetic(12),Except);//out of length
}

void AffxStringTest::isAlphanumericTest()
{
	Verbose::out(1, "AffxStringTest::isAlphanumericTest");
	AffxString s="	az9Gt 0AgZ\n";
	CPPUNIT_ASSERT(s.isAlphanumeric(0)==0);
	CPPUNIT_ASSERT(s.isAlphanumeric(1)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(2)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(3)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(6)==0);
	CPPUNIT_ASSERT(s.isAlphanumeric(7)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(8)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(10)!=0);
	CPPUNIT_ASSERT(s.isAlphanumeric(11)==0);
	NEGATIVE_TEST(s.isAlphanumeric(12),Except);//out of length
}

void AffxStringTest::isBlankIndexTest()
{
	Verbose::out(1, "AffxStringTest::isBlankIndexTest");
	AffxString s=" ac&gt	a cgt\n ";
	CPPUNIT_ASSERT(s.isBlank(0));
	CPPUNIT_ASSERT(!s.isBlank(1));
	CPPUNIT_ASSERT(!s.isBlank(3));
	CPPUNIT_ASSERT(!s.isBlank(6));
	CPPUNIT_ASSERT(s.isBlank(8));
	CPPUNIT_ASSERT(!s.isBlank(12));
	CPPUNIT_ASSERT(s.isBlank(13));
	NEGATIVE_TEST(s.isBlank(14),Except);//out of length
}

void AffxStringTest::isBlankTest()
{
	Verbose::out(1, "AffxStringTest::isBlankTest");
	AffxString s=" ac&gt	a cgt\n ";
	CPPUNIT_ASSERT(!s.isBlank());
	s="         ";
	CPPUNIT_ASSERT(s.isBlank());
	s="";
	CPPUNIT_ASSERT(s.isBlank());
}

void AffxStringTest::isBlankOrEmptyTest()
{
	Verbose::out(1, "AffxStringTest::isBlankOrEmptyTest");
	AffxString s=" ac&gt	a cgt\n ";
	CPPUNIT_ASSERT(!s.isBlankOrEmpty());
	s="         ";
	CPPUNIT_ASSERT(s.isBlankOrEmpty());
	s="";
	CPPUNIT_ASSERT(s.isBlankOrEmpty());
}

void AffxStringTest::isPunctuationTest()
{
	Verbose::out(1, "AffxStringTest::isPunctuationTest");
	AffxString s="	az$Gt _A-#\n";
	CPPUNIT_ASSERT(s.isPunctuation(0)==0);
	CPPUNIT_ASSERT(s.isPunctuation(11)==0);
	CPPUNIT_ASSERT(s.isPunctuation(1)==0);
	CPPUNIT_ASSERT(s.isPunctuation(3)!=0);
	CPPUNIT_ASSERT(s.isPunctuation(7)!=0);
	CPPUNIT_ASSERT(s.isPunctuation(8)==0);
	CPPUNIT_ASSERT(s.isPunctuation(9)!=0);
	CPPUNIT_ASSERT(s.isPunctuation(10)!=0);
	CPPUNIT_ASSERT(s.isPunctuation(11)==0);
	NEGATIVE_TEST(s.isPunctuation(12),Except);//out of length
}

void AffxStringTest::isUppercaseTest()
{
	Verbose::out(1, "AffxStringTest::isUppercaseTest");
	AffxString s="	az$Gt _A-#\n";
	CPPUNIT_ASSERT(s.isUppercase(0)==0);
	CPPUNIT_ASSERT(s.isUppercase(1)==0);
	CPPUNIT_ASSERT(s.isUppercase(3)==0);
	CPPUNIT_ASSERT(s.isUppercase(4)!=0);
	CPPUNIT_ASSERT(s.isUppercase(8)!=0);
	CPPUNIT_ASSERT(s.isUppercase(9)==0);
	CPPUNIT_ASSERT(s.isUppercase(10)==0);
	CPPUNIT_ASSERT(s.isUppercase(11)==0);
	NEGATIVE_TEST(s.isUppercase(12),Except);//out of length
	
}

void AffxStringTest::isLowercaseTest()
{
	Verbose::out(1, "AffxStringTest::isLowercaseTest");
	AffxString s="	az$Gt _A-#\n";
	CPPUNIT_ASSERT(s.isLowercase(0)==0);
	CPPUNIT_ASSERT(s.isLowercase(1)!=0);
	CPPUNIT_ASSERT(s.isLowercase(2)!=0);
	CPPUNIT_ASSERT(s.isLowercase(5)!=0);
	CPPUNIT_ASSERT(s.isLowercase(8)==0);
	CPPUNIT_ASSERT(s.isLowercase(9)==0);
	CPPUNIT_ASSERT(s.isLowercase(10)==0);
	CPPUNIT_ASSERT(s.isLowercase(11)==0);
	NEGATIVE_TEST(s.isLowercase(12),Except);//out of length
}

void AffxStringTest::isNumericTest()
{
	Verbose::out(1, "AffxStringTest::isNumericTest");
	AffxString s="	a0$G9 \n";
	CPPUNIT_ASSERT(s.isNumeric(0)==0);
	CPPUNIT_ASSERT(s.isNumeric(2)!=0);
	CPPUNIT_ASSERT(s.isNumeric(3)==0);
	CPPUNIT_ASSERT(s.isNumeric(5)!=0);
	CPPUNIT_ASSERT(s.isNumeric(7)==0);
	NEGATIVE_TEST(s.isNumeric(8),Except);//out of length
}

void AffxStringTest::isHexNumericTest()
{
	Verbose::out(1, "AffxStringTest::isHexNumericTest");
	AffxString s="0Ag0f9Ft8a";
	CPPUNIT_ASSERT(s.isHexNumeric(0)!=0); 
	CPPUNIT_ASSERT(s.isHexNumeric(1)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(2)==0);
	CPPUNIT_ASSERT(s.isHexNumeric(3)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(4)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(5)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(6)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(7)==0);
	CPPUNIT_ASSERT(s.isHexNumeric(8)!=0);
	CPPUNIT_ASSERT(s.isHexNumeric(9)!=0);
	NEGATIVE_TEST(s.isHexNumeric(10),Except);//out of length
	
}

void AffxStringTest::toUpperTest()
{
	Verbose::out(1, "AffxStringTest::toUpperTest");
	AffxString s =" acgT acgt ";
	s.toUpper(0);
	s.toUpper(2);
	s.toUpper(2);
	s.toUpper(4);
	s.toUpper(5);
	s.toUpper(7);
	s.toUpper(10);
    CPPUNIT_ASSERT(s==" aCgT aCgt "); 
}

void AffxStringTest::toLowerTest()
{
	Verbose::out(1, "AffxStringTest::toLowerTest");
	AffxString s =" ACGt ACGT ";
	s.toLower(0);
	s.toLower(2);
	s.toLower(2);
	s.toLower(4);
	s.toLower(5);
	s.toLower(7);
	s.toLower(10);
    CPPUNIT_ASSERT(s==" AcGt AcGT "); 
	
}

void AffxStringTest::getLengthTest()
{
	Verbose::out(1, "AffxStringTest::getLengthTest");
	AffxString s ="";
	CPPUNIT_ASSERT(s.getLength()==0); 
	s ="  ";
	CPPUNIT_ASSERT(s.getLength()==2); 
	s ="	";
	CPPUNIT_ASSERT(s.getLength()==1);
	s ="\t";
	CPPUNIT_ASSERT(s.getLength()==1);
	s ="\t\t\t";
	CPPUNIT_ASSERT(s.getLength()==3);
	s =" asdf fght ";
	CPPUNIT_ASSERT(s.getLength()==11);
}

void AffxStringTest::midTest()
{
	Verbose::out(1, "AffxStringTest::midTest");
	AffxString s ="";
	CPPUNIT_ASSERT(s.mid(0,5)=="");
	s ="     ";
	CPPUNIT_ASSERT(s.mid(0,5)=="     "); 
	s =" acgttgca ACGTTGCA ";
	CPPUNIT_ASSERT(s.mid(0,5)==" acgt");
	CPPUNIT_ASSERT(s.mid(0,19)==" acgttgca ACGTTGCA ");
	CPPUNIT_ASSERT(s.mid(0,29)==" acgttgca ACGTTGCA ");
	CPPUNIT_ASSERT(s.mid(6,4)=="gca ACGTTG");
}

void AffxStringTest::stripPrecedingBlanksTest()
{
	Verbose::out(1, "AffxStringTest::stripPrecedingBlanksTest");
	AffxString s ="     jhgdjaghdjasg";
	s.stripPrecedingBlanks();
	CPPUNIT_ASSERT(s=="jhgdjaghdjasg");
	s="jhgdj     aghdjasg";
	s.stripPrecedingBlanks();
	CPPUNIT_ASSERT(s=="jhgdj     aghdjasg");
	s="jhgdj     ";
	s.stripPrecedingBlanks();
	CPPUNIT_ASSERT(s=="jhgdj     ");
}

void AffxStringTest::stripTrailingBlanksTest()
{
	Verbose::out(1, "AffxStringTest::stripTrailingBlanksTest");
	AffxString s ="     jhgdjaghdjasg";
	s.stripTrailingBlanks();
	CPPUNIT_ASSERT(s=="     jhgdjaghdjasg");
	s="jhgdj     aghdjasg";
	s.stripTrailingBlanks();
	CPPUNIT_ASSERT(s=="jhgdj     aghdjasg");
	s="jhgdj     ";
	s.stripTrailingBlanks();
	CPPUNIT_ASSERT(s=="jhgdj");
}

void AffxStringTest::padBlanksTest()
{
	Verbose::out(1, "AffxStringTest::padBlanksTest");
	AffxString s="acgt"; //length=4
	s.padBlanks(9);//additional 5 spaces
	CPPUNIT_ASSERT(s=="acgt     ");
	s="";
	s.padBlanks(5);//5 spaces
	CPPUNIT_ASSERT(s=="     ");
	s="acgta";
	s.padBlanks(5);//identical
	CPPUNIT_ASSERT(s=="acgta");
	s.padBlanks(0);//identical
	CPPUNIT_ASSERT(s=="acgta");
}

void AffxStringTest::stripAndPadTest()
{
	Verbose::out(1, "AffxStringTest::stripAndPadTest");
	AffxString s="   acgt    "; 
	s.stripAndPad(12);
	CPPUNIT_ASSERT(s=="acgt        ");
	s="   acgt"; 
	s.stripAndPad(12);
	CPPUNIT_ASSERT(s=="acgt        ");
	s="acgt    "; 
	s.stripAndPad(12);
	CPPUNIT_ASSERT(s=="acgt        ");
	s="acgt    "; 
	s.stripAndPad(4);
	CPPUNIT_ASSERT(s=="acgt");
	s=""; 
	s.stripAndPad(4);
	CPPUNIT_ASSERT(s=="    ");
}

void AffxStringTest::stripTest()
{
	Verbose::out(1, "AffxStringTest::stripTest");
	AffxString s="   acgt    "; 
	s.strip();
	CPPUNIT_ASSERT(s=="acgt");
	s="   acgt"; 
	s.strip();
	CPPUNIT_ASSERT(s=="acgt");
	s="acgt    "; 
	s.strip();
	CPPUNIT_ASSERT(s=="acgt");
	s=""; 
	s.strip();
	CPPUNIT_ASSERT(s=="");
}

void AffxStringTest::trimTest()
{
	Verbose::out(1, "AffxStringTest::trimTest");
	AffxString s="   acgt    "; 
	CPPUNIT_ASSERT(s.trim()=="acgt");
	s="   acgt"; 
	CPPUNIT_ASSERT(s.trim()=="acgt");
	s="acgt    "; 
	CPPUNIT_ASSERT(s.trim()=="acgt");
	s=""; 
	CPPUNIT_ASSERT(s.trim()=="");
}

void AffxStringTest::intToStringTest()
{
	Verbose::out(1, "AffxStringTest::intToStringTest");
	CPPUNIT_ASSERT(AffxString::intToString(123456789,true)=="123,456,789");
	CPPUNIT_ASSERT(AffxString::intToString(123456789,false)=="123456789");
	CPPUNIT_ASSERT(AffxString::intToString(-123456789,true)=="-123,456,789");
	CPPUNIT_ASSERT(AffxString::intToString(-123456789,false)=="-123456789");
}

void AffxStringTest::doubleToStringTest()
{
	Verbose::out(1, "AffxStringTest::doubleToStringTest");
	#ifdef WIN32
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13145,3,true)=="123,456,789.131");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13155,3,true)=="123,456,789.132");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13165,3,true)=="123,456,789.132");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.1300,4,true)=="123,456,789.13");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.134,2,false)=="123456789.13");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789,3,false)=="123456789");
	CPPUNIT_ASSERT(AffxString::doubleToString(-123456789.13157,4,true)=="-123,456,789.1316");
    #else
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13145,3,true)=="123,456,789.131");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13155,3,true)=="123,456,789.132");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.13165,3,true)=="123,456,789.132");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.1300,4,true)=="123,456,789.13");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789.134,2,false)=="123456789.13");
	CPPUNIT_ASSERT(AffxString::doubleToString(123456789,3,false)=="123456789");
	CPPUNIT_ASSERT(AffxString::doubleToString(-123456789.13157,4,true)=="-123,456,789.1316");
    #endif
}

void AffxStringTest::setAtTest()
{
	Verbose::out(1, "AffxStringTest::setAtTest");
	AffxString s="acgt ACGT";
	s.setAt(0,'T');
	CPPUNIT_ASSERT(s=="Tcgt ACGT");
	s.setAt(4,'T');
	CPPUNIT_ASSERT(s=="TcgtTACGT");
	s.setAt(0,'c');
	CPPUNIT_ASSERT(s=="ccgtTACGT");
	s.setAt(8,'a');
	CPPUNIT_ASSERT(s=="ccgtTACGa");
	NEGATIVE_TEST(s.setAt(9,'c'),Except);//Out of Bounds exception
}

void AffxStringTest::getAtTest()
{
	Verbose::out(1, "AffxStringTest::getAtTest");
	AffxString s="acgt ACGT";
	CPPUNIT_ASSERT(ToStr(s.getAt(0))=="a");
	CPPUNIT_ASSERT(ToStr(s.getAt(4))==" ");
	CPPUNIT_ASSERT(ToStr(s.getAt(8))=="T");
	NEGATIVE_TEST(s.getAt(9),Except);//Out of Bounds exception
}
void AffxStringTest::charAtTest()
{
	Verbose::out(1, "AffxStringTest::charAtTest");
	AffxString s="acgt ACGT";
	CPPUNIT_ASSERT(ToStr(s.charAt(0))=="a");
	CPPUNIT_ASSERT(ToStr(s.charAt(4))==" ");
	CPPUNIT_ASSERT(ToStr(s.charAt(8))=="T");
	NEGATIVE_TEST(s.charAt(9),Except);//Out of Bounds exception
}
void AffxStringTest::toLowerCaseTest()
{
	Verbose::out(1, "AffxStringTest::toLowerCaseTest");
	AffxString s="ACGT AcGt";
	CPPUNIT_ASSERT(s.toLowerCase()=="acgt acgt");
	s=" ";
	CPPUNIT_ASSERT(s.toLowerCase()==" ");
}
void AffxStringTest::toUpperCaseTest()
{
	Verbose::out(1, "AffxStringTest::toUpperCaseTest");
	AffxString s="acGt acGt";
	CPPUNIT_ASSERT(s.toUpperCase()=="ACGT ACGT");
	s=" ";
	CPPUNIT_ASSERT(s.toUpperCase()==" ");
}
void AffxStringTest::substringTest()
{
	Verbose::out(1, "AffxStringTest::substringTest");;
	AffxString s="acGt acGt";
	CPPUNIT_ASSERT(s.substring(0)=="acGt acGt");
	CPPUNIT_ASSERT(s.substring(4)==" acGt");
	CPPUNIT_ASSERT(s.substring(8)=="t");
	CPPUNIT_ASSERT(s.substring(9)=="");
	//cout<<s.substring(0,4)<<endl;
	//cout<<s.substring(3,4)<<endl;
    CPPUNIT_ASSERT(s.substring(0,4)=="acGt"); 
	CPPUNIT_ASSERT(s.substring(3,4)=="t");   
    CPPUNIT_ASSERT(s.substring(0,9)=="acGt acGt");    
    CPPUNIT_ASSERT(s.substring(9,2)=="");
    CPPUNIT_ASSERT(s.substring(0,10)=="acGt acGt");
}
void AffxStringTest::replaceTest()
{
	Verbose::out(1, "AffxStringTest::replaceTest");
	AffxString s="acGt acGt";
	CPPUNIT_ASSERT(s.replace('a','A')=="AcGt AcGt");
	CPPUNIT_ASSERT(s.replace('A','T')=="TcGt TcGt");
	CPPUNIT_ASSERT(s.replace('G','a')=="Tcat Tcat");
}
void AffxStringTest::equalsTest()
{
	Verbose::out(1, "AffxStringTest::equalsTest");
	AffxString s="acGt acGt";
	CPPUNIT_ASSERT(s.equals("acGt acGt"));
	CPPUNIT_ASSERT(!s.equals("acGT acGt"));
	s="   ";
	CPPUNIT_ASSERT(s.equals("   "));
	CPPUNIT_ASSERT(!s.equals("  "));
}
void AffxStringTest::reverseComplementTest()
{
	Verbose::out(1, "AffxStringTest::reverseComplementTest");
	AffxString s="gggtttcccaaa";
	CPPUNIT_ASSERT(s.reverseComplement()=="tttgggaaaccc");
	s.assign("GGGTTTCCCAAA");
	CPPUNIT_ASSERT(s.reverseComplement()=="tttgggaaaccc");
    s.assign("GGGT_TCCCAAA");
	CPPUNIT_ASSERT(s.reverseComplement()=="tttggganaccc");
	s.assign("GGGT1TCCCAAA");
	CPPUNIT_ASSERT(s.reverseComplement()=="tttggganaccc");
	
}
void AffxStringTest::reverseTest()
{
	Verbose::out(1, "AffxStringTest::reverseTest");
	AffxString s="gggtttcccaaa";
	CPPUNIT_ASSERT(s.reverse()=="aaaccctttggg");
	s="GGGTTTCCCAAA";
	CPPUNIT_ASSERT(s.reverse()=="AAACCCTTTGGG");

}
void AffxStringTest::complementTest()
{
	Verbose::out(1, "AffxStringTest::complementTest");
	AffxString s="gggtttcccaaa";
	CPPUNIT_ASSERT(s.complement()=="cccaaagggttt");
	s="GGGTTTCCCAAA";
	CPPUNIT_ASSERT(s.complement()=="cccaaagggttt");
	
}
void AffxStringTest::convertAmbiguitiesTest()
{
	Verbose::out(1, "AffxStringTest::convertAmbiguitiesTest");
	AffxString s="nugtdtcwcayy";
	s.convertAmbiguities('t');
	CPPUNIT_ASSERT(s=="ttgtttctcatt");
	s="";
	s.convertAmbiguities('t');
	CPPUNIT_ASSERT(s=="");
	s="acgttgca";
	s.convertAmbiguities('t');
	CPPUNIT_ASSERT(s=="acgttgca");
}
void AffxStringTest::startsWithTest()
{
	Verbose::out(1, "AffxStringTest::startsWithTest");
	AffxString s="aCgTA";
	CPPUNIT_ASSERT(s.startsWith("a"));
	CPPUNIT_ASSERT(s.startsWith("aCgTA"));
    CPPUNIT_ASSERT(!s.startsWith("aCgTAT"));
	CPPUNIT_ASSERT(s.startsWith(""));
	s="";
	CPPUNIT_ASSERT(s.startsWith(""));
	CPPUNIT_ASSERT(!s.startsWith("a"));
}
void AffxStringTest::compareToTest()
{
	Verbose::out(1, "AffxStringTest::compareToTest");
	AffxString s="aCgTA";
	CPPUNIT_ASSERT(s.compareTo(AffxString("aCgTA"),0)==0);
	CPPUNIT_ASSERT(s.compareTo(AffxString("aCgTA"),1)==0);
	POSITIVE_TEST(s.compareTo(AffxString("aCgTA"),1));
	CPPUNIT_ASSERT(s.compareTo(AffxString("AaCgTA"),0)==1);
	CPPUNIT_ASSERT(s.compareTo(AffxString("caCgTA"),0)==-1);
}
void AffxStringTest::endsWithTest()
{
	Verbose::out(1, "AffxStringTest::endsWithTest");
	AffxString s="aCgTA";
	CPPUNIT_ASSERT(s.endsWith("TA"));
	CPPUNIT_ASSERT(s.endsWith("aCgTA"));
	CPPUNIT_ASSERT(!s.endsWith("ACgTA"));
	s="";
	CPPUNIT_ASSERT(s.endsWith(""));
	CPPUNIT_ASSERT(!s.endsWith(" "));
}
void AffxStringTest::indexOfTest()
{
	Verbose::out(1, "AffxStringTest::indexOfTest");
	AffxString s="aCgTA";
	CPPUNIT_ASSERT(s.indexOf("TA")==3);
	CPPUNIT_ASSERT(s.indexOf("aCgTA")==0);
	CPPUNIT_ASSERT(s.indexOf("ACgTA")==-1);
	s=" ";
	CPPUNIT_ASSERT(s.indexOf("")==0);
	CPPUNIT_ASSERT(s.indexOf(" ")==0);
}
void AffxStringTest::nextIndexOfTest()
{
	Verbose::out(1, "AffxStringTest::nextIndexOfTest");
	AffxString s="TACaTA";
	CPPUNIT_ASSERT(s.nextIndexOf("TA",1)==4);
	CPPUNIT_ASSERT(s.nextIndexOf("a",3)==-1);
	AffxString s1="hhhhhhhacgthhhhhhhhhacgt";
	CPPUNIT_ASSERT(s1.nextIndexOf("acgt",12)==20);
}
void AffxStringTest::lastIndexOfTest()
{
	Verbose::out(1, "AffxStringTest::lastIndexOfTest");
	AffxString s="TACaTA";
	CPPUNIT_ASSERT(s.lastIndexOf('A')==5);
	CPPUNIT_ASSERT(s.lastIndexOf('a')==3);
	CPPUNIT_ASSERT(s.lastIndexOf('g')==-1);
}
void AffxStringTest::operatorPlusTest()
{
	Verbose::out(1, "AffxStringTest::operatorPlusTest");
	AffxString s="TACaTA";
	CPPUNIT_ASSERT(s+'A'=="TACaTAA");
	CPPUNIT_ASSERT(s+' '=="TACaTAA ");
	CPPUNIT_ASSERT(s+'t'=="TACaTAA t");
	CPPUNIT_ASSERT(s+2=="TACaTAA t2");
	CPPUNIT_ASSERT(s+234567892=="TACaTAA t2234567892");
}
void AffxStringTest::removeSurroundingQuotesTest()
{
	Verbose::out(1, "AffxStringTest::removeSurroundingQuotesTest");
	AffxString s="\"TACaTA\"";
	s.removeSurroundingQuotes();
	CPPUNIT_ASSERT(s=="TACaTA");
	s="'acgt'";
	s.removeSurroundingQuotes();
	CPPUNIT_ASSERT(s=="'acgt'");
	s="TACaTA\"";
	s.removeSurroundingQuotes();
	CPPUNIT_ASSERT(s=="TACaTA\"");
	s="\"TACaTA";
	s.removeSurroundingQuotes();
	CPPUNIT_ASSERT(s=="\"TACaTA");
	s="TACaTA";
	s.removeSurroundingQuotes();
	CPPUNIT_ASSERT(s=="TACaTA");
}
void AffxStringTest::isBaseAtTest()
{
    Verbose::out(1, "AffxStringTest::isBaseAtTest");
	AffxString s="acgt T\tC&GA";
	CPPUNIT_ASSERT(s.isBaseAt(0));
	CPPUNIT_ASSERT(s.isBaseAt(1));
	CPPUNIT_ASSERT(s.isBaseAt(2));
	CPPUNIT_ASSERT(s.isBaseAt(3));
	CPPUNIT_ASSERT(!s.isBaseAt(4));
	CPPUNIT_ASSERT(s.isBaseAt(5));
	CPPUNIT_ASSERT(!s.isBaseAt(6));
	CPPUNIT_ASSERT(s.isBaseAt(7));
	CPPUNIT_ASSERT(!s.isBaseAt(8));
	CPPUNIT_ASSERT(s.isBaseAt(9));
	NEGATIVE_TEST(s.isBaseAt(11),Except);//Out of Bounds exception.
}

void AffxStringTest::isSnpAtTest()
{
    Verbose::out(1, "AffxStringTest::isSnpAtTest");
	AffxString s="rRBd Y\tw&hn";
	CPPUNIT_ASSERT(s.isSnpAt(0));
	CPPUNIT_ASSERT(s.isSnpAt(1));
	CPPUNIT_ASSERT(s.isSnpAt(2));
	CPPUNIT_ASSERT(s.isSnpAt(3));
	CPPUNIT_ASSERT(!s.isSnpAt(4));
	CPPUNIT_ASSERT(s.isSnpAt(5));
	CPPUNIT_ASSERT(!s.isSnpAt(6));
	CPPUNIT_ASSERT(s.isSnpAt(7));
	CPPUNIT_ASSERT(!s.isSnpAt(8));
	CPPUNIT_ASSERT(s.isSnpAt(9));
	NEGATIVE_TEST(s.isSnpAt(11),Except);//Out of Bounds exception.
}

void AffxStringTest::isIUPACCodeAtTest()
{
    Verbose::out(1, "AffxStringTest::isIUPACCodeAtTest");
	AffxString s="rRBd Y\tw&hn";
	CPPUNIT_ASSERT(s.isIUPACCodeAt(0));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(1));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(2));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(3));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(4));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(5));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(6));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(7));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(8));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(9));
	NEGATIVE_TEST(s.isIUPACCodeAt(11),Except);//Out of Bounds exception.
	s="acgt T\tC&GA";
	CPPUNIT_ASSERT(s.isIUPACCodeAt(0));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(1));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(2));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(3));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(4));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(5));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(6));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(7));
	CPPUNIT_ASSERT(!s.isIUPACCodeAt(8));
	CPPUNIT_ASSERT(s.isIUPACCodeAt(9));
	NEGATIVE_TEST(s.isIUPACCodeAt(11),Except);//Out of Bounds exception.
}

void AffxStringTest::getSnpStringAtTest()
{
    Verbose::out(1, "AffxStringTest::getSnpStringAtTest");
	AffxString s="aAcCgGtTrRyYmMkKwWsSbBdDhHvVnNF";
	CPPUNIT_ASSERT(s.getSnpStringAt(0)=="A");
	CPPUNIT_ASSERT(s.getSnpStringAt(1)=="A");
	CPPUNIT_ASSERT(s.getSnpStringAt(2)=="C");
	CPPUNIT_ASSERT(s.getSnpStringAt(3)=="C");
	CPPUNIT_ASSERT(s.getSnpStringAt(4)=="G");
	CPPUNIT_ASSERT(s.getSnpStringAt(5)=="G");
	CPPUNIT_ASSERT(s.getSnpStringAt(6)=="T");
	CPPUNIT_ASSERT(s.getSnpStringAt(7)=="T");
	CPPUNIT_ASSERT(s.getSnpStringAt(8)=="AG");
	CPPUNIT_ASSERT(s.getSnpStringAt(9)=="AG");
	CPPUNIT_ASSERT(s.getSnpStringAt(10)=="CT");
	CPPUNIT_ASSERT(s.getSnpStringAt(11)=="CT");
	CPPUNIT_ASSERT(s.getSnpStringAt(12)=="AC");
	CPPUNIT_ASSERT(s.getSnpStringAt(13)=="AC");
	CPPUNIT_ASSERT(s.getSnpStringAt(14)=="GT");
	CPPUNIT_ASSERT(s.getSnpStringAt(15)=="GT");
	CPPUNIT_ASSERT(s.getSnpStringAt(16)=="AT");
	CPPUNIT_ASSERT(s.getSnpStringAt(17)=="AT");
	CPPUNIT_ASSERT(s.getSnpStringAt(18)=="CG");
	CPPUNIT_ASSERT(s.getSnpStringAt(19)=="CG");
	CPPUNIT_ASSERT(s.getSnpStringAt(20)=="CGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(21)=="CGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(22)=="AGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(23)=="AGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(24)=="ACT");
	CPPUNIT_ASSERT(s.getSnpStringAt(25)=="ACT");
	CPPUNIT_ASSERT(s.getSnpStringAt(26)=="ACG");
	CPPUNIT_ASSERT(s.getSnpStringAt(27)=="ACG");
	CPPUNIT_ASSERT(s.getSnpStringAt(28)=="ACGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(29)=="ACGT");
	CPPUNIT_ASSERT(s.getSnpStringAt(30)=="");
	NEGATIVE_TEST(s.getSnpStringAt(31),Except); //Out of Bounds exception
}

void AffxStringTest::hashTest()
{
    Verbose::out(1, "AffxStringTest::hashTest");
	AffxString s="";
	CPPUNIT_ASSERT(s.hash()==0);
    s="a";
	CPPUNIT_ASSERT(s.hash()==97);
	s="ac";
	CPPUNIT_ASSERT(s.hash()==3300);
	s="acG";
	CPPUNIT_ASSERT(s.hash()==108971);
	s="acGT";
	CPPUNIT_ASSERT(s.hash()==3596127);
}

void AffxStringTest::constructorTest()
{
    Verbose::out(1, "AffxStringTest::constructorTest");
	//AffxString::AffxString(const std::string& stringSrc)
	string s="acgt";
    AffxString test= string("acgt");
    CPPUNIT_ASSERT(test==s);
	AffxString t1,t2,t3,t4;
	//AffxString::AffxString(char ch, int nRepeat)
	t1=AffxString('c',5);
	CPPUNIT_ASSERT(t1=="ccccc");
	//AffxString::AffxString(const char* pch, int nLength)
	t2= AffxString("acgtagcgtacgt",8);
    CPPUNIT_ASSERT(t2=="acgtagcg");
	//AffxString::AffxString(const AffxString& stringSrc)
    t3=AffxString(t2);
	CPPUNIT_ASSERT(t2==t3);
	//AffxString::AffxString(const char* psz)
    t4=AffxString("ACGTACGT");
	CPPUNIT_ASSERT(t4=="ACGTACGT");
}
