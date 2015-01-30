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

/**
 * @file   UtilTest.cpp
 * @author Chuck Sugnet
 * @date   Tue Mar 21 10:15:03 PST 2006
 *
 * @brief  Testing the util functions.
 *
 */

//

#include "util/AffxString.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
// Setup.h defines AffxString? therefore AffxString.h needs to come first for windows. 
#include "util/CPPTest/Setup.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
//
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#include <sys/file.h>  // for flock() function in Linux
#endif
using namespace std;

/**
 * @class UtilTest
 * @brief cppunit class for testing conversion functions.
 */
class UtilTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( UtilTest );
  CPPUNIT_TEST(testEscapeString);
  CPPUNIT_TEST(testDeEscapeString);
  CPPUNIT_TEST(schrageRandomTest );
  CPPUNIT_TEST(stringEndTest );
  CPPUNIT_TEST(stringChangeEndTest);
  //vliber
  CPPUNIT_TEST(compareStringTest);
  CPPUNIT_TEST(downcaseStringTest);
  CPPUNIT_TEST(roundTest);
  CPPUNIT_TEST(cloneStringTest);
  CPPUNIT_TEST(chompLastIfSepTest);
  CPPUNIT_TEST(chopSuffix);
  CPPUNIT_TEST(sameStringTest);
  CPPUNIT_TEST(chopStringTest);
  CPPUNIT_TEST(replaceStringTest);
  CPPUNIT_TEST(printStringWidthTest);
  CPPUNIT_TEST(nextWhiteSpaceTest);
  CPPUNIT_TEST(memInfo_win32Test);
  CPPUNIT_TEST(matrixDifferencesTest);
  // rsatin
  CPPUNIT_TEST(changeEndTest);
  CPPUNIT_TEST(memInfoTest);
  CPPUNIT_TEST(endsWithStrTest);
  CPPUNIT_TEST_SUITE_END();

public:
  void testEscapeString();
  void testDeEscapeString();
  /** Test the random number generator. */
  void schrageRandomTest();
  void addPrefixSuffixTest();
  void joinVectorStringTest();
  void stringEndTest();
  void stringChangeEndTest();
  //vliber
  void compareStringTest();
  void downcaseStringTest();
  void roundTest();
  void cloneStringTest();
  void chompLastIfSepTest();
  void chopSuffix();
  //void createDirTest();
  void sameStringTest();
  void chopStringTest();
    void replaceStringTest();
  void printStringWidthTest();
  void nextWhiteSpaceTest();
  void memInfo_win32Test();
  void matrixDifferencesTest();
  //rsatin
  void changeEndTest();
  void memInfoTest();
  
  void endsWithStrTest();

  void setUp() {
    if ( !Fs::dirExists(OUTPUT) ) {
      Fs::mkdirPath(OUTPUT);
    }
  };
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( UtilTest );

void UtilTest::schrageRandomTest() {
	Verbose::out(1, "");
	Verbose::out(1, "***UtilTest testcases***");
  Verbose::out(1, "UtilTest::schrageRandomTest");
  int32_t test = 42;
  int32_t gold[] = {705894,1126542223,1579310009,565444343,807934826,421520601,2095673201,1100194760,1139130650,552121545,229968128,1751246343,1933904666,970724117,526968160,531304892,404315618,694332618,222172928,1733822410,1147638727,1813450982,1582736250,168218361,1157513875,281338952,1852259217,997713207,1013554273,966378307,501383488,26451988,49447387,2133545567,1965890610,1687573175,1225826296,1651931201,1339106791,739215777,806666144,573618697,755349096,1376419055,781211901,113402349,1135284754,338656883,974568031,707121348};
  bool good = true;
  for(int i = 0; i < ArraySize(gold); i++) {
    int32_t result = Util::schrageRandom(&test);
    good &= (gold[i] == result);
  }
  CPPUNIT_ASSERT(good);
}

void UtilTest::testEscapeString() {
  string start = ",hello,here,we,are,";
  string gold = "\\,hello\\,here\\,we\\,are\\,";
  string result = Util::escapeString(start, ',');
  CPPUNIT_ASSERT(result == gold);
  start = ",hel\\lo,here,we,are,";
  gold = "\\,hel\\\\lo\\,here\\,we\\,are\\,";
  result = Util::escapeString(start, ',');
  CPPUNIT_ASSERT(result == gold);
}

void UtilTest::testDeEscapeString() {
  string start = ",hello,here,we,are,";
  string gold = "\\,hello\\,here\\,we\\,are\\,";
  string result = Util::escapeString(start, ',');
  string reversed = Util::deEscapeString(result);
  CPPUNIT_ASSERT(result == gold);
  CPPUNIT_ASSERT(reversed == start);
  start = ",hel\\lo,here,we,are,";
  gold = "\\,hel\\\\lo\\,here\\,we\\,are\\,";
  result = Util::escapeString(start, ',');
  reversed = Util::deEscapeString(result);
  CPPUNIT_ASSERT(reversed == start);
  CPPUNIT_ASSERT(result == gold);
}

void UtilTest::joinVectorStringTest()
{
  std::vector<std::string> myVec;
  myVec.push_back("One");
  myVec.push_back("Two");
  myVec.push_back("Three");
  myVec.push_back("Four");

  std::string myString = ", ";

  std::string test = Util::joinVectorString(myVec, myString);
  CPPUNIT_ASSERT(test == "One, Two, Three, Four, ");
}

void UtilTest::addPrefixSuffixTest()
{
  const char* myChar[] = {"AAAAAAAA", "BBBBBBBB", "CCCCCCCC", "DDDDDDDD", NULL}; 
  int csize = 4;

  string myString[] = {"One", "Two", "Three", "Four", "Five", ""};
  int ssize = 5;
  
  std::string prefix = "./testthis/";
  std::string suffix = ".cel";
    
  std::vector<std::string> myVec;
  myVec.push_back("First");
  myVec.push_back("Second");
  myVec.push_back("Third");
  myVec.push_back("Fourth");
  myVec.push_back("Fifth");
  
  std::vector<std::string> out;
  //Test Char* array input
  out=Util::addPrefixSuffix(myChar, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/AAAAAAAA.cel");
  CPPUNIT_ASSERT(out[1]=="./testthis/BBBBBBBB.cel");
  CPPUNIT_ASSERT(out[2]=="./testthis/CCCCCCCC.cel");
  CPPUNIT_ASSERT(out[3]=="./testthis/DDDDDDDD.cel");

  //Test Char* array input with size of array
  out=Util::addPrefixSuffix(myChar,csize, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/AAAAAAAA.cel");
  CPPUNIT_ASSERT(out[1]=="./testthis/BBBBBBBB.cel");
  CPPUNIT_ASSERT(out[2]=="./testthis/CCCCCCCC.cel");
  CPPUNIT_ASSERT(out[3]=="./testthis/DDDDDDDD.cel");

  //Test Char* array input
  out=Util::addPrefixSuffix(myChar, prefix);
  CPPUNIT_ASSERT(out[0]=="./testthis/AAAAAAAA");
  CPPUNIT_ASSERT(out[1]=="./testthis/BBBBBBBB");
  CPPUNIT_ASSERT(out[2]=="./testthis/CCCCCCCC");
  CPPUNIT_ASSERT(out[3]=="./testthis/DDDDDDDD");

  //Test String array input
  out=Util::addPrefixSuffix(myString, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/One.cel");
  CPPUNIT_ASSERT(out[1]=="./testthis/Two.cel");
  CPPUNIT_ASSERT(out[2]=="./testthis/Three.cel");
  CPPUNIT_ASSERT(out[3]=="./testthis/Four.cel");
  CPPUNIT_ASSERT(out[4]=="./testthis/Five.cel");  

  //Test String array input with size of array
  out=Util::addPrefixSuffix(myString, ssize, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/One.cel");
  CPPUNIT_ASSERT(out[1]=="./testthis/Two.cel");
  CPPUNIT_ASSERT(out[2]=="./testthis/Three.cel");
  CPPUNIT_ASSERT(out[3]=="./testthis/Four.cel");
  CPPUNIT_ASSERT(out[4]=="./testthis/Five.cel");  

  //Test string vector input
  out = Util::addPrefixSuffix(myVec, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/First.cel");
  CPPUNIT_ASSERT(out[1]=="./testthis/Second.cel");
  CPPUNIT_ASSERT(out[2]=="./testthis/Third.cel");
  CPPUNIT_ASSERT(out[3]=="./testthis/Fourth.cel");
  CPPUNIT_ASSERT(out[4]=="./testthis/Fifth.cel");  

  const char* myCharNull[]={NULL};
  csize = 0;
  out=Util::addPrefixSuffix(myCharNull, prefix, suffix);
  CPPUNIT_ASSERT(out.size()==0);

  out=Util::addPrefixSuffix(myCharNull, csize, prefix, suffix);
  CPPUNIT_ASSERT(out.size()==0);

  string myStringNull[] = {""};
  ssize = 0;
  out=Util::addPrefixSuffix(myStringNull, prefix, suffix);
  CPPUNIT_ASSERT(out.size()==0);

  out=Util::addPrefixSuffix(myStringNull, csize, prefix, suffix);
  CPPUNIT_ASSERT(out.size()==0);

  const char* myCharSmall[]={"solo", NULL};
  csize = 1;
  out=Util::addPrefixSuffix(myCharSmall, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/solo.cel");

  out=Util::addPrefixSuffix(myCharSmall, csize, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/solo.cel");

  string myStringSmall[]={"yay", ""};
  ssize = 1;
  out=Util::addPrefixSuffix(myStringSmall, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/yay.cel");

  out=Util::addPrefixSuffix(myStringSmall, ssize, prefix, suffix);
  CPPUNIT_ASSERT(out[0]=="./testthis/yay.cel");
}

void UtilTest::stringEndTest()
{
  Verbose::out(1, "UtilTest::stringEndTest");
  CPPUNIT_ASSERT(Util::stringEndsWith(""    ,""     )==true );
  CPPUNIT_ASSERT(Util::stringEndsWith(""    ,"a"    )==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("a"   ,"aa"   )==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("aa"  ,"aa"   )==true );
  CPPUNIT_ASSERT(Util::stringEndsWith("aaaa",""     )==true );
  CPPUNIT_ASSERT(Util::stringEndsWith("abb" ,"a"    )==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("abba","a"    )==true );
  CPPUNIT_ASSERT(Util::stringEndsWith("abba","b"    )==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("abba","c"    )==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("f.txt",".txt")==true );
  CPPUNIT_ASSERT(Util::stringEndsWith("f.bin",".txt")==false);
  //additional tests vliber
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_s","_s")==true);
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_a_at","_at")==true);
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_a_at","a_at")==true);
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_a_at","s_at")==false);
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_x_st","_st")==true);
  CPPUNIT_ASSERT(Util::stringEndsWith("abcd_x_st","_at")==false);
}

void UtilTest::stringChangeEndTest()
{
	Verbose::out(1, "UtilTest::stringChangeEndTest");
  std::string str;
  std::vector<std::string> str_vec;

  //
  str="str.aaa";
  Util::changeEnd(str,".aaa",".bbb");
  CPPUNIT_ASSERT(str=="str.bbb");
  Util::changeEnd(str,".ccc",".ddd");
  CPPUNIT_ASSERT(str=="str.bbb");

  //
  str="str.aaa.aaa";
  Util::changeEnd(str,".aaa",".bbb");
  CPPUNIT_ASSERT(str=="str.aaa.bbb");
  Util::changeEnd(str,".ccc",".ddd");
  CPPUNIT_ASSERT(str=="str.aaa.bbb");

  //
  str_vec.clear();
  str_vec.push_back("str.aaa");
  //
  Util::changeEnd(str_vec,".aaa",".bbb");
  CPPUNIT_ASSERT(str_vec[0]=="str.bbb");
  Util::changeEnd(str_vec,".ccc",".ddd");
  CPPUNIT_ASSERT(str_vec[0]=="str.bbb");
}
//test downcaseString() vliber
void UtilTest::compareStringTest()
{
	Util::ltstring lt;
	Verbose::out(1, "UtilTest::compareStringTest");
    CPPUNIT_ASSERT(lt("ABCDTH","aBCDTH"));//<
	CPPUNIT_ASSERT(!lt("ABCDTH","ABCDTH"));//>
	CPPUNIT_ASSERT(!lt("abfdh","abfdh"));//==
	CPPUNIT_ASSERT(!lt("aBCDTH","ABCDTH"));//>

	char a1[]="ABCDTH";
    char a2[]="AbCDTH";
	CPPUNIT_ASSERT(lt(a1,a2));//<
	char a3[]="AGHDTH";
    char a4[]="AGHDTH";
	CPPUNIT_ASSERT(!lt(a3,a4));//==
	char a5[]="AbHDTH";
    char a6[]="ABHDTH";
	CPPUNIT_ASSERT(!lt(a5,a6));//>

	string s1 ("ABCJKU");
	string s2 ("AbCJKU");
	CPPUNIT_ASSERT(lt(s1,s2));//<
	string s3 ("pou");
	string s4 ("pou");
	CPPUNIT_ASSERT(!lt(s3,s4));//==
	string s5 ("AbCJKU");
	string s6 ("ABCJKU");
	CPPUNIT_ASSERT(!lt(s5,s6));//>
} 
//test downcaseString() vliber
void UtilTest::downcaseStringTest()
{
	Verbose::out(1, "UtilTest::downcaseStringTest");
	string s1 ("ABCDTH");
	Util::downcaseString_inplace(s1);
	CPPUNIT_ASSERT_EQUAL(s1,string("abcdth"));
	string s2 ("");
	Util::downcaseString_inplace(s2);
	CPPUNIT_ASSERT_EQUAL(s2,std::string(""));
	string s3 ("abcd");
	Util::downcaseString_inplace(s3);
	CPPUNIT_ASSERT_EQUAL(s3,std::string("abcd"));
	
    string s4 ("3456^&_*%");
	Util::downcaseString_inplace(s4);
	CPPUNIT_ASSERT_EQUAL(s4,std::string("3456^&_*%"));
} 
//test round() vliber
void UtilTest::roundTest()
{
	Verbose::out(1, "UtilTest::roundTest");
	CPPUNIT_ASSERT_EQUAL(Util::round(4.0),4);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.01),4);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.49),4);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.50),5);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.51),5);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.99),5);
	CPPUNIT_ASSERT_EQUAL(Util::round(4.50),5);
/* FIXME: rsatin 26apr09 bug causes all following unit tests to fail, commented out until fixed
	double const NAN = std::numeric_limits<double>::quiet_NaN(); 
	double const INFINITY = std::numeric_limits<double>::infinity();
	CPPUNIT_ASSERT( Util::round( NAN ) == NAN );                     // NAN converts to extreme negative value INT_MIN
	CPPUNIT_ASSERT( Util::round( -INFINITY ) == -INFINITY );         // negative Infinity saturates output
	CPPUNIT_ASSERT( Util::round(  INFINITY ) ==  INFINITY );         // positive Infinity saturates output
	CPPUNIT_ASSERT_EQUAL(0.0,(double)Util::round((double)INT_MIN-2.0)-((double)INT_MIN-2.0));  // extreme negative values saturate
	CPPUNIT_ASSERT_EQUAL(0.0,(double)Util::round((double)INT_MAX+2.0)-((double)INT_MAX+2.0));  // extreme positive values overflow to extreme negative values
*/
} 
//test cloneString() vliber
void UtilTest::cloneStringTest()
{
  Verbose::out(1, "UtilTest::cloneStringTest");
  char input1[]="{[!@#$%^&-1234567890";
  char *input2=Util::cloneString(input1);
  CPPUNIT_ASSERT_EQUAL(string(input1),string(input2));
  int i=0;
  while(input1[i] != 0)
  {
	  CPPUNIT_ASSERT_EQUAL(input1[i], input2[i]);
	  i++;
  } 
  char input3[]="";
  char *input4=Util::cloneString(input3);
  CPPUNIT_ASSERT_EQUAL(*input3,*input4);
  //string
  string input1a="{[!@#$%^&-1234567890";
  string input2a=Util::cloneString(input1a.c_str());
  CPPUNIT_ASSERT(input1a==input2a);
  //literal
  CPPUNIT_ASSERT(input1a==Util::cloneString("{[!@#$%^&-1234567890"));
  
}
//test chompLastIfSep() /
void UtilTest::chompLastIfSepTest()
{
	Verbose::out(1, "UtilTest::chompLastIfSepTest");
	#ifdef WIN32
	  string s1 (".\\input\\");
	  Util::chompLastIfSep(s1);
	  CPPUNIT_ASSERT_EQUAL(s1,string(".\\input"));
	  string s2 ("./input/");
	  Util::chompLastIfSep(s2);
	  CPPUNIT_ASSERT_EQUAL(s2,string("./input/")); //todo vliber support only Windows back slash. What about forward slash
    #else
	  string s2 ("./input/");
	  Util::chompLastIfSep(s2);
	  CPPUNIT_ASSERT_EQUAL(s2,string("./input"));  //todo rsating should be same result on Windows for slash direction independence
    #endif
}

//test chompLastIfSep() /
void UtilTest::chopSuffix()
{
	Verbose::out(1, "UtilTest::chopSuffix");
        string chop("mycel.cel");
        string comp = Util::chopSuffix(chop);
        CPPUNIT_ASSERT_EQUAL(comp,string("mycel"));

        string chop2("mycel");
        string comp2 = Util::chopSuffix(chop2);
        CPPUNIT_ASSERT_EQUAL(comp2,string("mycel"));

        string chop3("my.cel.cel");
        string comp3 = Util::chopSuffix(chop3);
        CPPUNIT_ASSERT_EQUAL(comp3,string("my.cel"));

}

//test directoryReadable() vliber


//test sameString() vliber
void UtilTest::sameStringTest()
{
	Verbose::out(1, "UtilTest::sameStringTest");
  CPPUNIT_ASSERT(Util::sameString("ACGT","ACGT"));
  CPPUNIT_ASSERT(!Util::sameString("aCGT","ACGT"));
  char a[]="abcd";
  char b[]="abcd";
  char c[]="Abcd";
  CPPUNIT_ASSERT(Util::sameString(a,b));
  CPPUNIT_ASSERT(!Util::sameString(a,c));
  string a_s="abcd";
  string b_s="abcd";
  string c_s="Abcd";
  CPPUNIT_ASSERT(Util::sameString(a_s,b_s));
  CPPUNIT_ASSERT(!Util::sameString(a_s,c_s));
}


void UtilTest::replaceStringTest() {
    string s = "herex x we go";
    string from = "x ";
    string to = "yy";
    Util::replaceString(s, from, to);
    CPPUNIT_ASSERT(s == "hereyyyywe go");
}

//test chopString() vliber
void UtilTest::chopStringTest()
{
	Verbose::out(1, "UtilTest::chopStringTest");

	//positive
	string sInput1("aaaa/cccc/gggg/tttt/");
	char del ='/';
	std::vector<string> v;
	std::vector<string>::const_iterator i;
	Util::chopString(sInput1,del,v);
	string sOutput1;
	for(i=v.begin();i!=v.end();i++)
	{
		sOutput1=sOutput1+(*i)+'/';
	}
	v.clear();
	CPPUNIT_ASSERT_EQUAL(sInput1,sOutput1);
    string sInput2("aaaa\\cccc\\gggg\\tttt\\");
	del ='\\';
	Util::chopString(sInput2,del,v);
	string sOutput2;
	for(i=v.begin();i!=v.end();i++)
	{
	  sOutput2=sOutput2+(*i)+'\\';
	}
	CPPUNIT_ASSERT_EQUAL(sInput2,sOutput2);
	v.clear();

    //chopString==whole input string
	del =' ';
	POSITIVE_TEST(Util::chopString(sInput2,del,v));
	CPPUNIT_ASSERT_EQUAL(sInput2,v[0]);
    v.clear();

	string sInput3("");
	del ='/';
	POSITIVE_TEST(Util::chopString(sInput3,del,v));
	
}

//test printStringWidth() vliber
void UtilTest::printStringWidthTest()
{
	 Verbose::out(1, "UtilTest::printStringWidthTest");
   filebuf fb;
   AffxString s1=OUTPUT+"/test.txt";
   fb.open (s1.c_str(),ios::out);
   ostream output (&fb) ;
   string input ("aaaaattttt bbbbbkkkkk cccccjjjjj");
   Util::printStringWidth(output,input,0,0,40);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,0,10,40);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,0,20,40);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,0,30,40);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,0,0,0);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,5,0,0);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,5,0,40);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,5,0,30);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,5,0,20);
   output<<"\n"<<endl;
   Util::printStringWidth(output,input,5,0,10);
   output<<"\n"<<endl;
   string input1 ("kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk");
   Util::printStringWidth(output,input1,0,0,0);
   fb.close();
   //compare files
   AffxString s2=OUTPUT+"/test.txt";
   ifstream ifs1(s2.c_str());
   AffxString s3=INPUT+"/testGold.txt";
   ifstream ifs2(s3.c_str());
   string line1;
   string line2;
   while (!ifs2.eof() || !ifs1.eof())
   {
	  getline(ifs2,line2);
	  getline(ifs1,line1);
	  // drop trailing CR to ignore DOS vs Unix text file difference
	  int rf1 = line1.rfind("\r");
	  int rf2 = line2.rfind("\r");
	  if( rf1>=0 )
	    line1.resize(rf1);
	  if( rf2>=0 )
	    line2.resize(rf2);
	  CPPUNIT_ASSERT(line1==line2);
    }
	ifs1.close();
	ifs2.close();
} 

//test nextWhiteSpace() vliber
void UtilTest::nextWhiteSpaceTest()
{
	//char
   Verbose::out(1, "UtilTest::nextWhiteSpaceTest");
   const char *input1="aaaaa bbbbb ccccc ttttt";
   string s1 (Util::nextWhiteSpace(input1));
   string s2 (" bbbbb ccccc ttttt");
   CPPUNIT_ASSERT(s1==s2);
   
   char input2[]="aaaaabbbbbcccccttttt";
   s1=(Util::nextWhiteSpace(input2));
   s2=("");
   CPPUNIT_ASSERT(s1==s2);

   char input3[]="";
   s1=(Util::nextWhiteSpace(input3));
   s2=("");
   CPPUNIT_ASSERT(s1==s2);

   char input4[]=" aaaaabbbbbcccccttttt";
   s1=(Util::nextWhiteSpace(input4));
   s2=(" aaaaabbbbbcccccttttt");
   CPPUNIT_ASSERT(s1==s2);

   char input5[]="aaaaabbbbbcccccttttt ";
   s1=(Util::nextWhiteSpace(input5));
   s2=(" ");
   CPPUNIT_ASSERT(s1==s2);

   //string
   string input1_s="aaaaa bbbbb ccccc ttttt";
   string s1_s (Util::nextWhiteSpace(input1_s.c_str()));
   string s2_s (" bbbbb ccccc ttttt");
   CPPUNIT_ASSERT(s1_s==s2_s);
   
   string input2_s="aaaaabbbbbcccccttttt";
   s1_s=(Util::nextWhiteSpace(input2_s.c_str()));
   s2_s=("");
   CPPUNIT_ASSERT(s1_s==s2_s);

   string input3_s="";
   s1_s=(Util::nextWhiteSpace(input3_s.c_str()));
   s2_s=("");
   CPPUNIT_ASSERT(s1_s==s2_s);

   string input4_s=" aaaaabbbbbcccccttttt";
   s1_s=(Util::nextWhiteSpace(input4_s.c_str()));
   s2_s=(" aaaaabbbbbcccccttttt");
   CPPUNIT_ASSERT(s1_s==s2_s);

   string input5_s="aaaaabbbbbcccccttttt ";
   s1_s=(Util::nextWhiteSpace(input5_s.c_str()));
   s2_s=(" ");
   CPPUNIT_ASSERT(s1_s==s2_s);
   
   //literal
   string s1_l (Util::nextWhiteSpace("aaaaa bbbbb ccccc ttttt"));
   string s2_l (" bbbbb ccccc ttttt");
   CPPUNIT_ASSERT(s1_l==s2_l);

   s1_l=(Util::nextWhiteSpace("aaaaabbbbbcccccttttt "));
   s2_l= (" ");
   CPPUNIT_ASSERT(s1_l==s2_l);
}



//test memInfo_win32() vliber
void UtilTest::memInfo_win32Test()
{
	Verbose::out(1, "UtilTest::memInfo_win32Test");
   uint64_t free, total, swapAvail, memAvail;
   CPPUNIT_ASSERT(Util::memInfo(free, total, swapAvail, memAvail)==true);
}
//test matrixDifferences() vliber
void UtilTest::matrixDifferencesTest()
{
	Verbose::out(1, "UtilTest::matrixDifferencesTest");
   //negative empty paths
   //Goal is to receive a message: FATAL ERROR: Can't open file  to read. -OK
   NEGATIVE_TEST(Util::matrixDifferences("","",0,0,0.0,false,false),std::exception);
   NEGATIVE_TEST(Util::matrixDifferences("","",0,0,0.0,true,true),std::exception);
   NEGATIVE_TEST(Util::matrixDifferences("","",0,0,0.0,false,true),std::exception);
   NEGATIVE_TEST(Util::matrixDifferences("","",0,0,0.0,true,false),std::exception);
   //negative not double
   //Goal is to receive a message: FATAL ERROR: Could not convert 'X' to a double. -OK
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test_1lqGold.txt",INPUT+"/test_1lq.txt",0,0,0.0,true,true),std::exception);
   //Goal is to receive a message: FATAL ERROR: Could not convert 'AGGACTTGCCATCGTAGAACTGCTG' to a double. -OK
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test_1lqGold.txt",INPUT+"/test_1lq.txt",0,1,0.0,true,true),std::exception);
   //entry=null
   //Goal is to receive a message:  FATAL ERROR: Could not convert '' to a double. -OK
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test7Bad.txt",0,1,0.0,true,true),std::exception);
   //negative out of range
   //Goal is to receive a message: FATAL ERROR: RowFile::matrixFromFile() - Number of skipCols >= number of cols
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",7,1,0.0,true,true),std::exception);
   //Goal is to receive a message: FATAL ERROR: Nothing after header in file: .\input\test1Gold.txt
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,11,0.0,true,true),std::exception);
   //negative rowname
   //Goal is to receive a message: FATAL ERROR: Can't find rowname: 33 in matrix 1
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test2Bad.txt",0,1,0.0,true,true),std::exception);
   //Goal is to receive a message: FATAL ERROR: Entry: '33' has already been seen in file: .\input\test3Bad.txt in column 0.
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test3Bad.txt",0,1,0.0,true,true),std::exception);
   //Goal is to receive a message: FATAL ERROR: Entry: '8' has already been seen in file: .\input\test1Gold_Bad3.txt in column 0.
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold_Bad3.txt",INPUT+"/test1.txt",0,1,0.346,true,true),std::exception);
   //Goal is to receive a message: FATAL ERROR: Matrices are different sizes, not comparable. -OK
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test3.txt",0,1,0.0,true,true),std::exception);
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test3.txt",0,1,0.0,true,false),std::exception);
   //Goal is to receive a message: FATAL ERROR: Can't find rowname: 8 in matrix 1. -OK
   NEGATIVE_TEST(Util::matrixDifferences(INPUT+"/test1Gold_Bad2.txt",INPUT+"/test1.txt",0,1,0.0,true,true),std::exception);


   //positive identical
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1.txt",0,1,0.0,false,false)==0);
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1.txt",0,1,0.0,true,true)==0);
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1.txt",0,1,0.0,false,true)==0);
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1.txt",0,1,0.0,true,false)==0);
   //positive diff value; whole file
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,1,0.0,false,false)==4);
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,1,0.0,true,false)==4); //print output
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,1,0.0,false,true)==4);
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,1,0.0,true,true)==4); //print output
   
  
   //positive diff value; increase col & row
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",3,1,0.0,true,true)==3);//print output
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",6,1,0.0,true,true)==1);//print output
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,5,0.0,true,true)==3);//print output
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,8,0.0,true,true)==1);//print output
   //increase epsilon; no diff
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test1Bad.txt",0,1,2.0,true,true)==3); //print output
   //positive diff value; match rows
   //Goal is to receive a message: File: input\test2.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 60.
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test2.txt",0,1,0.0,true,false)==60); //print output
   //no errors in the output because matchRows=true
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test2.txt",0,1,0.0,true,true)==0); 

   //Goal is to receive a message: row: 8 col: 4 (10.3457 vs. 10) Diff: 0.345679 -OK
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test6Bad.txt",0,1,0.0,true,true)==1); //print output
   //increase epsilon; no diff
   CPPUNIT_ASSERT(Util::matrixDifferences(INPUT+"/test1Gold.txt",INPUT+"/test6Bad.txt",0,1,0.346,true,true)==0);
   
}

//test changeEnd() rsatin
void UtilTest::changeEndTest()
{
	Verbose::out(1, "UtilTest::changeEndTest");
	std::string tmpDir = INPUT+"/output_util_FileIO/";                             // temporary directory for testing
	// test with std:string input argument 
	std::string string_instance = tmpDir+"file1.txt";
    POSITIVE_TEST( Util::changeEnd( string_instance,".txt",".csv") );              // test matching pattern changed
	CPPUNIT_ASSERT( string_instance == tmpDir+"file1.csv" );
    POSITIVE_TEST( Util::changeEnd( string_instance,".xxx",".yyy") );              // test non-matching pattern not changed
	CPPUNIT_ASSERT( string_instance == tmpDir+"file1.csv" );
	// test with std::vector input argument 
	std::vector<std::string> string_vector;
	string_vector.insert( string_vector.end(), std::string(tmpDir+"file2.txt") );
	string_vector.insert( string_vector.end(), std::string(tmpDir+"file3.txt") );
	string_vector.insert( string_vector.end(), std::string(tmpDir+"file4.txt") );
    POSITIVE_TEST( Util::changeEnd( string_vector,"2.txt","X.txt") );              // test single matching pattern changed
	CPPUNIT_ASSERT( string_vector[0] == tmpDir+"fileX.txt" );
	CPPUNIT_ASSERT( string_vector[1] == tmpDir+"file3.txt" );
	CPPUNIT_ASSERT( string_vector[2] == tmpDir+"file4.txt" );
    POSITIVE_TEST( Util::changeEnd( string_vector,".txt",".csv") );                // test multiple matching pattern changed
	CPPUNIT_ASSERT( string_vector[0] == tmpDir+"fileX.csv" );
	CPPUNIT_ASSERT( string_vector[1] == tmpDir+"file3.csv" );
	CPPUNIT_ASSERT( string_vector[2] == tmpDir+"file4.csv" );
    POSITIVE_TEST( Util::changeEnd( string_vector,".xxx",".yyy") );                // test multiple matching pattern changed
	CPPUNIT_ASSERT( string_vector[0] == tmpDir+"fileX.csv" );
	CPPUNIT_ASSERT( string_vector[1] == tmpDir+"file3.csv" );
	CPPUNIT_ASSERT( string_vector[2] == tmpDir+"file4.csv" );
}

//test memInfo() rsatin
void UtilTest::memInfoTest()
{
	Verbose::out(1, "UtilTest::memInfoTest");
	uint64_t free;
	uint64_t total;
	uint64_t swapAvail;
	uint64_t memAvail;
	bool cap32bit = true;
	uint64_t oneGb = 1024*1024*1024;
    // BUG rsatin 01may09 ==> swapAvail is many orders of magnitude out of sync with "free", "total", and "memAvail"
	// test with 32-bit Windows OS cap in place (believable limits are rough estimates)
	CPPUNIT_ASSERT( Util::memInfo(free,total,swapAvail,memAvail,cap32bit) );
printf("memInfoTest: %ld,%ld,%ld,%ld\n",(long int)(free/1048576),(long int)(total/1048576),(long int)(swapAvail/1048576),(long int)(memAvail/1048576));
    //This test is fragile on a busy/loaded system
	//CPPUNIT_ASSERT(      free>oneGb/32 &&      free<=64*oneGb );         // believable values between 37Mb and 64 Gb
    CPPUNIT_ASSERT(     total>oneGb   &&     total<=256*oneGb );        // believable values between 1Gb and 256 Gb
    #ifdef NEVER    //todo rsatin - fails on Windows and Darwin ---v
    CPPUNIT_ASSERT( swapAvail>oneGb   && swapAvail<=256*oneGb );        // believable values between 1Gb and 256 Gb
    #endif          //todo rsatin - fails on Windows and Darwin ---^
    #ifdef WIN32    //todo rsatin - fails on Linux ---v and Windows
	  //CPPUNIT_ASSERT(  memAvail>oneGb/8 &&  memAvail<=MEMINFO_2GB_MAX );  // believable values between 128Mb and conservative limit
    #endif          //todo rsatin - fails on Linux ---^
	// test without 32-bit Windows OS cap in place (believable limits are rough estimates)
	cap32bit = false;
	CPPUNIT_ASSERT( Util::memInfo(free,total,swapAvail,memAvail,cap32bit) );
printf("memInfoTest: %ld,%ld,%ld,%ld\n",(long int)(free/1048576),(long int)(total/1048576),(long int)(swapAvail/1048576),(long int)(memAvail/1048576));
    //This test is fragile on a busy/loaded system
	//CPPUNIT_ASSERT(      free>oneGb/8 &&      free<=64*oneGb );         // believable values between 128Mb and 64 Gb
	CPPUNIT_ASSERT(     total>oneGb   &&     total<=256*oneGb );        // believable values between 1Gb and 256 Gb
    #ifdef NEVER    //todo rsatin - fails on Windows and Darwin ---v
	CPPUNIT_ASSERT( swapAvail>oneGb   && swapAvail<=256*oneGb );        // believable values between 1Gb and 256 Gb
    #endif          //todo rsatin - fails on Windows and Darwin ---^
    #ifdef WIN32    //todo rsatin - fails on Linux ---v and Windows
 	  //CPPUNIT_ASSERT(  memAvail>oneGb/8 &&  memAvail<=256*oneGb );        // believable values between 1Gb and 16 Gb
    #endif          //todo rsatin - fails on Linux ---^
    return;
}

void UtilTest::endsWithStrTest()
{
	cout<<"=============================="<<endl;

  Verbose::out(1,"UtilTest::endsWithStrTest");

  CPPUNIT_ASSERT(Util::endsWithStr("012345","5")==true);
  CPPUNIT_ASSERT(Util::endsWithStr("012345","6")==false);
  CPPUNIT_ASSERT(Util::endsWithStr("012345","4")==false);
  CPPUNIT_ASSERT(Util::endsWithStr("012345","")==true);
  CPPUNIT_ASSERT(Util::endsWithStr("062345","6")==false);

  //
  CPPUNIT_ASSERT(Util::endsWithStr("012345abc","5",3)==true);
  CPPUNIT_ASSERT(Util::endsWithStr("012345abc","6",3)==false);
  CPPUNIT_ASSERT(Util::endsWithStr("012345abc","4",3)==false);
  CPPUNIT_ASSERT(Util::endsWithStr("012345abc","4",3)==false);
  CPPUNIT_ASSERT(Util::endsWithStr("012345abc","",3)==true);
  CPPUNIT_ASSERT(Util::endsWithStr("062345abc","6",3)==false);

  //CPPUNIT_ASSERT(false);

	cout<<"=============================="<<endl;
}
