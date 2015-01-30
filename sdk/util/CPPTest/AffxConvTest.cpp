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

#include "util/AffxConv.h"
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
 * @class AffxConvTest
 * @brief cppunit class for testing AffxConv functions.
 * last change by vliber on 05/26/09
 */

class AffxConvTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(AffxConvTest);
    CPPUNIT_TEST(getIntTest);
	CPPUNIT_TEST(getDoubleTest);
	CPPUNIT_TEST(formatStringTest);
	CPPUNIT_TEST(roundDoubleTest);
    CPPUNIT_TEST(getMonthTest);
	CPPUNIT_TEST(likeStringTest);
    CPPUNIT_TEST_SUITE_END();

   public:
  
    void getIntTest();
    void getDoubleTest();
	void formatStringTest();
	void roundDoubleTest();
	void getMonthTest();
    void likeStringTest();
};
CPPUNIT_TEST_SUITE_REGISTRATION(AffxConvTest );

void AffxConvTest::getIntTest()
{
	cout<<endl;
	cout<<endl;
	Verbose::out(1, "***AffxConvTest testcases***");
    Verbose::out(1, "AffxConvTest::getIntTest");
	CPPUNIT_ASSERT(::getInt(-2147483647) == "-2147483647");
	CPPUNIT_ASSERT(::getInt(0,true)=="0");
    CPPUNIT_ASSERT(::getInt(2147483647) == "2147483647");
	CPPUNIT_ASSERT(::getInt(-1234567,true)=="-1,234,567");
	CPPUNIT_ASSERT(::getInt(1234567,false)=="1234567");
	CPPUNIT_ASSERT(::getUnsignedInt(4294967295,true)=="4,294,967,295");
	CPPUNIT_ASSERT(::getUnsignedInt(4294967295,false)=="4294967295");
	AffxString s="-2147483647";
	CPPUNIT_ASSERT(::getInt(s)==-2147483647);
	s="2147483647";
	CPPUNIT_ASSERT(::getInt(s)==2147483647);
	CPPUNIT_ASSERT(::getInt(1,(unsigned int)2)=="01"); 
	CPPUNIT_ASSERT(::getInt(1,(unsigned int)1)=="1"); 
	CPPUNIT_ASSERT(::getInt(11,(unsigned int)3)=="011"); 
	CPPUNIT_ASSERT(::getInt(11,(unsigned int)2)=="11"); 
}
void AffxConvTest::getDoubleTest()
{
	Verbose::out(1, "AffxConvTest::getDoubleTest");
	CPPUNIT_ASSERT(::getDouble(-4294967295.08) == "-4294967295.08");
	CPPUNIT_ASSERT(::getDouble(0,true)=="0");
    CPPUNIT_ASSERT(::getDouble(4294967295.09) == "4294967295.09");
	CPPUNIT_ASSERT(::getDouble(-1234567.0844,3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::getDouble(-1234567.8344,3,false)=="-1234567.834");
	CPPUNIT_ASSERT(::getDouble(1234567.8344,3,false)=="1234567.834");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,1,true)=="-1,234,567.8");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,2,true)=="-1,234,567.83");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,4,false)=="-1234567.8347");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,5,false)=="-1234567.83467");
//	#ifdef WIN32
	CPPUNIT_ASSERT(::getDouble(-1234567.0845,3,true)=="-1,234,567.085");
	CPPUNIT_ASSERT(::getDouble(-1234567.0846,3,true)=="-1,234,567.085");
	CPPUNIT_ASSERT(::getDouble(1234567.0844,3,true)=="1,234,567.084");
	CPPUNIT_ASSERT(::getDouble(1234567.0845,3,true)=="1,234,567.085");
	CPPUNIT_ASSERT(::getDouble(1234567.0846,3,true)=="1,234,567.085"); 
	CPPUNIT_ASSERT(::getDouble(-1234567.8345,3,false)=="-1234567.835");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346,3,false)=="-1234567.835");	
	CPPUNIT_ASSERT(::getDouble(1234567.8345,3,false)=="1234567.835");
	CPPUNIT_ASSERT(::getDouble(1234567.8346,3,false)=="1234567.835");
	CPPUNIT_ASSERT(::getDouble(-1234567.8356678,2,true)=="-1,234,567.84");
	CPPUNIT_ASSERT(::getDouble(-1234567.8366678,2,true)=="-1,234,567.84");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,3,false)=="-1234567.835");
/* 
	#else
    CPPUNIT_ASSERT(::getDouble(-1234567.0845,3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::getDouble(-1234567.0846,3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::getDouble(1234567.0845,3,true)=="1,234,567.084");
	CPPUNIT_ASSERT(::getDouble(1234567.0846,3,true)=="1,234,567.085");
	CPPUNIT_ASSERT(::getDouble(-1234567.8345,3,false)=="-1234567.834");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346,3,false)=="-1234567.834");	
	CPPUNIT_ASSERT(::getDouble(1234567.8345,3,false)=="1234567.834");
	CPPUNIT_ASSERT(::getDouble(1234567.8346,3,false)=="1234567.834");
	CPPUNIT_ASSERT(::getDouble(-1234567.8356678,2,true)=="-1,234,567.84");
	CPPUNIT_ASSERT(::getDouble(-1234567.8366678,2,true)=="-1,234,567.84");
	CPPUNIT_ASSERT(::getDouble(-1234567.8346678,3,false)=="-1234567.835");
	#endif
*/	
	AffxString s="123456789.0087";
	CPPUNIT_ASSERT(::getDouble(s)==123456789.0087);	
}


void AffxConvTest::formatStringTest()
{
	Verbose::out(1, "AffxConvTest::formatStringTest");
	CPPUNIT_ASSERT(::formatString("-1234567.08431",3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::formatString("-1234567.08451",3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::formatString("-1234567.08491",3,true)=="-1,234,567.084");
	CPPUNIT_ASSERT(::formatString("-1234567.08531",3,true)=="-1,234,567.085");
	CPPUNIT_ASSERT(::formatString("-1234567.83431",3,false)=="-1234567.834");
	CPPUNIT_ASSERT(::formatString("-1234567.83451",3,false)=="-1234567.834");
	CPPUNIT_ASSERT(::formatString("-1234567.83491",3,false)=="-1234567.834");
	CPPUNIT_ASSERT(::formatString("-1234567.83531",3,false)=="-1234567.835");
	CPPUNIT_ASSERT(::formatString("-1234567.83531",4,false)=="-1234567.8353");
	CPPUNIT_ASSERT(::formatString("-1234567.83535",4,false)=="-1234567.8353");
	CPPUNIT_ASSERT(::formatString("-1234567.83538",4,false)=="-1234567.8353");
	CPPUNIT_ASSERT(::formatString("-1234567.835317",5,false)=="-1234567.83531");
}
void AffxConvTest::roundDoubleTest()
{
	Verbose::out(1, "AffxConvTest::roundDoubleTest");
	//double d=647.0844567;
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0844567)==-2147483647);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0844567,1)==-2147483647.1);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0844567,1)==2147483647.1);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0845567,2)==-2147483647.08);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0845567,2)==2147483647.08);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846783,3)==-2147483647.085);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846783,3)==2147483647.085);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846783,4)==-2147483647.0847);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846783,4)==2147483647.0847);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846783,5)==-2147483647.08468);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846783,5)==2147483647.08468);

	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846784,6)==-2147483647.084678);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846785,6)==-2147483647.084678);
	CPPUNIT_ASSERT(::roundDouble(-2147483647.0846786,6)==-2147483647.084679);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846784,6)==2147483647.084678);
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846785,6)==2147483647.084678); //todo vliber why last digit 8 not 9
	CPPUNIT_ASSERT(::roundDouble(2147483647.0846786,6)==2147483647.084679);
}
void AffxConvTest::getMonthTest()
{
	Verbose::out(1, "AffxConvTest::getMonthTest");
	CPPUNIT_ASSERT(::getMonth(1)=="Jan");
	CPPUNIT_ASSERT(::getMonth(2)=="Feb");
	CPPUNIT_ASSERT(::getMonth(3)=="Mar");
	CPPUNIT_ASSERT(::getMonth(4)=="Apr");
	CPPUNIT_ASSERT(::getMonth(5)=="May");
	CPPUNIT_ASSERT(::getMonth(6)=="Jun");
	CPPUNIT_ASSERT(::getMonth(7)=="Jul");
	CPPUNIT_ASSERT(::getMonth(8)=="Aug");
	CPPUNIT_ASSERT(::getMonth(9)=="Sep");
	CPPUNIT_ASSERT(::getMonth(10)=="Oct");
	CPPUNIT_ASSERT(::getMonth(11)=="Nov");
	CPPUNIT_ASSERT(::getMonth(12)=="Dec");
	//
	CPPUNIT_ASSERT(::getMonth("")==0);
	CPPUNIT_ASSERT(::getMonth("Jan")==1);
    CPPUNIT_ASSERT(::getMonth("Feb")==2);
	CPPUNIT_ASSERT(::getMonth("Mar")==3);
	CPPUNIT_ASSERT(::getMonth("Apr")==4);
	CPPUNIT_ASSERT(::getMonth("May")==5);
	CPPUNIT_ASSERT(::getMonth("Jun")==6);
	CPPUNIT_ASSERT(::getMonth("Jul")==7);
	CPPUNIT_ASSERT(::getMonth("Aug")==8);
	CPPUNIT_ASSERT(::getMonth("Sep")==9);
	CPPUNIT_ASSERT(::getMonth("Oct")==10);
	CPPUNIT_ASSERT(::getMonth("Nov")==11);
	CPPUNIT_ASSERT(::getMonth("Dec")==12);
}
void AffxConvTest::likeStringTest()
{
	Verbose::out(1, "AffxConvTest::likeStringTest");
	std::string s="abcg'afgsk (jhkjh) k%hk_jh[kjh'";
	CPPUNIT_ASSERT(likeString(s)=="abcg''afgsk _jhkjh_ k[%]hk[_]jh[[]kjh''");
}

