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
 * @class AffxArrayTest
 * @brief cppunit class for testing AffxArray functions.
 * last change by vliber on 05/26/09
 */

class AffxArrayTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(AffxArrayTest);
    CPPUNIT_TEST(compareTest);
	CPPUNIT_TEST(compareDescendingTest);
	CPPUNIT_TEST(mainFuncionsTest);
	CPPUNIT_TEST(BinarySortingTest);
	CPPUNIT_TEST_SUITE_END();

   public:
  
    void compareTest();
	void compareDescendingTest();
	void mainFuncionsTest();
	void BinarySortingTest();
};
CPPUNIT_TEST_SUITE_REGISTRATION(AffxArrayTest);

void AffxArrayTest::compareTest()
{
	cout<<endl;
	cout<<endl;
	Verbose::out(1, "***AffxArrayTest testcases***");
    Verbose::out(1, "AffxArrayTest::compareTest");
    AffxArray<int> ar;
	//static int compare(AffxString& iThis, const AffxString& iThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareCase("abc","bbc")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareCase("abc","abc")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareCase("abc","Abc")==1);
	AffxString s="abc";
	CPPUNIT_ASSERT(AffxArray<int>::compareCase(s,"bbc")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareCase(s,"abc")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareCase(s,"Abc")==1);
    
	//static int compareNoCase(const AffxString& iThis, const AffxString& iThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareNoCase("abc","Bbc")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareNoCase("abc","abc")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareNoCase("cbc","abc")==1);

	// static int compare(bool bThis, bool bThat)
	CPPUNIT_ASSERT(AffxArray<int>::compare(true,true)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(false,true)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(true,false)==1);
	
	// static int compare(int iThis, int iThat)
    CPPUNIT_ASSERT(AffxArray<int>::compare(111,111)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(0,0)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-12,-12)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(111,112)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-112,-111)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-111,-112)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(112,111)==1);

	//static int compare(unsigned int iThis, unsigned int iThat)
	CPPUNIT_ASSERT(AffxArray<int>::compare(111,111)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(0,0)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(111,112)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(112,111)==1);

	// static int compare(double dThis, double dThat)
	CPPUNIT_ASSERT(AffxArray<int>::compare(111.55,111.55)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(0.0,0.0)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-12.4,-12.4)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compare(112.234,112.235)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-112.235,-112.234)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(-111.234,-111.235)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compare(112.235,112.234)==1);

	//static int compareAsInt(const AffxString& strThis, const AffxString& strThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("111","111")==0); 
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("0","0")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("-12","-12")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("111","112")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("-112","-111")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("-111","-112")==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsInt("112","111")==1);
}

void AffxArrayTest::compareDescendingTest()
{
	Verbose::out(1, "AffxArrayTest::compareDescendingTest");
    AffxArray<int> ar;

	//static int compareCaseDescending(AffxString& iThis, const AffxString& iThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending("abc","bbc")==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending("abc","abc")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending("abc","Abc")==-1);
	AffxString s="abc";
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending(s,"bbc")==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending(s,"abc")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareCaseDescending(s,"Abc")==-1);
	
	// static int compareDescending(bool bThis, bool bThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(false,true)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(true,true)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(true,false)==-1);
	
	// static int compareDescending(int iThis, int iThat)
    CPPUNIT_ASSERT(AffxArray<int>::compareDescending(111,111)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(0,0)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-12,-12)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(111,112)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-112,-111)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-111,-112)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(112,111)==-1);

	// static int compareDescending(double dThis, double dThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(111.55,111.55)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(0.0,0.0)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-12.4,-12.4)==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(112.234,112.235)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-112.235,-112.234)==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(-111.234,-111.235)==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareDescending(112.235,112.234)==-1);

	//static int compareAsIntDescending(const AffxString& strThis, const AffxString& strThat)
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("111","111")==0); 
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("0","0")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("-12","-12")==0);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("111","112")==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("-112","-111")==1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("-111","-112")==-1);
	CPPUNIT_ASSERT(AffxArray<int>::compareAsIntDescending("112","111")==-1);
}

void AffxArrayTest::mainFuncionsTest()
{
	Verbose::out(1, "AffxArrayTest::mainFuncionsTest");
	//empty
    AffxArray<double> ar;
    CPPUNIT_ASSERT(ar.size()==0);
	CPPUNIT_ASSERT(ar.getCount()==0);
	//add some data
	int i=0;
	double d=1.12;
	do{
	 d=d+0.01;
	 ar.add(new double(d));
     i++;
	}
	while(i<10);
	CPPUNIT_ASSERT(ToStr(*ar[9])=="1.22");
	CPPUNIT_ASSERT(ar.size()==10);
	CPPUNIT_ASSERT(ar.getCount()==10);
	//make empty using nullAll()
	ar.nullAll();
	CPPUNIT_ASSERT(ar.size()==0);
	CPPUNIT_ASSERT(ar.getCount()==0);
	//add data again
	i=0;
	do{
	 d=d+0.01;
     ar.add(new double(d));
     i++;
	}
	while(i<10);
	CPPUNIT_ASSERT(ToStr(*ar[9])=="1.32");
	CPPUNIT_ASSERT(ar.size()==10);
	CPPUNIT_ASSERT(ar.getCount()==10);
	//make empty using deleteAll()
	ar.deleteAll();
	CPPUNIT_ASSERT(ar.size()==0);
	CPPUNIT_ASSERT(ar.getCount()==0);
	//add data again
	i=0;
	do{
	 d=d+0.01;
     ar.add(new double(d));
     i++;
	}
	while(i<10);
	CPPUNIT_ASSERT(ToStr(*ar[9])=="1.42");
	CPPUNIT_ASSERT(ar.size()==10);
	CPPUNIT_ASSERT(ar.getCount()==10);

	//TYPE*& at(int iIndex)
	CPPUNIT_ASSERT(ToStr(*(ar.at(9)))=="1.42");

	//void setAt(int iIndex, TYPE* obj)
	double d1=1.51;
	ar.setAt(9,new double(d1));
	CPPUNIT_ASSERT(ToStr(*(ar.at(9)))=="1.51");
	POSITIVE_TEST(ar.setAt(20,new double(d1))); //out of bounds

	//TYPE* getAt(int iIndex)
    CPPUNIT_ASSERT(ToStr(*(ar.getAt(9)))=="1.51");

	//void swap(int iIndex1, int iIndex2)
	ar.swap(0,9);
	CPPUNIT_ASSERT(ToStr(*(ar.getAt(9)))=="1.33");
	CPPUNIT_ASSERT(ToStr(*(ar.getAt(0)))=="1.51");
 
	//void shuffle()
	POSITIVE_TEST(ar.shuffle());
	POSITIVE_TEST(ar.shuffle());
	CPPUNIT_ASSERT(ar.getCount()==10);

	//void removeAt(int iIndex) void deleteAt(int iIndex)
    POSITIVE_TEST(ar.removeAt(4));
	CPPUNIT_ASSERT(ar.getCount()==9);
	POSITIVE_TEST(ar.deleteAt(4));
	CPPUNIT_ASSERT(ar.getCount()==8);
}

void AffxArrayTest::BinarySortingTest()
{
	Verbose::out(1, "AffxArrayTest::BinarySortingTest");
	//empty
    AffxArray<AffxString> ar;
    CPPUNIT_ASSERT(ar.size()==0);
	CPPUNIT_ASSERT(ar.getCount()==0);
	//add some data
	int i=0;
	double d=1.12;
	
	do{
	 d=d+0.01;
	 ar.add(new AffxString(ToStr(d)));
     i++;
	}
	while(i<10);	
	AffxString test="1.16";
    CPPUNIT_ASSERT(ar.binarySearch(test,0)==3);

	POSITIVE_TEST(ar.shuffle());
	POSITIVE_TEST(ar.quickSort());
	CPPUNIT_ASSERT(ar.binarySearch(test,0)==3);
    //shuffle, sort, and search again
	POSITIVE_TEST(ar.shuffle());
	POSITIVE_TEST(ar.quickSort());
	CPPUNIT_ASSERT(ar.binarySearch(test,0)==3);
    //shuffle, sort, and search again
	POSITIVE_TEST(ar.shuffle());
	POSITIVE_TEST(ar.quickSort());
	CPPUNIT_ASSERT(ar.binarySearch(test,0)==3);

	//void quickSort(int iFrom, int iTo, int iCompareCode)
	ar.quickSort(0,ar.getCount()-1,0);
	POSITIVE_TEST(ar.shuffle());
	POSITIVE_TEST(ar.quickSort(3,6,0));
	//test that (6)>(5)>(4)>(3)
	CPPUNIT_ASSERT(*ar.getAt(3)<*ar.getAt(4));
	CPPUNIT_ASSERT(*ar.getAt(4)<*ar.getAt(5));
	CPPUNIT_ASSERT(*ar.getAt(5)<*ar.getAt(6));
}



	
	
	


