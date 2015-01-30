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



#include "util/AffxMultiDimensionalArray.h"
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
 * @class AffxMultiDimensionalArrayTest
 * @brief cppunit class for testing AffxMultiDimensionalArray functions.
 * last change by vliber on 05/26/09
 */

class AffxMultiDimensionalArrayTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(AffxMultiDimensionalArrayTest);
    CPPUNIT_TEST(constructormainFunctionTest);
	CPPUNIT_TEST(differentFunctionTest);
	CPPUNIT_TEST_SUITE_END();

   public:
  
    void constructormainFunctionTest();
	void differentFunctionTest();
    
};
CPPUNIT_TEST_SUITE_REGISTRATION(AffxMultiDimensionalArrayTest);

void AffxMultiDimensionalArrayTest::constructormainFunctionTest()
{
	cout<<endl;
	Verbose::out(1, "***AffxMultiDimensionalArrayTest testcases***");
    Verbose::out(1, "AffxMultiDimensionalArrayTest::constructormainFunctionTest");
	AffxMultiDimensionalArray<int> ar1;
	CPPUNIT_ASSERT(ar1.getCount()==1);
	CPPUNIT_ASSERT(ar1.getXDimension()==1);
	CPPUNIT_ASSERT(ar1.getYDimension()==1);
	CPPUNIT_ASSERT(ar1.getZDimension()==1);
	//
	AffxMultiDimensionalArray<int> ar2(2);
	CPPUNIT_ASSERT(ar2.getCount()==2);
	CPPUNIT_ASSERT(ar2.getXDimension()==2);
	CPPUNIT_ASSERT(ar2.getYDimension()==1);
	CPPUNIT_ASSERT(ar2.getZDimension()==1);
	//
	AffxMultiDimensionalArray<int> ar3(2,3);
	CPPUNIT_ASSERT(ar3.getCount()==6);
	CPPUNIT_ASSERT(ar3.getXDimension()==2);
	CPPUNIT_ASSERT(ar3.getYDimension()==3);
	CPPUNIT_ASSERT(ar3.getZDimension()==1);
	//
	AffxMultiDimensionalArray<int> ar4(2,3,4);
	CPPUNIT_ASSERT(ar4.getCount()==24);
	CPPUNIT_ASSERT(ar4.getXDimension()==2);
	CPPUNIT_ASSERT(ar4.getYDimension()==3);
	CPPUNIT_ASSERT(ar4.getZDimension()==4);
    //
    AffxMultiDimensionalArray<int> ar5(ar4);
	CPPUNIT_ASSERT(ar5.getCount()==24);
	CPPUNIT_ASSERT(ar5.getXDimension()==2);
	CPPUNIT_ASSERT(ar5.getYDimension()==3);
	CPPUNIT_ASSERT(ar5.getZDimension()==4);

	AffxMultiDimensionalArray<int> ar6(2,3,4);
	ar6.set(0,0,0,9);
    ar6.set(0,0,1,10);
	ar6.set(0,0,2,11);
	ar6.set(0,0,3,12);
	ar6.set(0,1,0,13);
    ar6.set(0,1,1,14);
	ar6.set(0,1,2,15);
	ar6.set(0,1,3,16);
	ar6.set(0,2,0,17);
    ar6.set(0,2,1,18);
	ar6.set(0,2,2,19);
	ar6.set(0,2,3,20);
	ar6.set(1,0,0,21);
    ar6.set(1,0,1,22);
	ar6.set(1,0,2,23);
	ar6.set(1,0,3,24);
    ar6.set(1,1,0,25);
	ar6.set(1,1,1,26);
	ar6.set(1,1,2,27);
	ar6.set(1,1,3,28);
    ar6.set(1,2,0,29);
	ar6.set(1,2,1,30);
	ar6.set(1,2,2,31);
	ar6.set(1,2,3,32);
	CPPUNIT_ASSERT(ar6.get(0,0,0)==9);
	CPPUNIT_ASSERT(ar6.get(0,0,1)==10);
	CPPUNIT_ASSERT(ar6.get(0,0,2)==11);
	CPPUNIT_ASSERT(ar6.get(0,0,3)==12);
	CPPUNIT_ASSERT(ar6.get(0,1,0)==13);
	CPPUNIT_ASSERT(ar6.get(0,1,1)==14);
	CPPUNIT_ASSERT(ar6.get(0,1,2)==15);
	CPPUNIT_ASSERT(ar6.get(0,1,3)==16);
	CPPUNIT_ASSERT(ar6.get(0,2,0)==17);
	CPPUNIT_ASSERT(ar6.get(0,2,1)==18);
	CPPUNIT_ASSERT(ar6.get(0,2,2)==19);
	CPPUNIT_ASSERT(ar6.get(0,2,3)==20);
	CPPUNIT_ASSERT(ar6.get(1,0,0)==21);
	CPPUNIT_ASSERT(ar6.get(1,0,1)==22);
	CPPUNIT_ASSERT(ar6.get(1,0,2)==23);
	CPPUNIT_ASSERT(ar6.get(1,0,3)==24);
	CPPUNIT_ASSERT(ar6.get(1,1,0)==25);
	CPPUNIT_ASSERT(ar6.get(1,1,1)==26);
	CPPUNIT_ASSERT(ar6.get(1,1,2)==27);
	CPPUNIT_ASSERT(ar6.get(1,1,3)==28);
	CPPUNIT_ASSERT(ar6.get(1,2,0)==29);
	CPPUNIT_ASSERT(ar6.get(1,2,1)==30);
	CPPUNIT_ASSERT(ar6.get(1,2,2)==31);
	CPPUNIT_ASSERT(ar6.get(1,2,3)==32);

    //void insertAt(int x, TYPE t)  TYPE* getPointer()
	ar2.set(0,555);
	ar2.set(1,666);
	CPPUNIT_ASSERT(ar2.get(0)==555);
	CPPUNIT_ASSERT(ar2.get(1)==666);
	ar2.insertAt(0,777);
	CPPUNIT_ASSERT(ar2.get(0)==777);
	CPPUNIT_ASSERT(*(ar2.getPointer())==777);
	CPPUNIT_ASSERT(*(ar2.getPointer()+1)==555);
	CPPUNIT_ASSERT(*(ar2.getPointer()+2)==666);
	CPPUNIT_ASSERT(ar2.getCount()==3);
	ar2.removeAt(1);
	CPPUNIT_ASSERT(ar2.getCount()==2);
	CPPUNIT_ASSERT(*(ar2.getPointer())==777);
	CPPUNIT_ASSERT(*(ar2.getPointer()+1)==666);
	ar2.removeAt(1);
	CPPUNIT_ASSERT(ar2.getCount()==1);
	CPPUNIT_ASSERT(*(ar2.getPointer())==777);
	ar2.removeAt(0);
	CPPUNIT_ASSERT(ar2.getCount()==0);

	//bool equals(AffxMultiDimensionalArray& that)
	AffxMultiDimensionalArray<int> ar7(ar6);
	CPPUNIT_ASSERT(ar6.equals(ar7));
	ar7.set(1,2,3,33);
	CPPUNIT_ASSERT(!ar6.equals(ar7));

	//int compareTo(AffxMultiDimensionalArray& that)
    CPPUNIT_ASSERT(ar6.compareTo(ar7)==-1);//difference by value
	AffxMultiDimensionalArray<int> ar8(3);
	ar8.set(0,555);
	ar8.set(1,666);
	ar8.set(2,777);
	CPPUNIT_ASSERT(ar8.getCount()==3);
	AffxMultiDimensionalArray<int> ar9(ar8);
	CPPUNIT_ASSERT(ar9.getCount()==3);
    CPPUNIT_ASSERT(ar8.equals(ar9));
	CPPUNIT_ASSERT(ar8.compareTo(ar9)==0);
	ar9.removeAt(2);
	CPPUNIT_ASSERT(ar9.getCount()==2);
	CPPUNIT_ASSERT(ar8.compareTo(ar9)==1);//difference by count

	//AffxMultiDimensionalArray<TYPE> rev()
    AffxMultiDimensionalArray<int> ar10=ar8.rev();
    CPPUNIT_ASSERT(ar10.get(0)==777);
	CPPUNIT_ASSERT(ar10.get(1)==666);
	CPPUNIT_ASSERT(ar10.get(2)==555);

	//AffxMultiDimensionalArray<TYPE> diff()
	ar10.set(2,999);
    AffxMultiDimensionalArray<int> ar11=ar10.diff();
	CPPUNIT_ASSERT(ar11.getCount()==2);
	CPPUNIT_ASSERT(ar11.get(0)==-111);
	CPPUNIT_ASSERT(ar11.get(1)==333);

	//void* getPointer(int x)
	CPPUNIT_ASSERT(*(int*)ar11.getPointer(0)==-111);//cust to *int and after that dereference
	CPPUNIT_ASSERT(*(int*)ar11.getPointer(1)==333);
}


void AffxMultiDimensionalArrayTest::differentFunctionTest()
{
	Verbose::out(1, "AffxMultiDimensionalArrayTest::differentFunctionTest");
	AffxMultiDimensionalArray<double> ar(10);
	ar.set(0,1.92);
	ar.set(1,1.12);
	ar.set(2,2.32);
	ar.set(3,1.42);
	ar.set(4,3.52);
	ar.set(5,1.62);
	ar.set(6,3.32);
	ar.set(7,1.22);
	ar.set(8,1.72);
	ar.set(9,1.62);
	CPPUNIT_ASSERT(ar.get(0)==1.92);
	CPPUNIT_ASSERT(ar.get(9)==1.62);
	
	//void quickSort()
	AffxMultiDimensionalArray<double> ar1(ar);
	ar1.quickSort();
	CPPUNIT_ASSERT(ar1.get(0)==1.12);
	CPPUNIT_ASSERT(ar1.get(9)==3.52);
	
	
	//double median()
    CPPUNIT_ASSERT(ar1.median()==1.67); 
	ar1.insertAt(0,2.42);
	CPPUNIT_ASSERT(ar1.median()==1.72);

	//double median(int iLength)
    AffxMultiDimensionalArray<double> ar2(ar);
	CPPUNIT_ASSERT(ar2.median(5)==1.92);
	CPPUNIT_ASSERT(ar2.median(6)==1.77);
    
	//double mean()
     AffxMultiDimensionalArray<double> ar3(ar);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(ar3.mean(),1.98, 0.000001);

}

