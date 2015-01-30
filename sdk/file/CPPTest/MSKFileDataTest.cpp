////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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


#include "file/CPPTest/MSKFileDataTest.h"
//
#include "file/MSKFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CMSKFileDataTest );

/* Test.MSK contents
		Hu6800
		[Comp]
		ABC_at
		XYZ_at
		[Call]
		A28102_at	17,18,19,20
		AB000381_s_at	7
		AC000064_cds2_at	18,19,20
		AC002086_at	1-10,15-20
		AC00001_at	1-20
*/
#define TEST_MSK_NUM_PROBE_SET_INDICIES_COUNT 5
#define TEST_MSK_NUM_PROBE_SET_COUNT 2
#define TEST_MSK_ARRAY_TYPE "Hu6800"

#ifdef _MSC_VER
#define TEST_MSK_FILE_NAME ".\\data\\test.msk"
#define NODATA_MSK_FILE_NAME ".\\data\\nodata.msk"
#else
#define TEST_MSK_FILE_NAME "./data/test.MSK"
#define NODATA_MSK_FILE_NAME "./data/nodata.MSK"
#endif

void CMSKFileDataTest::setUp()
{
}

void CMSKFileDataTest::tearDown()
{
}

void CMSKFileDataTest::testCreation()
{
	affxmsk::CMSKFileData msk;
	CPPUNIT_ASSERT( 1 );
}

void CMSKFileDataTest::testproperty_FileName()
{
	affxmsk::CMSKFileData msk;
	std::string path = "test";
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == msk.GetFileName() );
}

void CMSKFileDataTest::testmethod_ProbeSetIndiciesListCount()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();
	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == TEST_MSK_NUM_PROBE_SET_INDICIES_COUNT );
}

void CMSKFileDataTest::testproperty_ArrayType()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();
	std::string atype = msk.GetArrayType();
	CPPUNIT_ASSERT( atype == TEST_MSK_ARRAY_TYPE );
}

void CMSKFileDataTest::testmethod_Exists()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT( msk.Exists() == true );
}

void CMSKFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxmsk::CMSKFileData msk;
	std::string path = "test";
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT( msk.Exists() == false );
}

void CMSKFileDataTest::testmethod_Read()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	CPPUNIT_ASSERT( msk.Read() );
}

void CMSKFileDataTest::testmethod_ReadWhenNoDataIsInTheFile()
{
	affxmsk::CMSKFileData msk;
	std::string path = NODATA_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();
	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 0 );
	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == 0 );
}

void CMSKFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxmsk::CMSKFileData msk;
	std::string path = "test";
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT( msk.Read() == false );
}

void CMSKFileDataTest::testmethod_Clear()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();
	msk.Clear();
	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 0 );
	std::string atype = msk.GetArrayType();
	CPPUNIT_ASSERT( atype.length() == 0 );
}

void CMSKFileDataTest::testmethod_GetProbeSetIndiciesIterators()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();

	affxmsk::ProbeSetIndiciesType ind;
	std::list<int>::iterator index_iter;
	affxmsk::ProbeSetIndiciesListConstIt iter;
	affxmsk::ProbeSetIndiciesListConstIt begin;
	affxmsk::ProbeSetIndiciesListConstIt end;

	msk.GetProbeSetIndiciesIterators(begin, end);
	iter = begin;


	// Check the first entry. The line in the MSK file is: A28102_at	17,18,19,20
	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "A28102_at" );
	index_iter = ind.indicies.begin();
	CPPUNIT_ASSERT( *index_iter == 16 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 17 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 18 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 19 );
	++index_iter;
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );

	// Advance to the fourth entry. The line in the MSK file is: AC002086_at	1-10,15-20
	++iter;
	++iter;
	++iter;
	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "AC002086_at" );
	index_iter = ind.indicies.begin();
	for (int i=0; i<10; i++)
	{
		CPPUNIT_ASSERT( *index_iter == i );
		++index_iter;
	}
	for (int i=14; i<20; i++)
	{
		CPPUNIT_ASSERT( *index_iter == i );
		++index_iter;
	}
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );

	// Check the last entry. The line in the MSK file is: AC00001_at	1-20
	++iter;
	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "AC00001_at" );
	index_iter = ind.indicies.begin();
	for (int i=0; i<20; i++)
	{
		CPPUNIT_ASSERT( *index_iter == i );
		++index_iter;
	}
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );

	++iter;
	CPPUNIT_ASSERT( iter == end );
}

void CMSKFileDataTest::testmethod_GetProbeSetListCount()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();
	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == TEST_MSK_NUM_PROBE_SET_COUNT );
}

void CMSKFileDataTest::testmethod_GetProbeSetListIterators()
{
	affxmsk::CMSKFileData msk;
	std::string path = TEST_MSK_FILE_NAME;
	msk.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", msk.Exists() == true );
	msk.Read();

	affxmsk::ProbeSetListConstIt iter;
	affxmsk::ProbeSetListConstIt begin;
	affxmsk::ProbeSetListConstIt end;

	msk.GetProbeSetIterators(begin, end);

	iter = begin;
	CPPUNIT_ASSERT( *iter == "ABC_at" );

	++iter;
	CPPUNIT_ASSERT( *iter == "XYZ_at" );

	++iter;
	CPPUNIT_ASSERT( iter == end );
}


