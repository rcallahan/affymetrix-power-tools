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


#include "file/CPPTest/PSIFileDataTest.h"
//
#include "file/PSIFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CPSIFileDataTest );

#define TEST_PSI_FILE_NAME "./data/test.PSI"

void CPSIFileDataTest::setUp()
{
}

void CPSIFileDataTest::tearDown()
{
}

void CPSIFileDataTest::testCreation()
{
	affxpsi::CPSIFileData PSI;
	CPPUNIT_ASSERT( 1 );
}

void CPSIFileDataTest::testproperty_FileName()
{
	affxpsi::CPSIFileData PSI;
	std::string path = "test";
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == PSI.GetFileName() );
}

void CPSIFileDataTest::testmethod_Exists()
{
	affxpsi::CPSIFileData PSI;
	std::string path = TEST_PSI_FILE_NAME;
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT( PSI.Exists() == true );
}

void CPSIFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxpsi::CPSIFileData PSI;
	std::string path = "test";
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT( PSI.Exists() == false );
}

void CPSIFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxpsi::CPSIFileData PSI;
	std::string path = "test";
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT( PSI.Read() == false );
}

void CPSIFileDataTest::testmethod_Clear()
{
	affxpsi::CPSIFileData PSI;
	std::string path = TEST_PSI_FILE_NAME;
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", PSI.Exists() == true );
	PSI.Read();
	PSI.Clear();
	CPPUNIT_ASSERT( PSI.GetProbeSetCount() == 0 );
}

void CPSIFileDataTest::testmethod_Read()
{
	affxpsi::CPSIFileData PSI;
	std::string path = TEST_PSI_FILE_NAME;
	PSI.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", PSI.Exists() == true );
	CPPUNIT_ASSERT( PSI.Read() == true );
	CPPUNIT_ASSERT( PSI.GetProbeSetCount() == 5 );
	CPPUNIT_ASSERT( PSI.GetProbeSetName(0) == "one" );
	CPPUNIT_ASSERT( PSI.GetProbePairs(0) == 1 );
	CPPUNIT_ASSERT( PSI.GetProbeSetName(1) == "two" );
	CPPUNIT_ASSERT( PSI.GetProbePairs(1) == 2 );
	CPPUNIT_ASSERT( PSI.GetProbeSetName(2) == "three" );
	CPPUNIT_ASSERT( PSI.GetProbePairs(2) == 3 );
	CPPUNIT_ASSERT( PSI.GetProbeSetName(3) == "four" );
	CPPUNIT_ASSERT( PSI.GetProbePairs(3) == 4 );
	CPPUNIT_ASSERT( PSI.GetProbeSetName(4) == "five" );
	CPPUNIT_ASSERT( PSI.GetProbePairs(4) == 5 );

}
