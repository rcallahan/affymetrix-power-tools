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

#include "file/CPPTest/GRCFileDataTest.h"
//
#include "file/GRCFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CGRCFileDataTest );

#define TEST_GRC_FILE "./data/test.grc"

void CGRCFileDataTest::setUp()
{
}

void CGRCFileDataTest::tearDown()
{
}

void CGRCFileDataTest::testCreation()
{
	affxgrc::CGRCFileData grc;
	CPPUNIT_ASSERT( 1 );
}

void CGRCFileDataTest::testproperty_FileName()
{
	affxgrc::CGRCFileData grc;
	std::string path = TEST_GRC_FILE;
	grc.SetFileName(path.c_str());
	CPPUNIT_ASSERT( grc.GetFileName() == TEST_GRC_FILE );
	CPPUNIT_ASSERT( grc.Exists() );
}

void CGRCFileDataTest::testproperty_NumCols()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetColumns() == 2) ;
}

void CGRCFileDataTest::testproperty_NumRows()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetRows() == 2) ;
}

void CGRCFileDataTest::testproperty_NumB1()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetNumB1Probes() == 1) ;
}

void CGRCFileDataTest::testproperty_NumB2()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetNumB2Probes() == 2) ;
}

void CGRCFileDataTest::testproperty_NumNonSynth()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetNumNSProbes() == 1) ;
}

void CGRCFileDataTest::testmethod_Exists()
{
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	CPPUNIT_ASSERT( grc.Exists() == true );
	grc.SetFileName("test");
	CPPUNIT_ASSERT( grc.Exists() == false );
}

void CGRCFileDataTest::testmethod_Read()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	CPPUNIT_ASSERT( grc.Read(grid) == true);
	grc.SetFileName("test");
	CPPUNIT_ASSERT( grc.Read(grid) == false);
}

void CGRCFileDataTest::testmethod_GetB1()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetB1(0).x == 1);
	CPPUNIT_ASSERT( grid.GetB1(0).y == 0);
}

void CGRCFileDataTest::testmethod_GetB2()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetB2(0).x == 0);
	CPPUNIT_ASSERT( grid.GetB2(0).y == 0);
	CPPUNIT_ASSERT( grid.GetB2(1).x == 0);
	CPPUNIT_ASSERT( grid.GetB2(1).y == 1);
}

void CGRCFileDataTest::testmethod_GetNonSynth()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	grc.SetFileName(TEST_GRC_FILE);
	grc.Read(grid);
	CPPUNIT_ASSERT( grid.GetNS(0).x == 1);
	CPPUNIT_ASSERT( grid.GetNS(0).y == 1);
}

void CGRCFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxgrc::CGRCFileData grc;
	std::string path = "test";
	grc.SetFileName(path.c_str());
	CPPUNIT_ASSERT( grc.Exists() == false );
}

void CGRCFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affymetrix_grid_control::GridControlData grid;
	affxgrc::CGRCFileData grc;
	std::string path = "test";
	grc.SetFileName(path.c_str());
	CPPUNIT_ASSERT( grc.Read(grid) == false );
}
