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


#include "file/CPPTest/EXPFileDataTest.h"
//
#include "file/EXPFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CEXPFileDataTest );

#define FULL_EXP "./data/full.EXP"
#define BARE_EXP "./data/bare.EXP"
#define ABORT_EXP "./data/hyb_abort.EXP"
#define NODATA_EXP "./data/nodata.EXP"

void CEXPFileDataTest::setUp()
{
}

void CEXPFileDataTest::tearDown()
{
}

void CEXPFileDataTest::testCreation()
{
	affxexp::CEXPFileData EXP;
	CPPUNIT_ASSERT( 1 );
}

void CEXPFileDataTest::testproperty_FileName()
{
	affxexp::CEXPFileData EXP;
	std::string path = "test";
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == EXP.GetFileName() );
}

void CEXPFileDataTest::testproperty_ArrayType()
{
	affxexp::CEXPFileData EXP;
	std::string path = FULL_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", EXP.Exists() == true );
	EXP.Read();
	CPPUNIT_ASSERT( EXP.GetArrayType() == "Test3" );

	EXP.Clear();
	CPPUNIT_ASSERT( EXP.GetArrayType() == "" );

	EXP.SetArrayType("Test");
	CPPUNIT_ASSERT( EXP.GetArrayType() == "Test" );
}

void CEXPFileDataTest::testmethod_Exists()
{
	affxexp::CEXPFileData EXP;
	std::string path = FULL_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT( EXP.Exists() == true );
}

void CEXPFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxexp::CEXPFileData EXP;
	std::string path = "test";
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT( EXP.Exists() == false );
}

void CEXPFileDataTest::testmethod_Read_full()
{
	affxexp::CEXPFileData EXP;
	std::string path = FULL_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", EXP.Exists() == true );
	CPPUNIT_ASSERT( EXP.Read() );

	CPPUNIT_ASSERT( EXP.GetArrayType() == "Test3" );

	TagValuePairType param;
	TagValuePairTypeList::iterator it;
	TagValuePairTypeList &scan = EXP.GetScanParameters();
	CPPUNIT_ASSERT( scan.size() == 7);
	it = scan.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Pixel Size");
	CPPUNIT_ASSERT( param.Value == "3");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Filter");
	CPPUNIT_ASSERT( param.Value == "570");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Scan Temperature");
	CPPUNIT_ASSERT( param.Value == "200");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Scan Date");
	CPPUNIT_ASSERT( param.Value == "Aug 01 2002 08:37AM");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Scanner ID");
	CPPUNIT_ASSERT( param.Value == "123123");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Number of Scans");
	CPPUNIT_ASSERT( param.Value == "1");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Scanner Type");
	CPPUNIT_ASSERT( param.Value == "HP");



	TagValuePairTypeList &hyb = EXP.GetHybParameters();
	CPPUNIT_ASSERT( hyb.size() == 19);
	it = hyb.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Protocol");
	CPPUNIT_ASSERT( param.Value == "EukGE-WS1v4");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash A1 Recovery Mixes");
	CPPUNIT_ASSERT( param.Value == "0");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash A1 Temperature (C)");
	CPPUNIT_ASSERT( param.Value == "25");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Number of Wash A1 Cycles");
	CPPUNIT_ASSERT( param.Value == "10");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Mixes per Wash A1 Cycle");
	CPPUNIT_ASSERT( param.Value == "2");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash B Recovery Mixes");
	CPPUNIT_ASSERT( param.Value == "0");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash B Temperature (C)");
	CPPUNIT_ASSERT( param.Value == "50");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Number of Wash B Cycles");
	CPPUNIT_ASSERT( param.Value == "4");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Mixes per Wash B Cycle");
	CPPUNIT_ASSERT( param.Value == "15");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Stain Temperature (C)");
	CPPUNIT_ASSERT( param.Value == "25");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Stain Time (seconds)");
	CPPUNIT_ASSERT( param.Value == "1800");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash A2 Recovery Mixes");
	CPPUNIT_ASSERT( param.Value == "0");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Wash A2 Temperature (C)");
	CPPUNIT_ASSERT( param.Value == "25");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Number of Wash A2 Cycles");
	CPPUNIT_ASSERT( param.Value == "10");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Mixes per Wash A2 Cycle");
	CPPUNIT_ASSERT( param.Value == "4");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Holding Temperature (C)");
	CPPUNIT_ASSERT( param.Value == "25");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Station");
	CPPUNIT_ASSERT( param.Value == "1");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Module");
	CPPUNIT_ASSERT( param.Value == "1");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Hybridize Date");
	CPPUNIT_ASSERT( param.Value == "Oct 08 2002 03:31PM");





	TagValuePairTypeList &samp = EXP.GetSampleParameters();
	CPPUNIT_ASSERT( samp.size() == 8);
	it = samp.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Chip Lot");
	CPPUNIT_ASSERT( param.Value == "123");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Operator");
	CPPUNIT_ASSERT( param.Value == "ljevon");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Sample Type");
	CPPUNIT_ASSERT( param.Value == "s_type");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Description");
	CPPUNIT_ASSERT( param.Value == "Demo data");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Project");
	CPPUNIT_ASSERT( param.Value == "s_proj");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Comments");
	CPPUNIT_ASSERT( param.Value == "None");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Solution Type");
	CPPUNIT_ASSERT( param.Value == "sol type");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Solution Lot");
	CPPUNIT_ASSERT( param.Value == "123123");

}

void CEXPFileDataTest::testmethod_Read_bare()
{
	affxexp::CEXPFileData EXP;
	std::string path = BARE_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", EXP.Exists() == true );
	CPPUNIT_ASSERT( EXP.Read() );

	CPPUNIT_ASSERT( EXP.GetArrayType() == "Mapping10K_Xba131" );

	TagValuePairType param;
	TagValuePairTypeList::iterator it;
	TagValuePairTypeList &scan = EXP.GetScanParameters();
	CPPUNIT_ASSERT( scan.size() == 1);
	it = scan.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Scanner Type");
	CPPUNIT_ASSERT( param.Value == "M18");

	TagValuePairTypeList &hyb = EXP.GetHybParameters();
	CPPUNIT_ASSERT( hyb.size() == 0);

	TagValuePairTypeList &samp = EXP.GetSampleParameters();
	CPPUNIT_ASSERT( samp.size() == 1);
	it = samp.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Operator");
	CPPUNIT_ASSERT( param.Value == "sgelle");
}

void CEXPFileDataTest::testmethod_Read_hyb_abort()
{
	affxexp::CEXPFileData EXP;
	std::string path = ABORT_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", EXP.Exists() == true );
	CPPUNIT_ASSERT( EXP.Read() );

	CPPUNIT_ASSERT( EXP.GetArrayType() == "HG-U133A" );

	TagValuePairType param;
	TagValuePairTypeList::iterator it;
	TagValuePairTypeList &scan = EXP.GetScanParameters();
	CPPUNIT_ASSERT( scan.size() == 0);

	TagValuePairTypeList &hyb = EXP.GetHybParameters();
	CPPUNIT_ASSERT( hyb.size() == 7);

	it = hyb.begin();
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Protocol");
	CPPUNIT_ASSERT( param.Value == "EukGE-WS1v4");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "REMOVE VIAL");
	CPPUNIT_ASSERT( param.Value == "");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "            20°C");
	CPPUNIT_ASSERT( param.Value == "");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Aborted by user");
	CPPUNIT_ASSERT( param.Value == "");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Station");
	CPPUNIT_ASSERT( param.Value == "1");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Module");
	CPPUNIT_ASSERT( param.Value == "2");

	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Tag == "Hybridize Date");
	CPPUNIT_ASSERT( param.Value == "Mar 28 2003 09:52AM");



	TagValuePairTypeList &samp = EXP.GetSampleParameters();
	CPPUNIT_ASSERT( samp.size() == 0);
}

void CEXPFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxexp::CEXPFileData EXP;
	std::string path = "test";
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT( EXP.Read() == false );
}

void CEXPFileDataTest::testmethod_Clear()
{
	affxexp::CEXPFileData EXP;
	std::string path = FULL_EXP;
	EXP.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", EXP.Exists() == true );
	EXP.Read();
	EXP.Clear();
	CPPUNIT_ASSERT( EXP.GetArrayType().length() == 0 );
	CPPUNIT_ASSERT( EXP.GetScanParameters().size() == 0);
	CPPUNIT_ASSERT( EXP.GetHybParameters().size() == 0);
	CPPUNIT_ASSERT( EXP.GetSampleParameters().size() == 0);

}
