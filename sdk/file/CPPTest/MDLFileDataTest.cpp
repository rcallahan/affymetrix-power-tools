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


#include "file/CPPTest/MDLFileDataTest.h"
//
#include "file/MDLFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CMDLFileDataTest );

/* Test.MDL contents
*/
#define TEST_MDL_ANALYSIS_TYPE 0
#define TEST_MDL_NUM_PROBE_COUNT 40
#define TEST_MDL_CHIP_TYPE "HG_U95Av2"
#define TEST_MDL_MODIFIED_DATE "March 29 2005 01:48PM"
#define TEST_MDL_VERSION 1
#define TEST_MDL_AFFINITY_ONE 1.23456f
#define TEST_MDL_AFFINITY_TWO 7.89f

#ifdef _MSC_VER
#define TEST_MDL_FILE_NAME ".\\data\\test.mdl"
#define NODATA_MDL_FILE_NAME ".\\data\\nodata.mdl"
#define WRITE_MDL_FILE_NAME ".\\data\\write.mdl"
#else
#define TEST_MDL_FILE_NAME "./data/test.MDL"
#define NODATA_MDL_FILE_NAME "./data/nodata.MDL"
#define WRITE_MDL_FILE_NAME "./data/write.mdl"
#endif

void CMDLFileDataTest::setUp()
{
}

void CMDLFileDataTest::tearDown()
{
}

////////////////////////////////////////////////
// Test interface object creation
void CMDLFileDataTest::testCreation()
{
	affxmdl::CMDLFileData mdl_file;
	CPPUNIT_ASSERT( 1 );
}

void CMDLFileDataTest::testCreationHeader()
{
	affxmdl::CMDLFileHeader mdl_header;
	CPPUNIT_ASSERT( 1 );
}

void CMDLFileDataTest::testCreationData()
{
	affxmdl::CMDLProbeData mdl_data;
	CPPUNIT_ASSERT( 1 );
}

////////////////////////////////////////////////
// Test header interface properties
void CMDLFileDataTest::testproperty_HeaderVersion()
{
	affxmdl::CMDLFileHeader mdl_header;
	int version = TEST_MDL_VERSION;
	mdl_header.SetVersion(version);
	CPPUNIT_ASSERT( version == mdl_header.GetVersion() );
}

void CMDLFileDataTest::testproperty_HeaderNumProbes()
{
	affxmdl::CMDLFileHeader mdl_header;
	int numprobes = TEST_MDL_NUM_PROBE_COUNT;
	mdl_header.SetNumProbes(numprobes);
	CPPUNIT_ASSERT( numprobes == mdl_header.GetNumProbes() );
}

void CMDLFileDataTest::testproperty_HeaderAnalysisType()
{
	affxmdl::CMDLFileHeader mdl_header;
	int type = TEST_MDL_ANALYSIS_TYPE;
	mdl_header.SetAnalysisType(type);
	CPPUNIT_ASSERT( type == mdl_header.GetAnalysisType() );
}

void CMDLFileDataTest::testproperty_HeaderMasked()
{
	affxmdl::CMDLFileHeader mdl_header;
	bool type = true;
	mdl_header.SetMasked(type);
	CPPUNIT_ASSERT( type == mdl_header.GetMasked() );
}

void CMDLFileDataTest::testproperty_HeaderChipType()
{
	affxmdl::CMDLFileHeader mdl_header;
	std::string chiptype = TEST_MDL_CHIP_TYPE;
	mdl_header.SetChipType(chiptype.c_str());
	CPPUNIT_ASSERT( chiptype == mdl_header.GetChipType() );
}

void CMDLFileDataTest::testproperty_HeaderModifiedDate()
{
	affxmdl::CMDLFileHeader mdl_header;
	std::string date = TEST_MDL_MODIFIED_DATE;
	mdl_header.SetModifiedDate(date.c_str());
	CPPUNIT_ASSERT( date == mdl_header.GetModifiedDate() );
}

void CMDLFileDataTest::testproperty_HeaderEqualOperator()
{
	affxmdl::CMDLFileHeader mdl_header_right;
	affxmdl::CMDLFileHeader mdl_header_left;

	mdl_header_right.SetAnalysisType(TEST_MDL_ANALYSIS_TYPE);
	mdl_header_right.SetChipType(TEST_MDL_CHIP_TYPE);
	mdl_header_right.SetMasked(true);
	mdl_header_right.SetModifiedDate(TEST_MDL_MODIFIED_DATE);
	mdl_header_right.SetNumProbes(TEST_MDL_NUM_PROBE_COUNT);
	mdl_header_right.SetVersion(TEST_MDL_VERSION);

	mdl_header_left = mdl_header_right;

	CPPUNIT_ASSERT( mdl_header_right.GetAnalysisType() == mdl_header_left.GetAnalysisType() );
	CPPUNIT_ASSERT( mdl_header_right.GetChipType() == mdl_header_left.GetChipType() );
	CPPUNIT_ASSERT( mdl_header_right.GetMasked() == mdl_header_left.GetMasked() );
	CPPUNIT_ASSERT( mdl_header_right.GetModifiedDate() == mdl_header_left.GetModifiedDate() );
	CPPUNIT_ASSERT( mdl_header_right.GetNumProbes() == mdl_header_left.GetNumProbes() );
	CPPUNIT_ASSERT( mdl_header_right.GetVersion() == mdl_header_left.GetVersion() );
}

////////////////////////////////////////////////
// Test data interface properties
void CMDLFileDataTest::testproperty_DataAffinity()
{
	affxmdl::CMDLProbeData mdl_data;
	double value = TEST_MDL_AFFINITY_ONE;
	mdl_data.SetAffinity(value);
	CPPUNIT_ASSERT( value == mdl_data.GetAffinity() );
}

void CMDLFileDataTest::testproperty_DataRelativeBkg()
{
	affxmdl::CMDLProbeData mdl_data;
	float value = TEST_MDL_AFFINITY_ONE;
	mdl_data.SetRelativeBkg(value);
	CPPUNIT_ASSERT( value == mdl_data.GetRelativeBkg() );
}

void CMDLFileDataTest::testproperty_DataSaturation()
{
	affxmdl::CMDLProbeData mdl_data;
	float value = TEST_MDL_AFFINITY_ONE;
	mdl_data.SetSaturation(value);
	CPPUNIT_ASSERT( value == mdl_data.GetSaturation() );
}

void CMDLFileDataTest::testproperty_DataOffset()
{
	affxmdl::CMDLProbeData mdl_data;
	float value = TEST_MDL_AFFINITY_ONE;
	mdl_data.SetOffset(value);
	CPPUNIT_ASSERT( value == mdl_data.GetOffset() );
}

////////////////////////////////////////////////
// Test file interface properties
void CMDLFileDataTest::testproperty_Header()
{
	affxmdl::CMDLFileData mdl_file;
	affxmdl::CMDLFileHeader mdl_header;

	mdl_header.SetAnalysisType(TEST_MDL_ANALYSIS_TYPE);
	mdl_header.SetChipType(TEST_MDL_CHIP_TYPE);
	mdl_header.SetMasked(true);
	mdl_header.SetModifiedDate(TEST_MDL_MODIFIED_DATE);
	mdl_header.SetNumProbes(TEST_MDL_NUM_PROBE_COUNT);
	mdl_header.SetVersion(TEST_MDL_VERSION);

	mdl_file.SetHeader(mdl_header);

	CPPUNIT_ASSERT( mdl_header.GetAnalysisType() == mdl_file.GetHeader().GetAnalysisType() );
	CPPUNIT_ASSERT( mdl_header.GetChipType() == mdl_file.GetHeader().GetChipType() );
	CPPUNIT_ASSERT( mdl_header.GetMasked() == mdl_file.GetHeader().GetMasked() );
	CPPUNIT_ASSERT( mdl_header.GetModifiedDate() == mdl_file.GetHeader().GetModifiedDate() );
	CPPUNIT_ASSERT( mdl_header.GetNumProbes() == mdl_file.GetHeader().GetNumProbes() );
	CPPUNIT_ASSERT( mdl_header.GetVersion() == mdl_file.GetHeader().GetVersion() );
}

void CMDLFileDataTest::testproperty_FileName()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = "test";
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == mdl_file.GetFileName() );
}

////////////////////////////////////////////////
// Test file interface methods
void CMDLFileDataTest::testmethod_Exists()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = TEST_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( mdl_file.Exists() == true );
}

void CMDLFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = "test";
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( mdl_file.Exists() == false );
}

void CMDLFileDataTest::testmethod_Read()
{
	affxmdl::CMDLFileData mdl_file;

	std::string path = TEST_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", mdl_file.Exists() == true );
	CPPUNIT_ASSERT( mdl_file.Read() );
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetAnalysisType() == TEST_MDL_ANALYSIS_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetChipType() == TEST_MDL_CHIP_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetNumProbes() == TEST_MDL_NUM_PROBE_COUNT);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetVersion() == TEST_MDL_VERSION);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetModifiedDate() == TEST_MDL_MODIFIED_DATE);
}

void CMDLFileDataTest::testmethod_ReadHeader()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = TEST_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", mdl_file.Exists() == true );
	CPPUNIT_ASSERT( mdl_file.ReadHeader() );
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetAnalysisType() == TEST_MDL_ANALYSIS_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetChipType() == TEST_MDL_CHIP_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetNumProbes() == TEST_MDL_NUM_PROBE_COUNT);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetVersion() == TEST_MDL_VERSION);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetModifiedDate() == TEST_MDL_MODIFIED_DATE);
}

void CMDLFileDataTest::testmethod_ReadWhenNoDataIsInTheFile()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = NODATA_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", mdl_file.Exists() == true );
	mdl_file.Read();
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetNumProbes() == 0 );
}

void CMDLFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxmdl::CMDLFileData mdl_file;
	std::string path = "test";
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( mdl_file.Read() == false );
}

void CMDLFileDataTest::testmethod_InitializeForWriting()
{
	affxmdl::CMDLFileData mdl_file;
	affxmdl::CMDLProbeData mdl_data;

	mdl_file.InitializeForWriting(TEST_MDL_NUM_PROBE_COUNT, TEST_MDL_AFFINITY_ONE, TEST_MDL_ANALYSIS_TYPE);

	CPPUNIT_ASSERT( mdl_file.GetHeader().GetNumProbes() == TEST_MDL_NUM_PROBE_COUNT );
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetAnalysisType() == TEST_MDL_ANALYSIS_TYPE );

	// Get lower bound probe data
	mdl_file.GetProbeData(0, mdl_data);
	CPPUNIT_ASSERT( mdl_data.GetAffinity() == TEST_MDL_AFFINITY_ONE );

	// Get upper bound probe data
	mdl_file.GetProbeData(TEST_MDL_NUM_PROBE_COUNT-1, mdl_data);
	CPPUNIT_ASSERT( mdl_data.GetAffinity() == TEST_MDL_AFFINITY_ONE );
}

void CMDLFileDataTest::testmethod_SetProbeSetData()
{
	affxmdl::CMDLFileData mdl_file;
	affxmdl::CMDLProbeData mdl_data_get;
	affxmdl::CMDLProbeData mdl_data_set;

	mdl_file.InitializeForWriting(1, TEST_MDL_AFFINITY_ONE, TEST_MDL_ANALYSIS_TYPE);

	double value = TEST_MDL_AFFINITY_ONE;
	mdl_data_set.SetAffinity(value);
	mdl_file.SetProbeData(0, mdl_data_set);

	mdl_file.GetProbeData(0, mdl_data_get);
	CPPUNIT_ASSERT( mdl_data_get.GetAffinity() == TEST_MDL_AFFINITY_ONE );
}

void CMDLFileDataTest::testmethod_GetProbeSetData()
{
	affxmdl::CMDLFileData mdl_file;
	affxmdl::CMDLProbeData mdl_data;

	std::string path = WRITE_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", mdl_file.Exists() == true );
	mdl_file.Read();

	CPPUNIT_ASSERT( mdl_file.GetHeader().GetAnalysisType() == TEST_MDL_ANALYSIS_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetChipType() == TEST_MDL_CHIP_TYPE);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetNumProbes() == TEST_MDL_NUM_PROBE_COUNT);
	CPPUNIT_ASSERT( mdl_file.GetHeader().GetVersion() == TEST_MDL_VERSION);

	for (int i=0; i<mdl_file.GetHeader().GetNumProbes(); ++i)
	{
		mdl_file.GetProbeData(i, mdl_data);
		CPPUNIT_ASSERT_DOUBLES_EQUAL( mdl_data.GetAffinity(), TEST_MDL_AFFINITY_TWO, 0.1);
	}
}

void CMDLFileDataTest::testmethod_Write()
{
	affxmdl::CMDLFileData mdl_file;
	affxmdl::CMDLFileHeader mdl_header;
	affxmdl::CMDLProbeData mdl_data;

	mdl_header.SetAnalysisType(TEST_MDL_ANALYSIS_TYPE);
	mdl_header.SetChipType(TEST_MDL_CHIP_TYPE);
	mdl_header.SetVersion(TEST_MDL_VERSION);
	mdl_header.SetModifiedDate(TEST_MDL_MODIFIED_DATE);
	mdl_header.SetNumProbes(TEST_MDL_NUM_PROBE_COUNT);

	std::string path = WRITE_MDL_FILE_NAME;
	mdl_file.SetFileName(path.c_str());
	mdl_file.SetHeader(mdl_header);

	mdl_file.InitializeForWriting(TEST_MDL_NUM_PROBE_COUNT, TEST_MDL_AFFINITY_ONE, TEST_MDL_ANALYSIS_TYPE);
	for (int i=0; i<mdl_file.GetHeader().GetNumProbes(); ++i)
	{
		double value = TEST_MDL_AFFINITY_TWO;
		mdl_data.SetAffinity(value);
		mdl_file.SetProbeData(i, mdl_data);
	}
	CPPUNIT_ASSERT( mdl_file.Write() );
}

void CMDLFileDataTest::testSize()
{
	CPPUNIT_ASSERT( sizeof(affxmdl::MDLHeader) == 1172 );
}

