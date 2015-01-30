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


#include "file/CPPTest/CMSFileDataTest.h"
//
#include "file/CMSFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CCMSFileDataTest );

/* Test.CMS contents
[HEADER]
Version=1
Assay=Mapping
ArrayCount=3
SNPCount=10

[DATA]
Types=Mapping10K_Xba131	Mapping50K_Hind240	Mapping50K_Xba240	
SNP1=151	151	151
SNP2=152	152	152
SNP3=153	153	153
SNP4=154	154	154
SNP5=155	155	155
SNP6=156	156	156
SNP7=157	157	157
SNP8=158	158	158
SNP9=159	159	159
SNP10=160	160	160
*/
#define TEST_CMS_VERSION 1
#define TEST_CMS_ASSAY "Mapping"
#define TEST_CMS_ARRAY_COUNT 3
#define TEST_CMS_SNP_COUNT 31
#define TEST_CMS_ARRAY_10K_XBA131 "Mapping10K_Xba131"
#define TEST_CMS_ARRAY_50K_HIND240 "Mapping50K_Hind240"
#define TEST_CMS_ARRAY_50K_XBA240 "Mapping50K_Xba240"

#ifdef _MSC_VER
#define TEST_CMS_FILE_NAME ".\\data\\test.cms"
#define NODATA_CMS_FILE_NAME ".\\data\\nodata.cms"
#else
#define TEST_CMS_FILE_NAME "./data/test.cms"
#define NODATA_CMS_FILE_NAME "./data/nodata.cms"
#endif

void CCMSFileDataTest::setUp()
{
}

void CCMSFileDataTest::tearDown()
{
}

////////////////////////////////////////////////
// Test interface object creation
void CCMSFileDataTest::testCreation()
{
	affxcms::CCMSFileData cms_file;
	CPPUNIT_ASSERT( 1 );
}

void CCMSFileDataTest::testCreationHeader()
{
	affxcms::CCMSFileHeader cms_header;
	CPPUNIT_ASSERT( 1 );
}

////////////////////////////////////////////////
// Test header interface properties
void CCMSFileDataTest::testproperty_HeaderVersion()
{
	affxcms::CCMSFileHeader cms_header;
	int version = TEST_CMS_VERSION;
	cms_header.SetVersion(version);
	CPPUNIT_ASSERT( version == cms_header.GetVersion() );
}

void CCMSFileDataTest::testproperty_HeaderAssay()
{
	affxcms::CCMSFileHeader cms_header;
	std::string assay = TEST_CMS_ASSAY;
	cms_header.SetAssay(assay.c_str());
	CPPUNIT_ASSERT( assay == cms_header.GetAssay() );
}

void CCMSFileDataTest::testproperty_HeaderArrayCount()
{
	affxcms::CCMSFileHeader cms_header;
	int value = TEST_CMS_ARRAY_COUNT;
	cms_header.SetArrayCount(value);
	CPPUNIT_ASSERT( value == cms_header.GetArrayCount() );
}

void CCMSFileDataTest::testproperty_HeaderSNPCount()
{
	affxcms::CCMSFileHeader cms_header;
	int value = TEST_CMS_SNP_COUNT;
	cms_header.SetSNPCount(value);
	CPPUNIT_ASSERT( value == cms_header.GetSNPCount() );
}

////////////////////////////////////////////////
// Test file interface properties
void CCMSFileDataTest::testproperty_FileName()
{
	affxcms::CCMSFileData cms_file;
	std::string path = "test";
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == cms_file.GetFileName() );
}

////////////////////////////////////////////////
// Test file interface methods
void CCMSFileDataTest::testmethod_Exists()
{
	affxcms::CCMSFileData cms_file;
	std::string path = TEST_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( cms_file.Exists() == true );
}

void CCMSFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxcms::CCMSFileData cms_file;
	std::string path = "test";
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( cms_file.Exists() == false );
}

void CCMSFileDataTest::testmethod_Read()
{
	affxcms::CCMSFileData cms_file;

	std::string path = TEST_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", cms_file.Exists() == true );
	CPPUNIT_ASSERT( cms_file.Read() );
}

void CCMSFileDataTest::testmethod_ReadHeader()
{
	affxcms::CCMSFileData cms_file;
	std::string path = TEST_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", cms_file.Exists() == true );
	CPPUNIT_ASSERT( cms_file.ReadHeader() );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetVersion() == TEST_CMS_VERSION );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetAssay() == TEST_CMS_ASSAY );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetArrayCount() == TEST_CMS_ARRAY_COUNT );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetSNPCount() == TEST_CMS_SNP_COUNT );
}

void CCMSFileDataTest::testmethod_ReadWhenNoDataIsInTheFile()
{
	affxcms::CCMSFileData cms_file;
	std::string path = NODATA_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", cms_file.Exists() == true );
	cms_file.Read();
	CPPUNIT_ASSERT( cms_file.GetHeader().GetSNPCount() == 0 );
}

void CCMSFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxcms::CCMSFileData cms_file;
	std::string path = "test";
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT( cms_file.Read() == false );
}

void CCMSFileDataTest::testmethod_ArrayTypeInformation()
{
	affxcms::CCMSFileData cms_file;

	std::string path = TEST_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", cms_file.Exists() == true );
	cms_file.Read();
	CPPUNIT_ASSERT( cms_file.GetHeader().GetVersion() == TEST_CMS_VERSION );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetAssay() == TEST_CMS_ASSAY );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetArrayCount() == TEST_CMS_ARRAY_COUNT );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetSNPCount() == TEST_CMS_SNP_COUNT );

	std::map<int,std::string> mapArray = cms_file.ArrayTypeInformation();
	CPPUNIT_ASSERT( mapArray.size() == TEST_CMS_ARRAY_COUNT );
	CPPUNIT_ASSERT( mapArray[0] == TEST_CMS_ARRAY_10K_XBA131 );
	CPPUNIT_ASSERT( mapArray[1] == TEST_CMS_ARRAY_50K_HIND240 );
	CPPUNIT_ASSERT( mapArray[2] == TEST_CMS_ARRAY_50K_XBA240 );
}

void CCMSFileDataTest::testmethod_SNPInformation()
{
	affxcms::CCMSFileData cms_file;

	std::string path = TEST_CMS_FILE_NAME;
	cms_file.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", cms_file.Exists() == true );
	cms_file.Read();
	CPPUNIT_ASSERT( cms_file.GetHeader().GetVersion() == TEST_CMS_VERSION );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetAssay() == TEST_CMS_ASSAY );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetArrayCount() == TEST_CMS_ARRAY_COUNT );
	CPPUNIT_ASSERT( cms_file.GetHeader().GetSNPCount() == TEST_CMS_SNP_COUNT );

	affxcms::SNPTypeInfoList mapSNPs = cms_file.SNPInformation();
	CPPUNIT_ASSERT( mapSNPs.size() == TEST_CMS_SNP_COUNT );

	std::list<affxcms::SNPTypeInfo>::iterator iter = mapSNPs.begin();
	std::list<affxcms::SNPTypeInfo>::iterator end = mapSNPs.end();
	if (iter!=end)
	{
		std::string strCT = TEST_CMS_ARRAY_10K_XBA131;
		std::string snpID = iter->snpID;
		int unit = iter->chiptype_to_unit[strCT];
		CPPUNIT_ASSERT( unit == 151 );
	}

	++iter;
	if (iter!=end)
	{
		std::string strCT = TEST_CMS_ARRAY_50K_HIND240;
		std::string snpID = iter->snpID;
		int unit = iter->chiptype_to_unit[strCT];
		CPPUNIT_ASSERT( unit == 152 );
	}

	++iter;
	if (iter!=end)
	{
		std::string strCT = TEST_CMS_ARRAY_50K_XBA240;
		std::string snpID = iter->snpID;
		int unit = iter->chiptype_to_unit[strCT];
		CPPUNIT_ASSERT( unit == 153 );
	}
}






