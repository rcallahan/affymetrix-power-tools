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

#ifdef WIN32
#include <windows.h>
#endif

//
#include "calvin_files/converters/chp/test/CHPFileConverterTest.h"
//
#include "calvin_files/converters/chp/src/CHPFileConverter.h"
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
//
#include <cstdio>
#include <fstream>
//


using namespace std;
using namespace affymetrix_chp_converter;
using namespace affymetrix_calvin_io;

#ifndef WIN32
#define FALSE 0
#define TRUE 1
void CopyFile(const char *in, const char *out, bool foo) {
	std::ifstream ins(in);
	std::ofstream outs(out);
	outs << ins.rdbuf();
	ins.close();
	outs.close();
}
#endif

CPPUNIT_TEST_SUITE_REGISTRATION( CHPFileConverterTest );

void CHPFileConverterTest::setUp()
{
}

void CHPFileConverterTest::tearDown()
{
}

void CHPFileConverterTest::test_ConvertFile_to_Calvin()
{
	const char *cfile = "../data/Calvin-stat.CHP";
	const char *xfile = "../data/XDA-stat.CHP";
	//const char *origcfile = "../data/orig-Calvin-stat.CHP";
	const char *origxfile = "../data/orig-XDA-stat.CHP";
	const char *infile = "../data/testCalvin.CHP";
	const char *libPath = "../data/";
	ConverterRemoveFile(cfile);
	ConverterRemoveFile(xfile);
	ConverterRemoveFile(infile);
	CopyFile(origxfile, infile, FALSE);
	CHPFileConverter c;

	ConverterMoveFile("../data/HG-U133A.CDF", "../data/HG-U133A.CDF.txt");
	ConverterMoveFile("../data/HG-U133A.PSI", "../data/HG-U133A.PSI.txt");
	CPPUNIT_ASSERT(c.ConvertFile(infile, libPath, Calvin_Version1) == false);
	CPPUNIT_ASSERT(c.ErrorCode() == UnableToLoadProbeSetNames);

	ConverterMoveFile("../data/HG-U133A.CDF.txt", "../data/HG-U133A.CDF");
	ConverterMoveFile("../data/HG-U133A.PSI.txt", "../data/HG-U133A.PSI");
	CPPUNIT_ASSERT(c.ConvertFile(infile, libPath, Calvin_Version1) == true);

}

void CHPFileConverterTest::test_Convert_Mas5_File_to_Calvin()
{
	const char *bakFile = "../data/mas5/test3/Test3_Exp2.CHP.bak";
	const char *inFile = "../data/mas5/test3/Test3_Exp2.CHP";
	const char *libPath = "../data/mas5/test3/library/";
	const char *celFile = "../data/mas5/test3/Test3_Exp2.CEL";
	//const char *maskPath = "../data/maskPath";
	CHPFileConverter c;
	CPPUNIT_ASSERT(c.ConvertFile(inFile, libPath, celFile, NULL, Calvin_Version1, false));
	
	CHPFileReader* reader = new CHPFileReader();
	string file(inFile);
	reader->SetFilename(file);
	CHPData* data = new CHPData();
	reader->Read(*data);

	int32_t bgCnt = data->GetBackgroundZoneCnt();
	CHPBackgroundZoneVector zones;
	data->GetBackgroundZones(0, bgCnt, zones);
	size_t zoneCnt = zones.size();
	CPPUNIT_ASSERT(zoneCnt == 16);
	CPPUNIT_ASSERT(zones[0].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[0].GetBackground(), 195.582, .001);
	CPPUNIT_ASSERT(zones[1].GetCenterX() == 48);// && zones[1].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[1].GetBackground(), 180.391, .001);
	CPPUNIT_ASSERT(zones[2].GetCenterX() == 80);// && zones[2].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[2].GetBackground(), 174.345, .001);
	CPPUNIT_ASSERT(zones[3].GetCenterX() == 112);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[3].GetBackground(), 187.06, .001);
	CPPUNIT_ASSERT(zones[4].GetCenterX() == 16);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[4].GetBackground(), 201.033, .001);
	CPPUNIT_ASSERT(zones[5].GetCenterX() == 48);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[5].GetBackground(), 197.856, .001);
	CPPUNIT_ASSERT(zones[6].GetCenterX() == 80);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[6].GetBackground(), 196.96, .001);
	CPPUNIT_ASSERT(zones[7].GetCenterX() == 112);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[7].GetBackground(), 206.867, .001);
	CPPUNIT_ASSERT(zones[8].GetCenterX() == 16);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[8].GetBackground(), 235.307, .001);
	CPPUNIT_ASSERT(zones[9].GetCenterX() == 48);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[9].GetBackground(), 236.162, .001);
	CPPUNIT_ASSERT(zones[10].GetCenterX() == 80);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[10].GetBackground(), 232.188, .001);
	CPPUNIT_ASSERT(zones[11].GetCenterX() == 112);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[11].GetBackground(), 241.533, .001);
	CPPUNIT_ASSERT(zones[12].GetCenterX() == 16);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[12].GetBackground(), 235.535, .001);
	CPPUNIT_ASSERT(zones[13].GetCenterX() == 48);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[13].GetBackground(), 235.582, .001);
	CPPUNIT_ASSERT(zones[14].GetCenterX() == 80);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[14].GetBackground(), 237.438, .001);
	CPPUNIT_ASSERT(zones[15].GetCenterX() == 112);// && zones[3].GetCenterX() == 16);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(zones[15].GetBackground(), 244.48, .001);

	delete data;
	delete reader;
	ConverterMoveFile(bakFile, inFile, true);
}

void CHPFileConverterTest::test_ConvertFile_to_XDA()
{
	const char *cfile = "../data/Calvin-stat.CHP";
	const char *xfile = "../data/XDA-stat.CHP";
	const char *origcfile = "../data/orig-Calvin-stat.CHP";
	//const char *origxfile = "../data/orig-XDA-stat.CHP";
	const char *infile = "../data/testXDA.CHP";
	const char *libPath = "../data/";
	ConverterRemoveFile(cfile);
	ConverterRemoveFile(xfile);
	ConverterRemoveFile(infile);
	CopyFile(origcfile, infile, FALSE);
	CHPFileConverter c;
	CPPUNIT_ASSERT(c.ConvertFile(infile, libPath, GCOS_XDA_Version) == true);
}

void CHPFileConverterTest::test_ConvertFile_fail_to_open_file()
{
	CHPFileConverter c;
	CPPUNIT_ASSERT(c.ConvertFile("../data/no_file", "./", GCOS_XDA_Version) == false);
	CPPUNIT_ASSERT(c.ErrorCode() == FileDoesNotExist);

	CPPUNIT_ASSERT(c.ConvertFile("../data/invalid_file", "../data/", GCOS_XDA_Version) == false);
	CPPUNIT_ASSERT(c.ErrorCode() == InvalidChpFileFormat);
}

void CHPFileConverterTest::test_ConvertFile_already_in_format()
{
	const char *cfile = "../data/Calvin-stat.CHP";
	const char *xfile = "../data/XDA-stat.CHP";
	const char *origcfile = "../data/orig-Calvin-stat.CHP";
	const char *origxfile = "../data/orig-XDA-stat.CHP";
	const char *infile = "../data/test.CHP";
	const char *libPath = "../data/";
	ConverterRemoveFile(cfile);
	ConverterRemoveFile(xfile);
	ConverterRemoveFile(infile);
	CopyFile(origxfile, infile, FALSE);

	CHPFileConverter c;
	CPPUNIT_ASSERT(c.ConvertFile(infile, libPath, GCOS_XDA_Version) == true);

	ConverterRemoveFile(infile);
	CopyFile(origcfile, infile, FALSE);

	CPPUNIT_ASSERT(c.ConvertFile(infile, libPath, Calvin_Version1) == true);
}
