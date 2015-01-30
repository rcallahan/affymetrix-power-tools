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

#include "file/CPPTest/GRDFileWriterTest.h"
//
#include "file/FileIO.h"
#include "file/GRDFileWriter.h"
//
#include <cmath>
#include <cstring>
#include <string.h>
//

CPPUNIT_TEST_SUITE_REGISTRATION( CGRDFileWriterTest );

#define GRD_FILE "./data/test.out.grd"

void CGRDFileWriterTest::setUp()
{
}

void CGRDFileWriterTest::tearDown()
{
}

void CGRDFileWriterTest::WriteGRDFileTest()
{
	CGRDFileWriter writer;
	writer.SetFileName(GRD_FILE);

	// Initialize the GRD data object
	writer.GetHeader().SetVersion(2.0f);
	writer.GetHeader().SetCols(2);
	writer.GetHeader().SetRows(3);
	writer.GetHeader().SetFeaturePitchX(1.2f);
	writer.GetHeader().SetFeaturePitchY(1.3f);
	writer.GetHeader().SetFeatureSetbackX(1.4f);
	writer.GetHeader().SetFeatureSetbackY(1.5f);

	StrStrMap& props = writer.GetHeader().GetParameters();
	props[SZ_PARENT_DAT_PROP_NAME] = "test.out.dat";
	props[SZ_SCAN_DATE_TIME_PROP_NAME] = "2005-12-25T11:12:13Z";
	props[SZ_SCANNER_ID_PROP_NAME] = "Main Scanner";

	FRECT r;
	r.ul.fx = 1.1f;
	r.ul.fy = 1.2f;
	r.ur.fx = 10.1f;
	r.ur.fy = 1.3f;
	r.lr.fx = 10.2f;
	r.lr.fy = 10.3f;
	r.ll.fx = 1.4f;
	r.ll.fx = 10.4f;

	writer.GetHeader().AddOptSubgrid(&r);

	// Create the new file.
	CPPUNIT_ASSERT(writer.CreateNewFile());

	for (uint32_t row = 0; row < 3; ++row)
	{
		for (uint32_t col = 0; col < 2; ++col)
		{
			FCOORD coord;
			coord.fx = col;
			coord.fy = row;
			CPPUNIT_ASSERT(writer.WriteCenter(coord));
		}
	}

	CPPUNIT_ASSERT(writer.CloseNewFile());

	// Read the GRD file and check it.
	CGRDFileData grd;
	grd.SetFileName(GRD_FILE);
	CPPUNIT_ASSERT(grd.Read());

	double eps = 0.000001;

	CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0f, grd.GetHeader().GetVersion(), eps);
	CPPUNIT_ASSERT(grd.GetHeader().GetCols() == 2);
	CPPUNIT_ASSERT(grd.GetHeader().GetRows() == 3);
	CPPUNIT_ASSERT(grd.GetHeader().GetRows()*grd.GetHeader().GetCols() == grd.GetHeader().GetNumCells());

	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2, grd.GetHeader().GetFeaturePitchX(), eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3, grd.GetHeader().GetFeaturePitchY(), eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4, grd.GetHeader().GetFeatureSetbackX(), eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, grd.GetHeader().GetFeatureSetbackY(), eps);

	StrStrMap &params = grd.GetHeader().GetParameters();
	CPPUNIT_ASSERT(params.size() == 3);
	std::string val;

	val = params[SZ_PARENT_DAT_PROP_NAME];
	CPPUNIT_ASSERT( val == "test.out.dat" );

	val = params[SZ_SCAN_DATE_TIME_PROP_NAME];
	CPPUNIT_ASSERT( val == "2005-12-25T11:12:13Z" );

	val = params[SZ_SCANNER_ID_PROP_NAME];
	CPPUNIT_ASSERT( val == "Main Scanner");


	CPPUNIT_ASSERT( grd.GetHeader().GetNumSubgrids() == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ul.fx, grd.GetHeader().GetOptSubgrid(0).ul.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ul.fy, grd.GetHeader().GetOptSubgrid(0).ul.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ur.fx, grd.GetHeader().GetOptSubgrid(0).ur.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ur.fy, grd.GetHeader().GetOptSubgrid(0).ur.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ll.fx, grd.GetHeader().GetOptSubgrid(0).ll.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ll.fy, grd.GetHeader().GetOptSubgrid(0).ll.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.lr.fx, grd.GetHeader().GetOptSubgrid(0).lr.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.lr.fy, grd.GetHeader().GetOptSubgrid(0).lr.fy, eps);

	for (uint32_t row = 0; row < 3; ++row)
	{
		for (uint32_t col = 0; col < 2; ++col)
		{
			FCOORD center = grd.GetCenter(col, row);

			CPPUNIT_ASSERT_DOUBLES_EQUAL(col, center.fx, eps);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(row, center.fy, eps);
		}
	}

	grd.Close();


	// Check GRD file formatting
	std::ifstream ifs;
	ifs.open(GRD_FILE, std::ios::out | std::ios::binary);
	CPPUNIT_ASSERT(ifs);

	// Check the file identity
	char buf[8];
	ifs.read(buf, 8);
	CPPUNIT_ASSERT(memcmp(buf, "\211GRD\r\n\032\n", 8) == 0);

	// Check fixed header
	float version = 0.0f;
	ReadFloat_N(ifs, version);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0f, version, eps);
	uint32_t rows = 0;
	ReadUInt32_N(ifs, rows);
	CPPUNIT_ASSERT(rows == 3);
	uint32_t cols = 0;
	ReadUInt32_N(ifs, cols);
	CPPUNIT_ASSERT(cols == 2);
	float featurePitchX = 0.f;
	ReadFloat_N(ifs, featurePitchX);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2, featurePitchX, eps);
	float featurePitchY = 0.f;
	ReadFloat_N(ifs, featurePitchY);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3, featurePitchY, eps);
	float featureSetbackX = 0.f;
	ReadFloat_N(ifs, featureSetbackX);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4, featureSetbackX, eps);
	float featureSetbackY = 0.f;
	ReadFloat_N(ifs, featureSetbackY);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, featureSetbackY, eps);

	// Check the tag header section
	std::streampos posStart = ifs.tellg();
	uint32_t sectionBytes = 0;
	ReadUInt32_N(ifs, sectionBytes);
	uint32_t numParams = 0;
	ReadUInt32_N(ifs, numParams);
	CPPUNIT_ASSERT(numParams == 3);

	for (uint32_t param = 0; param < numParams; ++param)
	{
		std::string name, value;
		ReadString_N(ifs, name);
		ReadString_N(ifs, value);

		if (name == SZ_PARENT_DAT_PROP_NAME)
			CPPUNIT_ASSERT(value == "test.out.dat");
		else if (name == SZ_SCAN_DATE_TIME_PROP_NAME)
			CPPUNIT_ASSERT(value == "2005-12-25T11:12:13Z");
		else if (name == SZ_SCANNER_ID_PROP_NAME)
			CPPUNIT_ASSERT(value == "Main Scanner");
		else
			CPPUNIT_ASSERT_MESSAGE("Missing parameter", 0);
	}

	uint32_t diff = ifs.tellg() - posStart;
	CPPUNIT_ASSERT(diff == sectionBytes);

	// Go back to the start of the section
	ifs.seekg(posStart);
	// and skip over it.
	ifs.seekg(sectionBytes, std::ios_base::cur);

	// Check optimum corners section
	posStart = ifs.tellg();
	ReadUInt32_N(ifs, sectionBytes);
	uint32_t numSG = 0;
	ReadUInt32_N(ifs, numSG);
	CPPUNIT_ASSERT(numSG == 1);

	FRECT r2;
	ReadFloat_N(ifs, r2.ul.fx);
	ReadFloat_N(ifs, r2.ul.fy);
	ReadFloat_N(ifs, r2.ur.fx);
	ReadFloat_N(ifs, r2.ur.fy);
	ReadFloat_N(ifs, r2.ll.fx);
	ReadFloat_N(ifs, r2.ll.fy);
	ReadFloat_N(ifs, r2.lr.fx);
	ReadFloat_N(ifs, r2.lr.fy);

	CPPUNIT_ASSERT( grd.GetHeader().GetNumSubgrids() == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ul.fx, r2.ul.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ul.fy, r2.ul.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ur.fx, r2.ur.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ur.fy, r2.ur.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ll.fx, r2.ll.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.ll.fy, r2.ll.fy, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.lr.fx, r2.lr.fx, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( r.lr.fy, r2.lr.fy, eps);

	diff = ifs.tellg() - posStart;
	CPPUNIT_ASSERT(diff == sectionBytes);

	// Go back to the start of the section
	ifs.seekg(posStart);
	// and skip over it.
	ifs.seekg(sectionBytes, std::ios_base::cur);

	// Check a few centers
	float centerX = 0.f, centerY = 0.f;
	ReadFloat_N(ifs, centerX);
	ReadFloat_N(ifs, centerY);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0, centerX, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0, centerY, eps);
	ReadFloat_N(ifs, centerX);
	ReadFloat_N(ifs, centerY);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1, centerX, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0, centerY, eps);

}
