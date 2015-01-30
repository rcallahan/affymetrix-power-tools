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

#include "file/CPPTest/SMDFileDataTest.h"
//
#include "file/SMDFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( SMDFileDataTest );

#define TEST_SMD_FILE "./data/test.smd"
#define NONEXISTANT_SMD_FILE "nonexistant.smd"

void SMDFileDataTest::setUp()
{
}

void SMDFileDataTest::tearDown()
{
}

void SMDFileDataTest::CreationTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT( 1 );
}

void SMDFileDataTest::FileNameTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT_NO_THROW(smd.SetFileName(TEST_SMD_FILE));
	CPPUNIT_ASSERT(smd.GetFileName() == TEST_SMD_FILE);
}

void SMDFileDataTest::ExistsTest()
{
	affxsmd::SMDFileData smd;
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Exists());
	smd.SetFileName(NONEXISTANT_SMD_FILE);
	CPPUNIT_ASSERT(smd.Exists() == false);
}

void SMDFileDataTest::ReadTest()
{
	affxsmd::SMDFileData smd;
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());
}

void SMDFileDataTest::ReadFailedTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT(smd.Read() == false);
	smd.SetFileName(NONEXISTANT_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read() == false);
}

void SMDFileDataTest::FrameRowsMemberTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT(smd.frameRows == 0);
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());
	CPPUNIT_ASSERT(smd.frameRows == 12);
}

void SMDFileDataTest::FrameColsMemberTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT(smd.frameRows == 0);
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());
	CPPUNIT_ASSERT(smd.frameCols == 13);
}

void SMDFileDataTest::NumFramesTest()
{
	affxsmd::SMDFileData smd;
	CPPUNIT_ASSERT(smd.NumFrames() == 0);
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());
	CPPUNIT_ASSERT(smd.NumFrames() == 156);
}

void SMDFileDataTest::GetFrameTest()
{
	affxsmd::SMDFileData smd;
	affxsmd::SMDFrame frame1 = smd.GetFrame(1);
	CPPUNIT_ASSERT(frame1.frameIdx == -1);
	CPPUNIT_ASSERT(frame1.rows == -1);
	CPPUNIT_ASSERT(frame1.cols == -1);
	CPPUNIT_ASSERT(frame1.startRow == -1);
	CPPUNIT_ASSERT(frame1.startCol == -1);

	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());

	int expectedRows[12] = {126, 124, 126, 124, 126, 124, 126, 124, 126, 124, 126, 124};
	int expectedCols[13] = {126, 126, 126, 126, 126, 82, 176, 82, 126, 126, 126, 126, 126};
	int expectedStartRow[12] = {0, 124, 246, 370, 492, 616, 738, 862, 984, 1108, 1230, 1354};
	int expectedStartCol[13] = {0, 122, 244, 366, 488, 610, 688, 860, 938, 1060, 1182, 1304, 1426};

	int frameIdx = 0;
	int rows = smd.frameRows;
	int cols = smd.frameCols;

	for (int irow = 0; irow < rows; ++irow)
	{
		for (int icol = 0; icol < cols; ++icol, ++frameIdx)
		{
			affxsmd::SMDFrame frame = smd.GetFrame(frameIdx);

			CPPUNIT_ASSERT(frame.rows == expectedRows[irow]);
			CPPUNIT_ASSERT(frame.cols == expectedCols[icol]);
			CPPUNIT_ASSERT(frame.startRow == expectedStartRow[irow]);
			CPPUNIT_ASSERT(frame.startCol == expectedStartCol[icol]);
		}
	}
}

void SMDFileDataTest::GetFrameOutOfBoundsTest()
{
	affxsmd::SMDFileData smd;
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());

	affxsmd::SMDFrame frame = smd.GetFrame(200);
	CPPUNIT_ASSERT(frame.frameIdx == -1);
	CPPUNIT_ASSERT(frame.rows == -1);
	CPPUNIT_ASSERT(frame.cols == -1);
	CPPUNIT_ASSERT(frame.startRow == -1);
	CPPUNIT_ASSERT(frame.startCol == -1);
}

void SMDFileDataTest::FailedReadClearsTest()
{
	affxsmd::SMDFileData smd;
	smd.SetFileName(TEST_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read());
	smd.SetFileName(NONEXISTANT_SMD_FILE);
	CPPUNIT_ASSERT(smd.Read()==false);
	CPPUNIT_ASSERT(smd.NumFrames() == 0);
	CPPUNIT_ASSERT(smd.frameCols == 0);
	CPPUNIT_ASSERT(smd.frameRows == 0);
}
