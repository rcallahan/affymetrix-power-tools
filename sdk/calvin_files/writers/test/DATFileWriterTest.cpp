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

#include "calvin_files/writers/test/DATFileWriterTest.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/DATFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DATFileWriterTest );

void DATFileWriterTest::setUp() {}

void DATFileWriterTest::tearDown() {}

void DATFileWriterTest::testCreation()
{
	DATData fHdr("DAT_file");
	DATFileWriter* w = new DATFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void DATFileWriterTest::WriteTest()
{
	DATData data("DAT_file");
	data.SetPixelCount(10);
	data.SetStatsCount(10);
	DATFileWriter* writer = new DATFileWriter(data);
	u_int16_t stat1 = 16;
	u_int16_t stat2 = 22;
	Uint16Vector stats;
	stats.push_back(stat1);
	stats.push_back(stat2);
	writer->WriteStats(stats);

	u_int16_t p1 = 36;
	u_int16_t p2 = 3;
	Uint16Vector pixels;
	pixels.push_back(p1);
	pixels.push_back(p2);
	writer->WritePixels(pixels);
	CPPUNIT_ASSERT(1);
	delete writer;
}
