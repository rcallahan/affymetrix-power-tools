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
#include "calvin_files/writers/test/CDFCntrlProbeSetWriterTest.h"
//
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/writers/src/CDFCntrlProbeSetWriter.h"
//
#include "util/Fs.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CDFCntrlProbeSetWriterTest );

void CDFCntrlProbeSetWriterTest::setUp()
{

	std::string f = "cdf_cntrl_probeset_writer";
        Fs::aptOpen(os, f, std::ios::out|std::ios::binary|std::ios::trunc);
	if (!os.is_open())
	{
	CPPUNIT_ASSERT(0);
	}
	else 
	{
		hdr.SetName(L"xyz123");
		hdr.AddAsciiColumn(L"", 64);
		hdr.AddIntColumn(L"");
		hdr.SetRowCnt(3);
		probeSetWriter = new CDFCntrlProbeSetWriter(new DataSetWriter(&os, &hdr));
	}
}

void CDFCntrlProbeSetWriterTest::tearDown()
{
	delete probeSetWriter;
	os.close();
	hdr.Clear();
}

void CDFCntrlProbeSetWriterTest::testCreation()
{
	CDFCntrlProbeSetWriter w(new DataSetWriter(&os, &hdr));
	CPPUNIT_ASSERT(1);
}

void CDFCntrlProbeSetWriterTest::WriteTest()
{
	probeSetWriter->WriteHeader();
	CPPUNIT_ASSERT(1);

	probeSetWriter->Write(10, 10, 3, 56, 7);
	CPPUNIT_ASSERT(1);
}
