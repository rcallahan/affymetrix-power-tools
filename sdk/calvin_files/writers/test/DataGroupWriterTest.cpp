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

#include "calvin_files/writers/test/DataGroupWriterTest.h"
//
#include "calvin_files/writers/src/DataGroupWriter.h"
//
#include "util/Fs.h"
//
using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupWriterTest );

void DataGroupWriterTest::setUp()
{
    std::string f = "data dataGroup writer";
    Fs::aptOpen(os, f, std::ios::out|std::ios::binary|std::ios::trunc);
    if (!os.is_open())
	{
		CPPUNIT_ASSERT(0);
	}
    else 
    {
        hdr = new DataGroupHeader(L"data dataGroup writer");
        writer = new DataGroupWriter(&os, hdr);
    }
}

void DataGroupWriterTest::tearDown()
{
    delete writer;
    delete hdr;
}

void DataGroupWriterTest::testCreation()
{
    std::ofstream os;
    DataGroupHeader dcHdr(L"data dataGroup writer");
    DataGroupWriter w(&os, &dcHdr);
	CPPUNIT_ASSERT(1);
}

void DataGroupWriterTest::WriteTest()
{
	writer->WriteHeader();
    CPPUNIT_ASSERT(1);
}
