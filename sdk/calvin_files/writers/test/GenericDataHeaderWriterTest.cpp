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

#include "calvin_files/writers/test/GenericDataHeaderWriterTest.h"
//
#include "calvin_files/writers/src/GenericDataHeaderWriter.h"
//
#include "util/Fs.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataHeaderWriterTest );

void GenericDataHeaderWriterTest::setUp()
{
    writer = new GenericDataHeaderWriter();
}

void GenericDataHeaderWriterTest::tearDown()
{
    delete writer;
}

void GenericDataHeaderWriterTest::testCreation()
{
	GenericDataHeaderWriter w;
	CPPUNIT_ASSERT(1);
}

void GenericDataHeaderWriterTest::WriteTest()
{
	std::ofstream os;
    std::string f = "generic_data_header";
    Fs::aptOpen(os, f, std::ios::out|std::ios::binary|std::ios::trunc);
	if (!os)
	{
		CPPUNIT_ASSERT(0);
	}
    else
    {
	    GenericDataHeader g;
        g.SetFileId("mimi rogers");
	    writer->Write(os, g);
        CPPUNIT_ASSERT(1);
        os.close();
    }

    //TODO: use generic data header reader to read header and verify data
}
