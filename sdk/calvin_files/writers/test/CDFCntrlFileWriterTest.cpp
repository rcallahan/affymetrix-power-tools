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
#include "calvin_files/writers/test/CDFCntrlFileWriterTest.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CDFCntrlFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CDFCntrlFileWriterTest );

void CDFCntrlFileWriterTest::setUp()
{
	CDFData data("CDFCntrl_file");
	data.SetProbeSetCnt(10, Expression);
	writer = new CDFCntrlFileWriter(data);
}

void CDFCntrlFileWriterTest::tearDown()
{
	delete writer;
}

void CDFCntrlFileWriterTest::testCreation()
{
	CDFData c("CDFCntrl_file");
	c.SetProbeSetCnt(10, Expression);
	CDFCntrlFileWriter* w = new CDFCntrlFileWriter(c);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CDFCntrlFileWriterTest::WriteTest()
{
	CPPUNIT_ASSERT(1);
}
