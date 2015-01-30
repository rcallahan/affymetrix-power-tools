////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#include "calvin_files/fusion/src/FusionCHPMultiDataAccessor.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace std;
using namespace affymetrix_fusion_io;

class FusionCHPMultiDataAccessorTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( FusionCHPMultiDataAccessorTest );

	CPPUNIT_TEST ( testRead );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void testRead();
};

void FusionCHPMultiDataAccessorTest::setUp()
{
}

void FusionCHPMultiDataAccessorTest::tearDown()
{
}

void FusionCHPMultiDataAccessorTest::testRead()
{
    vector<string> chps;
    chps.push_back("../../parsers/data/CHP_MultiData_file");
    FusionCHPMultiDataAccessor accessor;
    CPPUNIT_ASSERT(accessor.Initialize(chps) == true);

    vector<string> snps;
    snps.push_back("gn1");
    snps.push_back("gn2");
    vector<vector<u_int8_t> > calls;
    vector<vector<float> > confs;
    accessor.ExtractData(snps, calls, confs);

    CPPUNIT_ASSERT(calls.size() == 1);
    CPPUNIT_ASSERT(confs.size() == 1);

    CPPUNIT_ASSERT(calls[0].size() == 2);
    CPPUNIT_ASSERT(confs[0].size() == 2);

    CPPUNIT_ASSERT(calls[0][0] == 1);
    CPPUNIT_ASSERT(calls[0][1] == 2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(confs[0][0], 11.0, 0.000001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(confs[0][1], 22.0, 0.000001f);
}
