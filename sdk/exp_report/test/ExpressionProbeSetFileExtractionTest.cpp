////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "exp_report/src/ExpressionProbeSetFileExtraction.h"
//
#include "util/Err.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class ExpressionProbeSetFileExtractionTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ExpressionProbeSetFileExtractionTest );

	CPPUNIT_TEST( testExtractParameters );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
    void testExtractParameters();
};

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionProbeSetFileExtractionTest );

void ExpressionProbeSetFileExtractionTest::setUp()
{
}

void ExpressionProbeSetFileExtractionTest::tearDown()
{
}

void ExpressionProbeSetFileExtractionTest::testExtractParameters()
{
    ProbeSetFileEntryMap entries;

    Err::setThrowStatus(true);    
    //CPPUNIT_ASSERT(ExpressionProbeSetFileExtraction::ExtractParameters("../data/nofile.txt", entries) == false);
    CPPUNIT_ASSERT(ExpressionProbeSetFileExtraction::ExtractParameters("../data/probesetlist.txt", entries) == true);
    CPPUNIT_ASSERT(entries.size() == 3);
    CPPUNIT_ASSERT(entries["10"] == "AFFX-BioB-3_at");
    CPPUNIT_ASSERT(entries["40"] == "AFFX-BioB-5_at");
    CPPUNIT_ASSERT(entries["60"] == "AFFX-BioB-M_at");

}

