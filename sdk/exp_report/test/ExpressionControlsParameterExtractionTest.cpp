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

#include "exp_report/src/ExpressionControlsParameterExtraction.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class ExpressionControlsParameterExtractionTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ExpressionControlsParameterExtractionTest );

	CPPUNIT_TEST( testExtractParameters );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
    void testExtractParameters();
};

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionControlsParameterExtractionTest );

void ExpressionControlsParameterExtractionTest::setUp()
{
}

void ExpressionControlsParameterExtractionTest::tearDown()
{
}

void ExpressionControlsParameterExtractionTest::testExtractParameters()
{
    ExpressionReport::ExpressionControls controls;
    int thr=0;

    CPPUNIT_ASSERT(ExpressionControlsParameterExtraction::ExtractParameters("nofile", thr, controls) == false);
    CPPUNIT_ASSERT(ExpressionControlsParameterExtraction::ExtractParameters("../data/Hu6800.default.report_controls", thr, controls) == true);

    CPPUNIT_ASSERT(thr == 10);

    CPPUNIT_ASSERT(controls.ProbeArrayType() == "Hu6800");

    ExpressionReport::ExpressionControlList::iterator it = controls.SpikeControls().begin();

    CPPUNIT_ASSERT(controls.SpikeControls().size() == 1);

    CPPUNIT_ASSERT(it->Name() == "BioB");
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::MIDDLE_PROBE_SET) == false);

    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == 1);
    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == 3);



    it = controls.HousekeepingControls().begin();

    CPPUNIT_ASSERT(controls.HousekeepingControls().size() == 2);

    CPPUNIT_ASSERT(it->Name() == "GAPDH");
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::MIDDLE_PROBE_SET) == true);

    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == 11);
    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == 13);
    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::MIDDLE_PROBE_SET) == 12);


    ++it;
    CPPUNIT_ASSERT(it->Name() == "Lys");
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == true);
    CPPUNIT_ASSERT(it->HasValue(ExpressionReport::ExpressionControl::MIDDLE_PROBE_SET) == true);

    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::FIVE_PRIME_PROBE_SET) == 310);
    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::THREE_PRIME_PROBE_SET) == 330);
    CPPUNIT_ASSERT(it->GetProbeSetIndex(ExpressionReport::ExpressionControl::MIDDLE_PROBE_SET) == 320);

}

