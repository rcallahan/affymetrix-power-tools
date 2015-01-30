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

#include "exp_report/test/ProbeSetStatsTest.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ProbeSetStatsTest );

void ProbeSetStatsTest::setUp()
{
}

void ProbeSetStatsTest::tearDown()
{
}

void ProbeSetStatsTest::testClass()
{
	ProbeSetStats stats;
	CPPUNIT_ASSERT( stats.NumSets() == 0 );
	CPPUNIT_ASSERT( stats.MarginalCalls().Count() == 0);
	CPPUNIT_ASSERT( stats.AbsentCalls().Count() == 0);
	CPPUNIT_ASSERT( stats.PresentCalls().Count() == 0);

	stats.AddSet();
	stats.AddSet();
	CPPUNIT_ASSERT( stats.NumSets() == 2 );

	stats.Clear();
	CPPUNIT_ASSERT( stats.NumSets() == 0 );


	stats.PresentCalls().IncrementCount();
	CPPUNIT_ASSERT( stats.PresentCalls().Count() == 1);

	stats.AbsentCalls().IncrementCount();
	CPPUNIT_ASSERT( stats.AbsentCalls().Count() == 1);

	stats.MarginalCalls().IncrementCount();
	CPPUNIT_ASSERT( stats.MarginalCalls().Count() == 1);

	stats.Clear();
	CPPUNIT_ASSERT( stats.NumSets() == 0 );
	CPPUNIT_ASSERT( stats.MarginalCalls().Count() == 0);
	CPPUNIT_ASSERT( stats.AbsentCalls().Count() == 0);
	CPPUNIT_ASSERT( stats.PresentCalls().Count() == 0);

}
