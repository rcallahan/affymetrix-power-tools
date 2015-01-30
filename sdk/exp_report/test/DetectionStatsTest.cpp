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

#include "exp_report/test/DetectionStatsTest.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( DetectionStatsTest );

void DetectionStatsTest::setUp()
{
}

void DetectionStatsTest::tearDown()
{
}

void DetectionStatsTest::testClass()
{
	DetectionStats stats;

	CPPUNIT_ASSERT( stats.Count() == 0 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.Signal(), 0, 0.0001 );

	stats.IncrementCount();

	CPPUNIT_ASSERT( stats.Count() == 1 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.Signal(), 0, 0.0001 );

	stats.AddSignal(0.05f);

	CPPUNIT_ASSERT( stats.Count() == 1 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.Signal(), 0.05, 0.0001 );

	stats.IncrementCount();
	stats.AddSignal(1);

	CPPUNIT_ASSERT( stats.Count() == 2 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.Signal(), 1.05, 0.0001 );

	stats.Clear();

	CPPUNIT_ASSERT( stats.Count() == 0 );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.Signal(), 0, 0.0001 );




}
