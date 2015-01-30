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

#include "exp_report/test/ChangeStatsTest.h"
//
#include "exp_report/src/ChangeStats.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ChangeStatsTest );

void ChangeStatsTest::setUp()
{
}

void ChangeStatsTest::tearDown()
{
}

void ChangeStatsTest::testClass()
{
	ChangeStats stats;

	CPPUNIT_ASSERT( stats.ChangeCount() == 0 );
	CPPUNIT_ASSERT( stats.DetectionChangeCount()==0);
	CPPUNIT_ASSERT( stats.DetectionPresentCount()==0);
	CPPUNIT_ASSERT( stats.DetectionAbsentCount()==0);
	CPPUNIT_ASSERT( stats.ModerateCount()==0);
	for (int i=0; i<NUMBER_FOLD_CHANGE_BINS; i++)
		CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.FoldChangeCount(i), 0, 0.001);

	float l,u;

	ChangeStats::GetBinRange(0, l, u);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(l, 0, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(u, 1, 0.001);

	ChangeStats::GetBinRange(1, l, u);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(l, 1, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(u, 2, 0.001);

	ChangeStats::GetBinRange(2, l, u);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(l, 2, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(u, 3, 0.001);

	ChangeStats::GetBinRange(3, l, u);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(l, 3, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(u, 4, 0.001);

	ChangeStats::GetBinRange(4, l, u);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(l, 4, 0.001);
	CPPUNIT_ASSERT(u > 99999999.0f);



	stats.IncrementChangeCount();
	CPPUNIT_ASSERT(stats.ChangeCount() == 1);

	stats.IncrementDetectionChangeCount();
	CPPUNIT_ASSERT(stats.DetectionChangeCount() == 1);

	stats.IncrementDetectionPresentCount();
	CPPUNIT_ASSERT(stats.DetectionPresentCount() == 1);

	stats.IncrementDetectionAbsentCount();
	CPPUNIT_ASSERT(stats.DetectionAbsentCount() == 1);

	stats.IncrementModerateCount();
	CPPUNIT_ASSERT(stats.ModerateCount() == 1);

	stats.IncrementFoldChangeCount(0);
	CPPUNIT_ASSERT(stats.FoldChangeCount(0) == 1);

	stats.IncrementFoldChangeCount(1);
	CPPUNIT_ASSERT(stats.FoldChangeCount(1) == 1);

	stats.IncrementFoldChangeCount(2);
	CPPUNIT_ASSERT(stats.FoldChangeCount(2) == 1);

	stats.IncrementFoldChangeCount(3);
	CPPUNIT_ASSERT(stats.FoldChangeCount(3) == 1);

	stats.IncrementFoldChangeCount(4);
	CPPUNIT_ASSERT(stats.FoldChangeCount(4) == 1);




}
