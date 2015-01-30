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

#include "exp_report/test/ExpressionControlResultTest.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionControlResultTest );

void ExpressionControlResultTest::setUp()
{
}

void ExpressionControlResultTest::tearDown()
{
}

void ExpressionControlResultTest::testClass()
{
	ExpressionControlResult stats;

	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.GetThreeFiveRatio(), -1, 0.001 );
	stats.SetThreeFiveRatio(0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.GetThreeFiveRatio(), 0.001, 0.001 );
	CPPUNIT_ASSERT( stats.GetName() == "" );
	stats.SetName("abc");
	CPPUNIT_ASSERT( stats.GetName() == "abc" );
	

	for (int i=0; i<NUM_SETS_PER_CONTROL; i++)
	{
		ExpressionControl::ControlValueType ii = (ExpressionControl::ControlValueType) i;
		CPPUNIT_ASSERT(stats.HasControlResult(ii) == false);
		CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.GetControlSignalResult(ii), -1, 0.0001 );
		CPPUNIT_ASSERT( stats.GetControlDetectionResult(ii) == ReportDataAccessor::DetectionNoCall);
		stats.SetControlSignalResult(ii, 1.0f);
		stats.SetControlDetectionResult(ii, ReportDataAccessor::DetectionPresent);
		CPPUNIT_ASSERT_DOUBLES_EQUAL( stats.GetControlSignalResult(ii), 1, 0.0001 );
		CPPUNIT_ASSERT( stats.GetControlDetectionResult(ii) == ReportDataAccessor::DetectionPresent);
		CPPUNIT_ASSERT(stats.HasControlResult(ii) == true);
	}
}
