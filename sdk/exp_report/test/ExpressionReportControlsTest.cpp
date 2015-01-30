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

#include "exp_report/test/ExpressionReportControlsTest.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionControlsTest );

void ExpressionControlsTest::setUp()
{
}

void ExpressionControlsTest::tearDown()
{
}

void ExpressionControlsTest::testExpressionControl()
{
	ExpressionControl c;
	ExpressionControl c2;

	CPPUNIT_ASSERT( c.Name() == "" );
	c.Name() = "abc";
	CPPUNIT_ASSERT( c.Name() == "abc" );
	
	CPPUNIT_ASSERT( c == "abc" );
	CPPUNIT_ASSERT( c != "aabc" );
	c2.Name() = "aabc";
	CPPUNIT_ASSERT( c != c2 );
	c2.Name() = "abc";
	CPPUNIT_ASSERT( c == c2 );

	for (int i=0; i<NUM_SETS_PER_CONTROL; i++)
	{
		ExpressionControl::ControlValueType ii = (ExpressionControl::ControlValueType) i;
		CPPUNIT_ASSERT(c.HasValue(ii) == false);
		CPPUNIT_ASSERT(c.GetProbeSetIndex(ii) == -1);
		c.SetProbeSetIndex(ii, 1);
		CPPUNIT_ASSERT(c.GetProbeSetIndex(ii) == 1);
		CPPUNIT_ASSERT(c.HasValue(ii) == true);
	}
}

void ExpressionControlsTest::testExpressionControls()
{
	ExpressionControls c1;
	ExpressionControls c2;

	c1.ProbeArrayType() = "test3";
	CPPUNIT_ASSERT(c1.ProbeArrayType() == "test3");

	c2.ProbeArrayType() = "test4";
	CPPUNIT_ASSERT(c2.ProbeArrayType() == "test4");

	CPPUNIT_ASSERT(c1 != c2);
	c1.ProbeArrayType() = "test4";
	CPPUNIT_ASSERT(c1 == c2);

	ExpressionControls c3("test4");
	CPPUNIT_ASSERT(c3 == c2);


}
