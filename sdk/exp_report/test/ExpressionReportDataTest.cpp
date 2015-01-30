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

#include "exp_report/test/ExpressionReportDataTest.h"
//
#include "exp_report/src/ExpressionReportData.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionReportDataTest );

void ExpressionReportDataTest::setUp()
{
}

void ExpressionReportDataTest::tearDown()
{
}

void ExpressionReportDataTest::testClear()
{
}

void ExpressionReportDataTest::testProperties()
{
	ExpressionReportData data;

	data.ArrayType() = "test3";
	CPPUNIT_ASSERT( data.ArrayType() == "test3" );

	data.AlgName() = "alg";
	CPPUNIT_ASSERT( data.AlgName() == "alg" );

	data.CHPFileName() = "f";
	CPPUNIT_ASSERT( data.CHPFileName() == "f" );

	data.Date() = "Date";
	CPPUNIT_ASSERT( data.Date() == "Date" );

	data.ProbePairThreshold() = 10;
	CPPUNIT_ASSERT( data.ProbePairThreshold() == 10 );

	data.AntiSenseControls() = false;
	CPPUNIT_ASSERT( data.AntiSenseControls() == false);

	data.AntiSenseControls() = true;
	CPPUNIT_ASSERT( data.AntiSenseControls() == true);

	data.BackgroundStats().avg = 11.0f;
	data.BackgroundStats().max = 21.0f;
	data.BackgroundStats().min = 31.0f;
	data.BackgroundStats().std = 41.0f;

	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().avg, 11.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().max, 21.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().min, 31.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().std, 41.0f, 0.0001f);

	data.NoiseStats().avg = 111.0f;
	data.NoiseStats().max = 121.0f;
	data.NoiseStats().min = 131.0f;
	data.NoiseStats().std = 141.0f;

	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().avg, 111.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().max, 121.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().min, 131.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().std, 141.0f, 0.0001f);

	NameValuePair p;
	p.name = "n1";
	p.value = "v1";
	data.AlgParams().push_back(p);
	p.name = "n2";
	p.value = "v2";
	data.AlgParams().push_back(p);

	CPPUNIT_ASSERT(data.AlgParams().size() == 2);

	NameValuePairList::iterator it = data.AlgParams().begin();
	CPPUNIT_ASSERT(it->name == "n1");
	CPPUNIT_ASSERT(it->value == "v1");
	++it;
	CPPUNIT_ASSERT(it->name == "n2");
	CPPUNIT_ASSERT(it->value == "v2");

}
void ExpressionReportDataTest::testNameValuePair()
{
	NameValuePair n;
	n.name = "n";
	n.value = "v";
	n.Clear();
	CPPUNIT_ASSERT(n.name == "");
	CPPUNIT_ASSERT(n.value == "");
}

void ExpressionReportDataTest::testAvgStdvMinMax()
{
	AvgStdvMinMax a;
	a.avg = 1.0f;
	a.std = 2.0f;
	a.min = 3.0f;
	a.max = 4.0f;
	a.Clear();
	CPPUNIT_ASSERT(a.avg == 0.0f);
	CPPUNIT_ASSERT(a.std == 0.0f);
	CPPUNIT_ASSERT(a.min == 0.0f);
	CPPUNIT_ASSERT(a.max == 0.0f);
}

void ExpressionReportDataTest::testNameAvgCount()
{
	NameAvgCount n;
	n.name = "n";
	n.avg = 1.0f;
	n.count = 10;
	n.Clear();
	CPPUNIT_ASSERT(n.name == "");
	CPPUNIT_ASSERT(n.avg == 0.0f);
	CPPUNIT_ASSERT(n.count == 0);
}

