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

#include "rawq/test/RawQTest.h"
//
#include "rawq/src/RawQ.h"
//
#include <cstring>
#include <string>
//

CPPUNIT_TEST_SUITE_REGISTRATION( RawQTest );

void RawQTest::setUp()
{
}

void RawQTest::tearDown()
{
}

void RawQTest::testProperties()
{
	const double eps = 1e-10;
	CRawQ r;

	CPPUNIT_ASSERT(r.GetVerticalZones() == 4);
	CPPUNIT_ASSERT(r.GetHorizontalZones() == 4);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetPercentBG(), 2.0, eps);

	r.SetVerticalZones(1);
	r.SetHorizontalZones(2);
	r.SetPercentBG(3.0f);

	CPPUNIT_ASSERT(r.GetVerticalZones() == 1);
	CPPUNIT_ASSERT(r.GetHorizontalZones() == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetPercentBG(), 3.0, eps);

	r.SetDefaults();

	CPPUNIT_ASSERT(r.GetVerticalZones() == 4);
	CPPUNIT_ASSERT(r.GetHorizontalZones() == 4);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetPercentBG(), 2.0, eps);
}

void RawQTest::testCompute()
{
	const double eps = 1e-6;
	CRawQ r;
	affymetrix_fusion_io::FusionCELData cel;
	cel.SetFileName("./data/Test3.CEL");
	CPPUNIT_ASSERT(cel.Read() == true);

	affymetrix_fusion_io::FusionCDFData cdf;
	cdf.SetFileName("./data/Test3.CDF");
	CPPUNIT_ASSERT(cdf.Read() == true);

	float rq = r.ComputeRawQ(cel, cdf);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(rq, 6.3773966, eps);
}
