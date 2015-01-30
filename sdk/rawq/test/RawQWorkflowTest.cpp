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

#include "rawq/test/RawQWorkflowTest.h"
//
#include "rawq/src/RawQWorkflow.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( RawQWorkflowTest );

void RawQWorkflowTest::setUp()
{
}

void RawQWorkflowTest::tearDown()
{
}

void RawQWorkflowTest::testProperties()
{
	const double eps = 1e-10;
	CRawQWorkflow r;

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

void RawQWorkflowTest::testCompute()
{
	const double eps = 1e-6;
	CRawQWorkflow r;
	CPPUNIT_ASSERT(r.ComputeRawQ("./data/Test3.CEL", "./data") == CRawQWorkflow::NoError);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetRawQ(), 6.3773966, eps);
}

void RawQWorkflowTest::testCompute_no_CEL_file()
{
	CRawQWorkflow r;
	CPPUNIT_ASSERT(r.ComputeRawQ("./data/no_file.CEL", "./data") == CRawQWorkflow::UnableToReadCELFile);
}

void RawQWorkflowTest::testCompute_no_CDF_file()
{
	CRawQWorkflow r;
	CPPUNIT_ASSERT(r.ComputeRawQ("./data/Test3.CEL", "./no_data_path") == CRawQWorkflow::UnableToReadCDFFile);
}
