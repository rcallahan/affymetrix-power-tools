////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


#include "calvin_files/array/test/PATAssignmentTest.h"
//
#include "calvin_files/array/src/PATAssignment.h"
//

using namespace affymetrix_calvin_array;
using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( PATAssignmentMethodTest );

void PATAssignmentMethodTest::setUp()
{
}

void PATAssignmentMethodTest::tearDown()
{
}

void PATAssignmentMethodTest::test_MediaToString()
{
	CPPUNIT_ASSERT(PATAssignmentMethodToString(NoAssignment) == L"None");
	CPPUNIT_ASSERT(PATAssignmentMethodToString(AffyBarcodeAssignment) == L"AffyBarcode");
	CPPUNIT_ASSERT(PATAssignmentMethodToString(UserSelectedAssignment) == L"UserSelected");
	CPPUNIT_ASSERT(PATAssignmentMethodToString(OtherAssignment) == L"Other");
}

void PATAssignmentMethodTest::test_MediaFromString()
{
	CPPUNIT_ASSERT(PATAssignmentMethodFromString(L"None") == NoAssignment);
	CPPUNIT_ASSERT(PATAssignmentMethodFromString(L"AffyBarcode") == AffyBarcodeAssignment);
	CPPUNIT_ASSERT(PATAssignmentMethodFromString(L"UserSelected") == UserSelectedAssignment);
	CPPUNIT_ASSERT(PATAssignmentMethodFromString(L"Other") == OtherAssignment);
}

