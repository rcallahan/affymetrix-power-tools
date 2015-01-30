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


#include "calvin_files/array/test/CreateStepTest.h"
//
#include "calvin_files/array/src/CreateStep.h"
//

using namespace affymetrix_calvin_array;
using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( CreateStepTest );

void CreateStepTest::setUp()
{
}

void CreateStepTest::tearDown()
{
}

void CreateStepTest::test_MediaToString()
{
	CPPUNIT_ASSERT(CreateStepToString(NoStep) == L"None");
	CPPUNIT_ASSERT(CreateStepToString(ArrayRegistrationStep) == L"ArrayRegistration");
	CPPUNIT_ASSERT(CreateStepToString(ScanningStep) == L"Scanning");
	CPPUNIT_ASSERT(CreateStepToString(GriddingStep) == L"Gridding");
	CPPUNIT_ASSERT(CreateStepToString(CELAnalysisStep) == L"CELAnalysis");
	CPPUNIT_ASSERT(CreateStepToString(OtherStep) == L"Other");
	CPPUNIT_ASSERT(CreateStepToString(FromStep) == L"From");
	CPPUNIT_ASSERT(CreateStepToString(JobOrderServerStep) == L"JobOrderServer");
	CPPUNIT_ASSERT(CreateStepToString(FileIndexerStep) == L"FileIndexer");


}

void CreateStepTest::test_MediaFromString()
{
	CPPUNIT_ASSERT(CreateStepFromString(L"None") == NoStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"ArrayRegistration") ==ArrayRegistrationStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"Scanning") == ScanningStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"Gridding") == GriddingStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"CELAnalysis") == CELAnalysisStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"Other") == OtherStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"From") == FromStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"JobOrderServer") == JobOrderServerStep);
	CPPUNIT_ASSERT(CreateStepFromString(L"FileIndexer") == FileIndexerStep);
}

