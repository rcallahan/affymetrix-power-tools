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
#pragma once

#include <cppunit/extensions/HelperMacros.h>
//
#include <cstring>
#include <string>
//

using namespace std;

// This test fixture tests the FusionCHPData class
// with GCOS CHP data.

class FusionCHPDataCalvinTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (FusionCHPDataCalvinTest);

	CPPUNIT_TEST (TestAll);
	CPPUNIT_TEST (TestFileId);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void TestAll();

	void CreationTest(string fileName);
	void ReadHeaderTest(string fileName);
	void RowsColsTest(string fileName);
	void VersionTest(string fileName);
	void AlgorithmNameAndVersionTest(string fileName);
	void ChipTypeTest(string fileName);
    void ParentCelFileName(string fileName);
    void ProgIDTest(string fileName);
	void AlgParamsTest(string fileName);
	void SummaryParamsTest(string fileName);
	void NumberOfSetsTest(string fileName);
	void ReadHeaderAndVerifyAlgParamsTest(string fileName);
	void ReadHeaderAndVerifySummaryParamsTest(string fileName);
	void ReadTest(string fileName);
	void CallMethodsWhenObjectIsNotReadyTest(string fileName);
	void DataCompareTest(string fileName);
	void BackgroundDataTest(string fileName);
	void TestFileId();
};
