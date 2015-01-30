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

//
#include "calvin_files/data/test/ProbeSetMultiDataDataTest.h"
//
#include "calvin_files/data/src/ProbeSetMultiDataData.h"
//

using namespace std;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( ProbeSetMultiDataDataTest );

void ProbeSetMultiDataDataTest::setUp()
{
}

void ProbeSetMultiDataDataTest::tearDown()
{
}

void ProbeSetMultiDataDataTest::test_ChromosomeFromString()
{
	CPPUNIT_ASSERT( ChromosomeFromString("X") == X_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("Y") == Y_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("MT") == MT_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("--") == NO_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("-") == NO_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("") == NO_CHR);
	CPPUNIT_ASSERT( ChromosomeFromString("1") == 1);
	CPPUNIT_ASSERT( ChromosomeFromString("2") == 2);
	CPPUNIT_ASSERT( ChromosomeFromString("3") == 3);
	CPPUNIT_ASSERT( ChromosomeFromString("4") == 4);
	CPPUNIT_ASSERT( ChromosomeFromString("5") == 5);
	CPPUNIT_ASSERT( ChromosomeFromString("6") == 6);
	CPPUNIT_ASSERT( ChromosomeFromString("7") == 7);
	CPPUNIT_ASSERT( ChromosomeFromString("8") == 8);
	CPPUNIT_ASSERT( ChromosomeFromString("9") == 9);
	CPPUNIT_ASSERT( ChromosomeFromString("10") == 10);
}

void ProbeSetMultiDataDataTest::test_ChromosomeToString()
{
	CPPUNIT_ASSERT(ChromosomeToString(X_CHR) == "X");
	CPPUNIT_ASSERT(ChromosomeToString(Y_CHR) == "Y");
	CPPUNIT_ASSERT(ChromosomeToString(MT_CHR) == "MT");
	CPPUNIT_ASSERT(ChromosomeToString(NO_CHR) == "-");
	CPPUNIT_ASSERT(ChromosomeToString(1) == "1");
	CPPUNIT_ASSERT(ChromosomeToString(2) == "2");
	CPPUNIT_ASSERT(ChromosomeToString(3) == "3");
	CPPUNIT_ASSERT(ChromosomeToString(4) == "4");
	CPPUNIT_ASSERT(ChromosomeToString(5) == "5");
	CPPUNIT_ASSERT(ChromosomeToString(6) == "6");
	CPPUNIT_ASSERT(ChromosomeToString(7) == "7");
	CPPUNIT_ASSERT(ChromosomeToString(8) == "8");
	CPPUNIT_ASSERT(ChromosomeToString(9) == "9");
	CPPUNIT_ASSERT(ChromosomeToString(10) == "10");
}

void ProbeSetMultiDataDataTest::test_CytoCallToString()
{
	CPPUNIT_ASSERT(CytoCallToString(CYTO_ABSENT_CALL) == "A");
	CPPUNIT_ASSERT(CytoCallToString(CYTO_PRESENT_CALL) == "P");
	CPPUNIT_ASSERT(CytoCallToString(CYTO_NO_CALL) == "NC");
}

void ProbeSetMultiDataDataTest::test_CytoCallFromString()
{
	CPPUNIT_ASSERT(CytoCallFromString("A") == CYTO_ABSENT_CALL);
	CPPUNIT_ASSERT(CytoCallFromString("P") == CYTO_PRESENT_CALL);
    CPPUNIT_ASSERT(CytoCallFromString("NC") == CYTO_NO_CALL);
}
