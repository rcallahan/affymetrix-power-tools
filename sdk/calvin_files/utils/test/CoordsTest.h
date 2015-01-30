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

#ifndef __COORDSTEST_H_
#define __COORDSTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CoordsTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (CoordsTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (CastFRegionToFGridCoordsTest);
	CPPUNIT_TEST (CastRegionToGridCoordsTest);
	CPPUNIT_TEST (FGridCoordsIsEmptyTest);
	CPPUNIT_TEST (GridCoordsIsEmptyTest);
	CPPUNIT_TEST (CastFGridCoordsToFRegionTest);
	CPPUNIT_TEST (CastGridCoordsToRegionTest);
	CPPUNIT_TEST (CastFRegionTooSmallToFGridCoordsTest);
	CPPUNIT_TEST (CastRegionTooSmallToGridCoordsTest);
	CPPUNIT_TEST (FRegionClearTest);
	CPPUNIT_TEST (RegionClearTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void CastFRegionToFGridCoordsTest();
	void CastRegionToGridCoordsTest();
	void FGridCoordsIsEmptyTest();
	void GridCoordsIsEmptyTest();
	void CastFGridCoordsToFRegionTest();
	void CastGridCoordsToRegionTest();
	void CastFRegionTooSmallToFGridCoordsTest();
	void CastRegionTooSmallToGridCoordsTest();
	void FRegionClearTest();
	void RegionClearTest();
};

#endif // __COORDSTEST_H_
