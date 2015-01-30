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

#include "calvin_files/utils/test/CoordsTest.h"
//
#include "calvin_files/utils/src/Coords.h"
//

using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( CoordsTest );

void CoordsTest::setUp()
{
}

void CoordsTest::tearDown()
{
}

void CoordsTest::CreationTest()
{
	FPoint fpt;
    // quell compiler warning of unused var
    fpt.x = 8.0f;
	FRegion frgn;
	FGridCoords fgrid;
	Region rgn;
	Point pt;
    // quell compiler warning of unused var
    pt.x = 8;
	CPPUNIT_ASSERT( 1 );
}

void CoordsTest::CastFRegionToFGridCoordsTest()
{
	FRegion rgn;
	FPoint pt = {8.0f, 7.0f};
	rgn.pts.push_back(pt);
	pt.x = 3.0f;
	rgn.pts.push_back(pt);
	pt.y = 2.0f;
	rgn.pts.push_back(pt);
	pt.x = 11.0f;
	rgn.pts.push_back(pt);

	FGridCoords grid;
	grid = (FGridCoords)rgn;

	CPPUNIT_ASSERT(grid.upperleft.x == 8.0f);
	CPPUNIT_ASSERT(grid.upperleft.y == 7.0f);
	CPPUNIT_ASSERT(grid.upperright.x == 3.0f);
	CPPUNIT_ASSERT(grid.upperright.y == 7.0f);
	CPPUNIT_ASSERT(grid.lowerright.x == 3.0f);
	CPPUNIT_ASSERT(grid.lowerright.y == 2.0f);
	CPPUNIT_ASSERT(grid.lowerleft.x == 11.0f);
	CPPUNIT_ASSERT(grid.lowerleft.y == 2.0f);
}

void CoordsTest::CastRegionToGridCoordsTest()
{
	Region rgn;
	Point pt = {8, 7};
	rgn.pts.push_back(pt);
	pt.x = 3;
	rgn.pts.push_back(pt);
	pt.y = 2;
	rgn.pts.push_back(pt);
	pt.x = 11;
	rgn.pts.push_back(pt);

	GridCoords grid;
	grid = (GridCoords)rgn;

	CPPUNIT_ASSERT(grid.upperleft.x == 8);
	CPPUNIT_ASSERT(grid.upperleft.y == 7);
	CPPUNIT_ASSERT(grid.upperright.x == 3);
	CPPUNIT_ASSERT(grid.upperright.y == 7);
	CPPUNIT_ASSERT(grid.lowerright.x == 3);
	CPPUNIT_ASSERT(grid.lowerright.y == 2);
	CPPUNIT_ASSERT(grid.lowerleft.x == 11);
	CPPUNIT_ASSERT(grid.lowerleft.y == 2);
}

void CoordsTest::FGridCoordsIsEmptyTest()
{
	FGridCoords grid;
	CPPUNIT_ASSERT(grid.IsEmpty());
	grid.lowerleft.x = 1.0f;
	grid.lowerleft.y = 1.0f;
	grid.lowerright = grid.upperright = grid.upperleft = grid.lowerleft;
	CPPUNIT_ASSERT(grid.IsEmpty());
}
void CoordsTest::GridCoordsIsEmptyTest()
{
	GridCoords grid;
	CPPUNIT_ASSERT(grid.IsEmpty());
	grid.lowerleft.x = 1;
	grid.lowerleft.y = 1;
	grid.lowerright = grid.upperright = grid.upperleft = grid.lowerleft;
	CPPUNIT_ASSERT(grid.IsEmpty());
}

void CoordsTest::CastFGridCoordsToFRegionTest()
{
	FGridCoords grid;
	grid.lowerleft.x = 1.0f;
	grid.lowerleft.y = 1.2f;
	grid.lowerright.x = 11.0f;
	grid.lowerright.y = 14.2f;
	grid.upperright.x = 111.0f;
	grid.upperright.y = 114.2f;
	grid.upperleft.x = 1001.0f;
	grid.upperleft.y = 1200.2f;

	FRegion r = (FRegion)grid;
	CPPUNIT_ASSERT(r.pts.size() == 4);
	CPPUNIT_ASSERT(r.pts[LowerLeft].x == 1.0f);
	CPPUNIT_ASSERT(r.pts[LowerLeft].y == 1.2f);
	CPPUNIT_ASSERT(r.pts[LowerRight].x == 11.0f);
	CPPUNIT_ASSERT(r.pts[LowerRight].y == 14.2f);
	CPPUNIT_ASSERT(r.pts[UpperRight].x == 111.0f);
	CPPUNIT_ASSERT(r.pts[UpperRight].y == 114.2f);
	CPPUNIT_ASSERT(r.pts[UpperLeft].x == 1001.0f);
	CPPUNIT_ASSERT(r.pts[UpperLeft].y == 1200.2f);
}

void CoordsTest::CastGridCoordsToRegionTest()
{
	GridCoords grid;
	grid.lowerleft.x = 1;
	grid.lowerleft.y = 1;
	grid.lowerright.x = 11;
	grid.lowerright.y = 14;
	grid.upperright.x = 111;
	grid.upperright.y = 114;
	grid.upperleft.x = 1001;
	grid.upperleft.y = 1200;

	Region r = (Region)grid;
	CPPUNIT_ASSERT(r.pts.size() == 4);
	CPPUNIT_ASSERT(r.pts[LowerLeft].x == 1);
	CPPUNIT_ASSERT(r.pts[LowerLeft].y == 1);
	CPPUNIT_ASSERT(r.pts[LowerRight].x == 11);
	CPPUNIT_ASSERT(r.pts[LowerRight].y == 14);
	CPPUNIT_ASSERT(r.pts[UpperRight].x == 111);
	CPPUNIT_ASSERT(r.pts[UpperRight].y == 114);
	CPPUNIT_ASSERT(r.pts[UpperLeft].x == 1001);
	CPPUNIT_ASSERT(r.pts[UpperLeft].y == 1200);
}

void CoordsTest::CastFRegionTooSmallToFGridCoordsTest()
{
	FRegion rgn;
	FPoint pt = {8.0f, 7.0f};
	rgn.pts.push_back(pt);
	pt.x = 3.0f;
	rgn.pts.push_back(pt);
	pt.y = 2.0f;
	rgn.pts.push_back(pt);

	FGridCoords grid = (FGridCoords)rgn;
	CPPUNIT_ASSERT(grid.IsEmpty());
}

void CoordsTest::CastRegionTooSmallToGridCoordsTest()
{
	Region rgn;
	Point pt = {8, 7};
	rgn.pts.push_back(pt);
	pt.x = 3;
	rgn.pts.push_back(pt);
	pt.y = 2;
	rgn.pts.push_back(pt);

	GridCoords grid = (GridCoords)rgn;
	CPPUNIT_ASSERT(grid.IsEmpty());
}

void CoordsTest::FRegionClearTest()
{
	FRegion rgn;
	FPoint pt = {8.0f, 7.0f};
	rgn.pts.push_back(pt);
	CPPUNIT_ASSERT(rgn.pts.size() == 1);
	rgn.Clear();
	CPPUNIT_ASSERT(rgn.pts.size() == 0);
}

void CoordsTest::RegionClearTest()
{
	Region rgn;
	Point pt = {8, 7};
	rgn.pts.push_back(pt);
	CPPUNIT_ASSERT(rgn.pts.size() == 1);
	rgn.Clear();
	CPPUNIT_ASSERT(rgn.pts.size() == 0);
}
