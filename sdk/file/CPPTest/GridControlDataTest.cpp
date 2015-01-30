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


#include "file/CPPTest/GridControlDataTest.h"
//
#include "file/GridControlData.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_grid_control;

CPPUNIT_TEST_SUITE_REGISTRATION( GridControlDataTest );

void GridControlDataTest::setUp()
{
}

void GridControlDataTest::tearDown()
{
}

void GridControlDataTest::testCreation()
{
	GridControlData *grid = new GridControlData;
	CPPUNIT_ASSERT(1);
	delete grid;
}

void GridControlDataTest::testproperty_Rows()
{
	GridControlData grid;
	grid.SetRows(10);
	CPPUNIT_ASSERT( grid.GetRows() == 10);
}

void GridControlDataTest::testproperty_Columns()
{
	GridControlData grid;
	grid.SetColumns(10);
	CPPUNIT_ASSERT( grid.GetColumns() == 10);
}

void GridControlDataTest::testproperty_B1()
{
	GridControlData grid;
	grid.ResizeB1(3);
	CPPUNIT_ASSERT( grid.GetNumB1Probes() == 3);

	FeatureCoordinate coord;

	coord.x = 1;
	coord.y = 2;
	grid.SetB1(0, coord);
	coord.x = 3;
	coord.y = 4;
	grid.SetB1(1, coord);
	coord.x = 5;
	coord.y = 6;
	grid.SetB1(2, coord);

	CPPUNIT_ASSERT( grid.GetB1(0).x == 1);
	CPPUNIT_ASSERT( grid.GetB1(0).y == 2);
	CPPUNIT_ASSERT( grid.GetB1(1).x == 3);
	CPPUNIT_ASSERT( grid.GetB1(1).y == 4);
	CPPUNIT_ASSERT( grid.GetB1(2).x == 5);
	CPPUNIT_ASSERT( grid.GetB1(2).y == 6);
}

void GridControlDataTest::testproperty_B2()
{
	GridControlData grid;
	grid.ResizeB2(3);
	CPPUNIT_ASSERT( grid.GetNumB2Probes() == 3);

	FeatureCoordinate coord;

	coord.x = 1;
	coord.y = 2;
	grid.SetB2(0, coord);
	coord.x = 3;
	coord.y = 4;
	grid.SetB2(1, coord);
	coord.x = 5;
	coord.y = 6;
	grid.SetB2(2, coord);

	CPPUNIT_ASSERT( grid.GetB2(0).x == 1);
	CPPUNIT_ASSERT( grid.GetB2(0).y == 2);
	CPPUNIT_ASSERT( grid.GetB2(1).x == 3);
	CPPUNIT_ASSERT( grid.GetB2(1).y == 4);
	CPPUNIT_ASSERT( grid.GetB2(2).x == 5);
	CPPUNIT_ASSERT( grid.GetB2(2).y == 6);
}

void GridControlDataTest::testproperty_NS()
{
	GridControlData grid;
	grid.ResizeNS(3);
	CPPUNIT_ASSERT( grid.GetNumNSProbes() == 3);

	FeatureCoordinate coord;

	coord.x = 1;
	coord.y = 2;
	grid.SetNS(0, coord);
	coord.x = 3;
	coord.y = 4;
	grid.SetNS(1, coord);
	coord.x = 5;
	coord.y = 6;
	grid.SetNS(2, coord);

	CPPUNIT_ASSERT( grid.GetNS(0).x == 1);
	CPPUNIT_ASSERT( grid.GetNS(0).y == 2);
	CPPUNIT_ASSERT( grid.GetNS(1).x == 3);
	CPPUNIT_ASSERT( grid.GetNS(1).y == 4);
	CPPUNIT_ASSERT( grid.GetNS(2).x == 5);
	CPPUNIT_ASSERT( grid.GetNS(2).y == 6);
}

void GridControlDataTest::testmethod_Clear()
{
	GridControlData grid;
	grid.SetRows(10);
	grid.SetColumns(11);
	grid.ResizeNS(3);
	grid.ResizeB1(10);
	grid.ResizeB2(99);
	grid.Clear();
	CPPUNIT_ASSERT( grid.GetRows() == 0);
	CPPUNIT_ASSERT( grid.GetColumns() == 0);
	CPPUNIT_ASSERT( grid.GetNumB1Probes() == 0);
	CPPUNIT_ASSERT( grid.GetNumB2Probes() == 0);
	CPPUNIT_ASSERT( grid.GetNumNSProbes() == 0);
}
