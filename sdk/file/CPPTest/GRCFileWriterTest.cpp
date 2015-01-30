////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include "file/CPPTest/GRCFileWriterTest.h"
//
#include "file/1LQFileData.h"
#include "file/GRCFileWriter.h"
//

using namespace affxgrc;
using namespace affx1lq;
using namespace std;
using namespace affymetrix_grid_control;

CPPUNIT_TEST_SUITE_REGISTRATION( CGRCFileWriterTest );

#define INPUT_FILE "./data/grcinput.1lq"
#define OUTPUT_FILE "./data/grcoutput.grc"
#define OPPOUTPUT_FILE "./data/grcoppoutput.grc"

void CGRCFileWriterTest::setUp()
{
}

void CGRCFileWriterTest::tearDown()
{
}

void CGRCFileWriterTest::testCreation()
{
	CGRCFileWriter grc;
	CPPUNIT_ASSERT( 1 );
}

void CGRCFileWriterTest::testmethod_Write()
{
	CGRCFileWriter grc;
	grc.SetFileName(OUTPUT_FILE);
	CPPUNIT_ASSERT(grc.GetFileName() == OUTPUT_FILE);

	C1LQFileData lq;
	lq.SetFileName(INPUT_FILE);
	lq.Read();
	CPPUNIT_ASSERT(grc.Write(lq) == true);

	GridControlData grid;
	grc.SetFileName(OUTPUT_FILE);
	CPPUNIT_ASSERT(grc.GetFileName() == OUTPUT_FILE);
	CPPUNIT_ASSERT(grc.Read(grid) == true);

	CPPUNIT_ASSERT( grid.GetNumNSProbes() == 1) ;
	CPPUNIT_ASSERT( grid.GetNumB1Probes() == 2) ;
	CPPUNIT_ASSERT( grid.GetNumB2Probes() == 3) ;

	CPPUNIT_ASSERT( grid.GetNS(0).x == 1);
	CPPUNIT_ASSERT( grid.GetNS(0).y == 2);

	CPPUNIT_ASSERT( grid.GetB1(0).x == 2);
	CPPUNIT_ASSERT( grid.GetB1(0).y == 3);
	CPPUNIT_ASSERT( grid.GetB1(1).x == 3);
	CPPUNIT_ASSERT( grid.GetB1(1).y == 4);

	CPPUNIT_ASSERT( grid.GetB2(0).x == 4);
	CPPUNIT_ASSERT( grid.GetB2(0).y == 5);
	CPPUNIT_ASSERT( grid.GetB2(1).x == 5);
	CPPUNIT_ASSERT( grid.GetB2(1).y == 6);
	CPPUNIT_ASSERT( grid.GetB2(2).x == 6);
	CPPUNIT_ASSERT( grid.GetB2(2).y == 7);

	// write opposite.  Resulting in B2 and B1 swapped.
	grc.SetFileName(OPPOUTPUT_FILE);
	CPPUNIT_ASSERT(grc.GetFileName() == OPPOUTPUT_FILE);
	CPPUNIT_ASSERT(grc.Write(true) == true);

	GridControlData gridopp;
	CPPUNIT_ASSERT(grc.Read(gridopp) == true);

	CPPUNIT_ASSERT( grid.GetNumNSProbes() == gridopp.GetNumNSProbes()) ;
	CPPUNIT_ASSERT( grid.GetNumB1Probes() == gridopp.GetNumB2Probes()) ;
	CPPUNIT_ASSERT( grid.GetNumB2Probes() == gridopp.GetNumB1Probes()) ;

	CPPUNIT_ASSERT( grid.GetNS(0).x == gridopp.GetNS(0).x);
	CPPUNIT_ASSERT( grid.GetNS(0).y == gridopp.GetNS(0).y);

	CPPUNIT_ASSERT( grid.GetB1(0).x == gridopp.GetB2(0).x);
	CPPUNIT_ASSERT( grid.GetB1(0).y == gridopp.GetB2(0).y);
	CPPUNIT_ASSERT( grid.GetB1(1).x == gridopp.GetB2(1).x);
	CPPUNIT_ASSERT( grid.GetB1(1).y == gridopp.GetB2(1).y);

	CPPUNIT_ASSERT( grid.GetB2(0).x == gridopp.GetB1(0).x);
	CPPUNIT_ASSERT( grid.GetB2(0).y == gridopp.GetB1(0).y);
	CPPUNIT_ASSERT( grid.GetB2(1).x == gridopp.GetB1(1).x);
	CPPUNIT_ASSERT( grid.GetB2(1).y == gridopp.GetB1(1).y);
	CPPUNIT_ASSERT( grid.GetB2(2).x == gridopp.GetB1(2).x);
	CPPUNIT_ASSERT( grid.GetB2(2).y == gridopp.GetB1(2).y);
	
}
