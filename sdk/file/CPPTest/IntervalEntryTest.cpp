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

#include "file/CPPTest/IntervalEntryTest.h"
//
#include "file/IntervalEntry.h"
//
#include <cmath>
//

CPPUNIT_TEST_SUITE_REGISTRATION( IntervalEntryTest );

void IntervalEntryTest::setUp()
{
}

void IntervalEntryTest::tearDown()
{
}

void IntervalEntryTest::testSize()
{
	IntervalEntry i;
	i.start = 1;
	i.stop = 10;
	CPPUNIT_ASSERT(i.size() == 9 );
}

void IntervalEntryTest::testAssignmentOperator()
{
	IntervalEntry r1;
	r1.start = 1;
	r1.stop = 2;
	r1.seq = "abc";
	r1.probeSetName = "name";
	r1.overlap = 12.1f;
	r1.strand = '+';
	
	IntervalEntry res;
	res = r1;
	CPPUNIT_ASSERT(res.overlap == 12.1f);
	CPPUNIT_ASSERT(res.start == 1);
	CPPUNIT_ASSERT(res.stop == 2);
	CPPUNIT_ASSERT(res.seq == "abc");
	CPPUNIT_ASSERT(res.probeSetName == "name");
	CPPUNIT_ASSERT(res.strand == '+');
}
