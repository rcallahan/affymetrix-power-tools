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

#include "file/CPPTest/TagValuePairTest.h"
//
#include "file/TagValuePair.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( TagValuePairTest );


void TagValuePairTest::setUp()
{
}

void TagValuePairTest::tearDown()
{
}

void TagValuePairTest::testmethod_TagValuePairType_assignment_operator()
{
	TagValuePairType orig;
	TagValuePairType copy;
	orig.Tag = "tag";
	orig.Value = "value";
	copy = orig;
	CPPUNIT_ASSERT(copy.Tag == orig.Tag);
	CPPUNIT_ASSERT(copy.Value == orig.Value);
}

void TagValuePairTest::testmethod_TagValuePairType_equality_operators()
{
	TagValuePairType orig;
	TagValuePairType copy;
	orig.Tag = "tag";
	orig.Value = "value";
	copy = orig;
	CPPUNIT_ASSERT( copy == orig );
	CPPUNIT_ASSERT( copy == "tag" );
}
