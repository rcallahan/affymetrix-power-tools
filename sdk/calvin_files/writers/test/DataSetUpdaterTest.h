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

#ifndef __DATASETUPDATERTEST_H_
#define __DATASETUPDATERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class DataSetUpdaterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DataSetUpdaterTest );

	CPPUNIT_TEST ( testUpdateFail );
	CPPUNIT_TEST ( testUpdateInt8 );
	CPPUNIT_TEST ( testUpdateInt16 );
	CPPUNIT_TEST ( testUpdateInt32 );
	CPPUNIT_TEST ( testUpdateUInt8 );
	CPPUNIT_TEST ( testUpdateUInt16 );
	CPPUNIT_TEST ( testUpdateUInt32 );
	CPPUNIT_TEST ( testUpdateFloat );
	CPPUNIT_TEST ( testUpdateString8 );
	CPPUNIT_TEST ( testUpdateString16 );
	CPPUNIT_TEST_SUITE_END();

	void CreateReferenceFile();

public:
	void setUp();
	void tearDown();

	void testUpdateFail();
	void testUpdateInt8();
	void testUpdateInt16();
	void testUpdateInt32();
	void testUpdateUInt8();
	void testUpdateUInt16();
	void testUpdateUInt32();
	void testUpdateFloat();
	void testUpdateString8();
	void testUpdateString16();
};

#endif // __DATASETUPDATERTEST_H_
