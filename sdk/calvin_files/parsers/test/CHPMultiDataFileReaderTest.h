////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#ifndef __CHPGenericGenotypeFILEREADERTEST_H_
#define __CHPGenericGenotypeFILEREADERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CHPMultiDataFileReaderTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( CHPMultiDataFileReaderTest );

	CPPUNIT_TEST ( testCreation );
    CPPUNIT_TEST ( testReadCN );
	CPPUNIT_TEST ( testRead );
    CPPUNIT_TEST ( testReadCNV );
    CPPUNIT_TEST ( testReadCyto );


	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void testCreation();
	void testRead();
    void testReadCN();
    void testReadCNV();
	void testReadCyto();
};

#endif // __CHPGenericGenotypeFILEREADERTEST_H_
