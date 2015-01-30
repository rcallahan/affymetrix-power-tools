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

#ifndef __CDFDATATEST_H_
#define __CDFDATATEST_H_

//
#include "calvin_files/data/src/CDFData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class CDFDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CDFDataTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( FilenameTest );
	CPPUNIT_TEST ( DataTypeIdTest );
	CPPUNIT_TEST (ProbeSetCntTest);
	CPPUNIT_TEST (ArrayRowsTest);
	CPPUNIT_TEST (ArrayColsTest);
	CPPUNIT_TEST (RefSequenceTest);
	CPPUNIT_TEST (FileHeaderTest);
	CPPUNIT_TEST (GetGenericDataTest);

	CPPUNIT_TEST_SUITE_END();

private:

	CDFData *data;

public:
	void setUp();
	void tearDown();

	void testCreation();
	void FilenameTest();
	void DataTypeIdTest();
	void ProbeSetCntTest();
	void ArrayRowsTest();
	void ArrayColsTest();
	void RefSequenceTest();
	void FileHeaderTest();
	void GetGenericDataTest();
};

#endif // __CDFDATATEST_H_
