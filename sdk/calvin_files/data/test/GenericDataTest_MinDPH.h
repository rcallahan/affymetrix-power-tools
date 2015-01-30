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

#ifndef __GENERICDATATEST_MINDPH_H_
#define __GENERICDATATEST_MINDPH_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests the GenericData class when the minimal DataSet
 * information has been read.
 */
class GenericDataTest_MinDPH : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericDataTest_MinDPH );

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (DataSetByIndexTest);
	CPPUNIT_TEST (DataSetByNameTest);
	CPPUNIT_TEST (FileIdentifierTest);
	CPPUNIT_TEST (ArrayFileIdentifierTest);
	CPPUNIT_TEST (HeaderTest);
	CPPUNIT_TEST (DataGroupNamesTest);
	CPPUNIT_TEST (DataSetCntTest);
	CPPUNIT_TEST (DataSetNamesTest);

	CPPUNIT_TEST_SUITE_END();

public:
	GenericDataTest_MinDPH();
	~GenericDataTest_MinDPH();

	void setUp();
	void tearDown();

	void CreationTest();
	void DataSetByIndexTest();
	void DataSetByNameTest();
	void FileIdentifierTest();
	void ArrayFileIdentifierTest();
	void HeaderTest();
	void DataGroupNamesTest();
	void DataSetCntTest();
	void DataSetNamesTest();

private:
	affymetrix_calvin_io::GenericData* data;

};

// tests to write:
// FileIdentifier
// ArrayFileIdentifier
// Header
// DataGroupNames
// DataSetCnt
// DataSetNames
// 
// Multiset access when only one set has been read.

#endif // __GENERICDATATEST_MINDPH_H_
