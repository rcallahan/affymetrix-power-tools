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

#ifndef __GENERICDATATEST_H_
#define __GENERICDATATEST_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * These are the general GenericData tests.
 */
class GenericDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericDataTest );

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (FileIdentifierTest);
	CPPUNIT_TEST (ArrayFileIdentifierTest);
	CPPUNIT_TEST (HeaderTest);
	CPPUNIT_TEST (DataGroupCntTest);
	CPPUNIT_TEST (DataGroupNamesTest);
	CPPUNIT_TEST (DataSetCntByIndexTest);
	CPPUNIT_TEST (DataSetCntByIndexErrorTest);
	CPPUNIT_TEST (DataSetCntByNameTest);
	CPPUNIT_TEST (DataSetCntByNameErrorTest);
	CPPUNIT_TEST (DataSetNamesByGroupIndexTest);
	CPPUNIT_TEST (DataSetNamesByGroupIndexErrorTest);
	CPPUNIT_TEST (DataSetNamesByGroupNameTest);
	CPPUNIT_TEST (DataSetNamesByGroupNameErrorTest);
	CPPUNIT_TEST (DataSetByIndexUsingMMTest);
	CPPUNIT_TEST (DataSetByNameUsingMMTest);
	CPPUNIT_TEST (DataSetByIndexUsingMMErrorTest);
	CPPUNIT_TEST (DataSetByNameUsingMMErrorTest);
	CPPUNIT_TEST (DataSetByIndexUsingFStreamTest);
	CPPUNIT_TEST (DataSetByNameUsingFStreamTest);
	CPPUNIT_TEST (DataSetByIndexUsingFStreamWithLoadEntireDataSetTest);
	CPPUNIT_TEST (DataSetByNameUsingFStreamWithLoadEntireDataSetTest);
	CPPUNIT_TEST (DataSetByIndexUsingFStreamErrorTest);
	CPPUNIT_TEST (DataSetByNameUsingFStreamErrorTest);
	CPPUNIT_TEST (DataSetByIndexUsingMMAndFStreamTest);
	CPPUNIT_TEST (DataSetByNameUsingMMAndFStreamTest);
	CPPUNIT_TEST (DataGroupUsingMMTest);
	CPPUNIT_TEST (DataGroupUsingFStreamTest);
	CPPUNIT_TEST (DataGroupUsingFStreamLoadEntireTest);
	CPPUNIT_TEST (DataGroupUsingMMAndFStreamTest);
	CPPUNIT_TEST (FindDataGroupHeaderByNameTest);
	CPPUNIT_TEST (FindDataGroupHeaderByIndexTest);
	CPPUNIT_TEST (FindDataSetHeaderByNameTest);
	CPPUNIT_TEST (FindDataSetHeaderByIndexTest);
	CPPUNIT_TEST (ClearTest);

	CPPUNIT_TEST_SUITE_END();

public:
	GenericDataTest();
	~GenericDataTest();

	void setUp();
	void tearDown();

	void CreationTest();
	void FileIdentifierTest();
	void ArrayFileIdentifierTest();
	void HeaderTest();
	void DataGroupCntTest();
	void DataGroupNamesTest();
	void DataSetCntByIndexTest();
	void DataSetCntByIndexErrorTest();
	void DataSetCntByNameTest();
	void DataSetCntByNameErrorTest();
	void DataSetNamesByGroupIndexTest();
	void DataSetNamesByGroupIndexErrorTest();
	void DataSetNamesByGroupNameTest();
	void DataSetNamesByGroupNameErrorTest();
	void DataSetByIndexUsingMMTest();
	void DataSetByNameUsingMMTest();
	void DataSetByIndexUsingMMErrorTest();
	void DataSetByNameUsingMMErrorTest();
	void DataSetByIndexUsingFStreamTest();
	void DataSetByNameUsingFStreamTest();
	void DataSetByIndexUsingFStreamWithLoadEntireDataSetTest();
	void DataSetByNameUsingFStreamWithLoadEntireDataSetTest();
	void DataSetByIndexUsingFStreamErrorTest();
	void DataSetByNameUsingFStreamErrorTest();
	void DataSetByIndexUsingMMAndFStreamTest();
	void DataSetByNameUsingMMAndFStreamTest();
	void DataGroupUsingMMTest();
	void DataGroupUsingFStreamTest();
	void DataGroupUsingFStreamLoadEntireTest();
	void DataGroupUsingMMAndFStreamTest();
	void FindDataGroupHeaderByNameTest();
	void FindDataGroupHeaderByIndexTest();
	void FindDataSetHeaderByNameTest();
	void FindDataSetHeaderByIndexTest();
	void ClearTest();

private:
	affymetrix_calvin_io::GenericData* data;
};


#endif // __GENERICDATATEST_H_
