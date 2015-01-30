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

#ifndef __DATAGROUPTEST_H_
#define __DATAGROUPTEST_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests the basic functions of the DataGroup class.
 */
class DataGroupTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataGroupTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (HeaderTest);
	CPPUNIT_TEST (DataSetByIndexTest);
	CPPUNIT_TEST (DataSetByNameTest);

	CPPUNIT_TEST_SUITE_END();

public:
	DataGroupTest();
	~DataGroupTest();

	void setUp();
	void tearDown();

	void CreationTest();
	void HeaderTest();
	void DataSetByIndexTest();
	void DataSetByNameTest();

// Support methods
	void* MapFile(const std::string& fileName);
	void UnmapFile();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataGroup* dc;

#ifdef _MSC_VER
	HANDLE fileHandle;
	HANDLE fileMapHandle;
#endif

};

#endif // __DATAGROUPTEST_H_
