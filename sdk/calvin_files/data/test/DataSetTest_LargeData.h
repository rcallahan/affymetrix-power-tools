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
#pragma once

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests reading data from a large data set.
 */
class DataSetTest_LargeData : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataSetTest_LargeData);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (GetDataRaw300MBTest);
	CPPUNIT_TEST (GetData300MBTest);

	CPPUNIT_TEST_SUITE_END();

public:
	DataSetTest_LargeData();
	~DataSetTest_LargeData();

	void setUp();
	void tearDown();

	void CreationTest();
	void GetDataRaw300MBTest();
	void GetData300MBTest();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataSet* dataPlane;
};
