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

#ifndef __DATASETTEST_ALLCOLUMNTYPES_H_
#define __DATASETTEST_ALLCOLUMNTYPES_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests every that every column type can be accessed.
 */
class DataSetTest_AllColumnTypes : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataSetTest_AllColumnTypes);

	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testOpenDataSet);
	CPPUNIT_TEST (testAccessColumnValues);

	CPPUNIT_TEST_SUITE_END();

public:
	DataSetTest_AllColumnTypes();
	~DataSetTest_AllColumnTypes();

	void setUp();
	void tearDown();

	void testCreation();
	void testOpenDataSet();
	void testAccessColumnValues();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataSet* dataPlane;
};

#endif // __DATASETTEST_ALLCOLUMNTYPES_H_
