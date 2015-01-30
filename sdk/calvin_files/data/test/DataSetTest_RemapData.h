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
 * This tests reading data from a large data set.  On Windows
 * this means that the view will be re-mapped.
 */
class DataSetTest_RemapData : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataSetTest_RemapData);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (RemapTest);
	CPPUNIT_TEST (RemapBackwardTest);
	CPPUNIT_TEST (RemapSizeTest);

	CPPUNIT_TEST_SUITE_END();

public:
	DataSetTest_RemapData();
	~DataSetTest_RemapData();

	void setUp();
	void tearDown();

	void CreationTest();
	void RemapTest();
	void RemapBackwardTest();
	void RemapSizeTest();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataSet* dataPlane;
};
