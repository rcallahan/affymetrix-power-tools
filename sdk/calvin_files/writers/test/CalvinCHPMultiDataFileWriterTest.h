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

#ifndef __CALVINCHPMultiDataFILEWRITERTEST_H_
#define __CALVINCHPMultiDataFILEWRITERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CHPMultiDataFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CHPMultiDataFileWriterTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( WriteTestGeno );
	CPPUNIT_TEST ( WriteTestCN );
	CPPUNIT_TEST ( WriteTestExp );
	CPPUNIT_TEST ( WriteTestAll );
	CPPUNIT_TEST ( WriteTestCNV );
	CPPUNIT_TEST ( WriteTestDMET );
    CPPUNIT_TEST ( WriteTestChrSummary );
    CPPUNIT_TEST ( WriteTestChrSegment );
    CPPUNIT_TEST ( WriteTestChrSegmentEx );
    CPPUNIT_TEST ( WriteTestFamilial );
	CPPUNIT_TEST ( WriteTestAllelePeaks );
	CPPUNIT_TEST ( WriteTestMarkerABSignals );
	CPPUNIT_TEST ( WriteTestCytoGeno );
	CPPUNIT_TEST_SUITE_END();

public:

	void setUp();
	void tearDown();
	void testCreation();
	void WriteTestGeno();
	void WriteTestExp();
	void WriteTestAll();
	void WriteTestCN();
	void WriteTestCNV();
	void WriteTestDMET();
    void WriteTestChrSummary();
    void WriteTestChrSegment();
    void WriteTestChrSegmentEx();
    void WriteTestFamilial();
	void WriteTestAllelePeaks();
	void WriteTestMarkerABSignals();
	void WriteTestCytoGeno();
};

#endif // __CALVINCHPMultiDataFILEWRITERTEST_H_
