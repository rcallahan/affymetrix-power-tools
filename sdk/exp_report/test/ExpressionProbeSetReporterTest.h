////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#ifndef _ExpressionProbeSetReporterTEST_HEADER_
#define _ExpressionProbeSetReporterTEST_HEADER_

#include <cppunit/extensions/HelperMacros.h>

class ExpressionProbeSetReporterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ExpressionProbeSetReporterTest );

	CPPUNIT_TEST( testDetection );
	CPPUNIT_TEST( testChange );
	CPPUNIT_TEST( testControls );
	CPPUNIT_TEST( testControlProbeSet );
	CPPUNIT_TEST( testProbeSetSignals );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testDetection();
	void testChange();
	void testControls();
    void testControlProbeSet();
    void testProbeSetSignals();
};


#endif
