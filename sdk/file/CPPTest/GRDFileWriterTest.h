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

#ifndef __GRDFILEWRITERTEST_H_
#define __GRDFILEWRITERTEST_H_

#include "file/GRDFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affxgrd;

class CGRDFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CGRDFileWriterTest );
	
	CPPUNIT_TEST(WriteGRDFileTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void WriteGRDFileTest();

private:
};
#endif // __GRDFILEWRITERTEST_H_
