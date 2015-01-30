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

#ifndef __GENERICFILEWRITERTEST_H_
#define __GENERICFILEWRITERTEST_H_

#include "calvin_files/writers/src/GenericFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class GenericFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericFileWriterTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( WriteTest );

	CPPUNIT_TEST_SUITE_END();

private:
    
    GenericFileWriter *writer;
    FileHeader* hdr;

public:

    void setUp();
	void tearDown();
	void testCreation();
	void WriteTest();
};

#endif // __GENERICFILEWRITERTEST_H_
