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

#ifndef __CALVINCHPFILEWRITERTEST_H_
#define __CALVINCHPFILEWRITERTEST_H_

#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/writers/src/CalvinCHPFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class CHPFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CHPFileWriterTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( WriteExpressionTest );
	CPPUNIT_TEST ( WriteGenotypeTest );
	CPPUNIT_TEST ( WriteUniversalTest );
	CPPUNIT_TEST ( WriteReseqTest );

	CPPUNIT_TEST_SUITE_END();

private:

	CHPFileWriter* writer;

public:

	void setUp();
	void tearDown();
	void testCreation();
	void WriteExpressionTest();
	void WriteGenotypeTest();
	void WriteUniversalTest();
	void WriteReseqTest();
    void WriteCopyNumberVariationTest();
};

#endif // __CALVINCHPFILEWRITERTEST_H_
