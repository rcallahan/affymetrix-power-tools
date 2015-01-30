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

#ifndef __FUSIONGRCFILETEST_H_
#define __FUSIONGRCFILETEST_H_

#include "calvin_files/fusion/src/FusionGrcFileReader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_fusion_io;

class FusionGrcFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FusionGrcFileReaderTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testmethod_Read_invalid_file );
	CPPUNIT_TEST ( testmethod_Read_missing_file );
	CPPUNIT_TEST ( testmethod_Read_calvin_grc_file );
	CPPUNIT_TEST ( testmethod_Read_gcos_grc_file );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testmethod_Read_invalid_file();
	void testmethod_Read_missing_file();
	void testmethod_Read_calvin_grc_file();
	void testmethod_Read_gcos_grc_file();
};

#endif // __FUSIONGRCFILETEST_H_
