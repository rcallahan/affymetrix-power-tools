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

#ifndef __GENERICDATAHEADERTEST_H_
#define __GENERICDATAHEADERTEST_H_

#include "calvin_files/data/src/GenericDataHeader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class GenericDataHeaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericDataHeaderTest );

	CPPUNIT_TEST ( testCreation );
    CPPUNIT_TEST ( FileTypeIdTest );
    CPPUNIT_TEST ( FileIdTest );
    CPPUNIT_TEST ( FileCreationTimeTest );
    CPPUNIT_TEST ( LocaleTest );
    CPPUNIT_TEST ( NameValTest );
    CPPUNIT_TEST ( AddParentEntryTest );
    CPPUNIT_TEST ( GetNumParentsTest );
    CPPUNIT_TEST ( GetNameValPairTest );
    CPPUNIT_TEST ( GetNameValPairCntTest );
    CPPUNIT_TEST ( UpdateNameValPairTest );
    CPPUNIT_TEST ( FindNameValPairTest );
		CPPUNIT_TEST ( GetNameValParamsBeginsWithTest );
		CPPUNIT_TEST ( FindParentByFileTypeIdTest );

	CPPUNIT_TEST_SUITE_END();

private:
    
    GenericDataHeader *header;

public:
	void setUp();
	void tearDown();

	void testCreation();
    void FileTypeIdTest();
    void FileIdTest();
    void FileCreationTimeTest();
    void LocaleTest();
    void NameValTest();
    void AddParentEntryTest();
    void GetNumParentsTest();
    void GetNameValPairTest();
    void GetNameValPairCntTest();
    void UpdateNameValPairTest();
    void FindNameValPairTest();
		void GetNameValParamsBeginsWithTest();
		void FindParentByFileTypeIdTest();
};

#endif // __GENERICDATAHEADERTEST_H_
