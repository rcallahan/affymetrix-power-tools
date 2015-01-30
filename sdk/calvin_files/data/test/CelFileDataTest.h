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

#ifndef __CELFILEDATATEST_H_
#define __CELFILEDATATEST_H_

#include "calvin_files/data/src/CELData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class CelFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CelFileDataTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( FilenameTest );
	CPPUNIT_TEST ( SetIntensityCountTest );
	CPPUNIT_TEST ( SetStdDevCountTest );
	CPPUNIT_TEST ( SetPixelCountTest );
	CPPUNIT_TEST ( SetOutlierCountTest );
	CPPUNIT_TEST ( SetMaskCountTest );
	CPPUNIT_TEST (RowsTest);
	CPPUNIT_TEST (ColsTest);
	CPPUNIT_TEST (ArrayTypeTest);
	CPPUNIT_TEST (MasterFileTest);
	CPPUNIT_TEST (LibraryPackageTest);
	CPPUNIT_TEST (AlgorithmNameTest);
	CPPUNIT_TEST (AddAlgorithmParametersTest);
	CPPUNIT_TEST (GetAlgorithmParametersTest);
	CPPUNIT_TEST (FindAlgorithmParameterTest);

	CPPUNIT_TEST_SUITE_END();

private:
    
	CelFileData *cfMeta;

public:
	void setUp();
	void tearDown();

	void testCreation();
	void FilenameTest();
	void SetIntensityCountTest();
	void SetStdDevCountTest();
	void SetPixelCountTest();
	void SetOutlierCountTest();
	void SetMaskCountTest();
	void RowsTest();
	void ColsTest();
	void ArrayTypeTest();
    void MasterFileTest();
	void LibraryPackageTest();
	void AlgorithmNameTest();
	void AddAlgorithmParametersTest();
	void GetAlgorithmParametersTest();
	void FindAlgorithmParameterTest();
};

#endif // __CELFILEDATATEST_H_
