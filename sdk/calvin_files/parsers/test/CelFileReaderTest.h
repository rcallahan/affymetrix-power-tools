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

#ifndef __CELFILEREADERTEST_H_
#define __CELFILEREADERTEST_H_

#include "calvin_files/data/src/CELData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class CelFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CelFileReaderTest );

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (ReadSmallCelFileNoOutlierNoMaskTest);
	CPPUNIT_TEST (ReadSmallCelFileTest);
	CPPUNIT_TEST (ReadSmallCelFileNoStdevTest);
	//CPPUNIT_TEST (ReadLargeCelFileCheckHeaderTest);
	//CPPUNIT_TEST (ReadLargeCelFileSingleDataAccessTest);
	//CPPUNIT_TEST (ReadLargeCelFileVectorDataAccessTest);
	CPPUNIT_TEST (ReadSmallCelFileGetOutlierCoordsTest);
	CPPUNIT_TEST (ReadSmallCelFileGetMaskedCoordsTest);
	CPPUNIT_TEST (ReadSmallUInt16CelFile);
	CPPUNIT_TEST (ReadMultiChannelCelFileTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void ReadSmallCelFileNoOutlierNoMaskTest();
	void ReadSmallCelFileTest();
	void ReadSmallCelFileNoStdevTest();
	void ReadLargeCelFileCheckHeaderTest();
	void ReadLargeCelFileSingleDataAccessTest();
	void ReadLargeCelFileVectorDataAccessTest();
	void ReadSmallCelFileGetOutlierCoordsTest();
	void ReadSmallCelFileGetMaskedCoordsTest();
	void ReadSmallUInt16CelFile();
	void ReadMultiChannelCelFileTest();

// support methods
	void CheckOutlier(int32_t cell, bool outlier);
	void CheckMasked(int32_t cell, bool masked);
	void CheckSmallCelFileHeader(affymetrix_calvin_io::CelFileData& data);

};

#endif // __CELFILEREADERTEST_H_
