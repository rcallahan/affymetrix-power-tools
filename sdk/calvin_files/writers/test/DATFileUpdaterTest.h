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

#ifndef __DATFILEUPDATERTEST_H_
#define __DATFILEUPDATERTEST_H_

#include "calvin_files/data/src/DATData.h"
#include "calvin_files/utils/src/Coords.h"
#include "calvin_files/writers/src/DATFileUpdater.h"
#include "calvin_files/writers/src/DATFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class DATFileUpdaterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DATFileUpdaterTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( GlobalGridTest );
	CPPUNIT_TEST ( GlobalGridIITest );
	CPPUNIT_TEST ( GlobalGridIIITest );
	CPPUNIT_TEST ( SubGridTest );
	CPPUNIT_TEST ( SubGridIITest );
	CPPUNIT_TEST ( SubGridIIITest );
	CPPUNIT_TEST ( UpdateFileIdTest );

	CPPUNIT_TEST_SUITE_END();

public:

	void setUp();
	void tearDown();
	void testCreation();
	void GlobalGridTest();
	void GlobalGridIITest();
	void GlobalGridIIITest();
	void SubGridTest();
	void SubGridIITest();
	void SubGridIIITest();
	void UpdateFileIdTest();

private:

	void WriteDatFile(const std::string& filename);

	void WriteGlobalGridDatFile(const std::string& filename, float increment);

	void WriteSubGridDatFile(const std::string& filename, float increment);

	void AddGlobalGridData(DATData& data, float increment);

	void AddSubGridData(DATData& data, float increment);

	void ReadGlobalGridData(const std::string& filename, FRegion& region, u_int32_t& status, ParameterNameValueTypeVector& params);

	void ReadSubGridData(const std::string& filename, FRegionVector& regions, Uint32Vector& status, ParameterNameValueTypeVector& params);

	bool IsEqual(const FRegionVector& r1, const FRegionVector& r2);

	void GetSubGridData(const DATData& data, FRegionVector& regions, Uint32Vector& status);

	AffymetrixGuidType ReadFileId(const std::string& filename);

	void AddParameters(DATData& data);

	void CompareParameterNameValueTypeVectors(ParameterNameValueTypeVector& left, ParameterNameValueTypeVector& right);


};

#endif // __DATFILEUPDATERTEST_H_
