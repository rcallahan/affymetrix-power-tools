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

//
#include "calvin_files/utils/test/AffymetrixGuidTest.h"
//
#include "calvin_files/utils/src/AffymetrixGuid.h"
//
#include <map>
//

using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( AffymetrixGuidTest );

void AffymetrixGuidTest::setUp()
{
}

void AffymetrixGuidTest::tearDown()
{
}

void AffymetrixGuidTest::testCreation()
{
	AffymetrixGuid guid;
	CPPUNIT_ASSERT(1);
}

void AffymetrixGuidTest::testmethod_GenerateNewGuid()
{
	AffymetrixGuidType guid;
	std::map<AffymetrixGuidType, bool> guids;
	int nguids = 10000;
	size_t ncount;
	AffymetrixGuid guidgen;
	for (int i=0; i<nguids; i++)
	{
		guid = guidgen.GenerateNewGuid();
		CPPUNIT_ASSERT ( guid.length() == 54 );
		ncount = guids.count(guid);
		if (ncount == 0)
		{
			guids[guid] = true;
		}
		CPPUNIT_ASSERT( ncount == 0);
	}
}
