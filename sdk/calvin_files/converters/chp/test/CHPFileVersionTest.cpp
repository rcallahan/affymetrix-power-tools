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
#include "calvin_files/converters/chp/test/CHPFileVersionTest.h"
//
#include "calvin_files/converters/chp/src/CHPFileVersion.h"
//
#include <cstdio>
#include <fstream>
//

using namespace std;
using namespace affymetrix_chp_converter;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPFileVersionTest );

void CHPFileVersionTest::setUp()
{
}

void CHPFileVersionTest::tearDown()
{
}

void CHPFileVersionTest::test_DetermineCHPFileVersion()
{
	CHPFileVersionType ver;

	ver = CHPFileVersion::DetermineCHPFileVersion("../data/orig-XDA-stat.CHP");
	CPPUNIT_ASSERT(ver == GCOS_XDA_Version);

	ver = CHPFileVersion::DetermineCHPFileVersion("../data/orig-Calvin-stat.CHP");
	CPPUNIT_ASSERT(ver == Calvin_Version1);

	ver = CHPFileVersion::DetermineCHPFileVersion("./CHPFileVersionTest.cpp");
	CPPUNIT_ASSERT(ver == Unknown_Version);
}
