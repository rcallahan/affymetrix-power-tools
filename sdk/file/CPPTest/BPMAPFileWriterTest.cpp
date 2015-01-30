////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include "file/CPPTest/BPMAPFileWriterTest.h"
//
#include "file/BPMAPFileWriter.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

CPPUNIT_TEST_SUITE_REGISTRATION( CBPMAPFileWriterTest );

#ifdef _MSC_VER
#define TPMAP_V3_IN_FILE "TilingArrays\\test.v3.tpmap"
#define V3_OUT_FILE "TilingArrays\\test.v3.bpmap"
#else
#define TPMAP_V3_IN_FILE "TilingArrays/test.v3.tpmap"
#define V3_OUT_FILE "TilingArrays/test.v3.bpmap"
#endif

using namespace affxbpmap;
using namespace affxbpmapwriter;

static std::string GetDataPath(const char *relPath)
{
//	extern std::string externalDataPath; // The path to the data.
//	std::string path = externalDataPath;
//#ifdef _MSC_VER
//	path += "\\";
//#else
//	path += "/";
//#endif
//	path += relPath;
//	return path;
  return Fs::join(externalDataPath,relPath);
}

void CBPMAPFileWriterTest::setUp()
{
}

void CBPMAPFileWriterTest::tearDown()
{
}

void CBPMAPFileWriterTest::testWrite()
{
	CBPMAPFileWriter bpmap;

	std::string path=GetDataPath(TPMAP_V3_IN_FILE);
	bpmap.SetTpmapFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE("TPMAP file does not exist", bpmap.TpmapExists() == true);

	path=GetDataPath(V3_OUT_FILE);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.ReadTpmap() );
	CPPUNIT_ASSERT( bpmap.WriteBpmap() );
}
