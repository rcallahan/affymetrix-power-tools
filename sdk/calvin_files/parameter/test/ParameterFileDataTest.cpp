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
#include "calvin_files/parameter/src/ParameterFileData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_parameter;

class ParameterFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (ParameterFileDataTest);

	CPPUNIT_TEST (testAccess);

	CPPUNIT_TEST_SUITE_END();
public:
	void setUp();
	void tearDown();

    void testAccess();
};


CPPUNIT_TEST_SUITE_REGISTRATION( ParameterFileDataTest );

void ParameterFileDataTest::setUp()
{
}

void ParameterFileDataTest::tearDown()
{
}

void ParameterFileDataTest::testAccess()
{
	ParameterFileData d;

    d.ImplementationAttributes().description = L"d";
    d.ImplementationAttributes().executableFileName = L"e";
    d.ImplementationAttributes().name = L"n";
    d.ImplementationAttributes().version = L"v";

    CPPUNIT_ASSERT(d.ImplementationAttributes().description == L"d");
    CPPUNIT_ASSERT(d.ImplementationAttributes().executableFileName == L"e");
    CPPUNIT_ASSERT(d.ImplementationAttributes().name == L"n");
    CPPUNIT_ASSERT(d.ImplementationAttributes().version == L"v");

    d.ParameterFileAttributes().company = L"c";
    d.ParameterFileAttributes().userName = L"u";
    d.ParameterFileAttributes().contentVersion = L"cv";

    CPPUNIT_ASSERT(d.ParameterFileAttributes().company == L"c");
    CPPUNIT_ASSERT(d.ParameterFileAttributes().userName == L"u");
    CPPUNIT_ASSERT(d.ParameterFileAttributes().contentVersion == L"cv");

    ParameterType p;

    p.name = L"n";
    p.index = L"i";
    p.displayName = L"d";
    p.category = L"c";
    p.isEditable = L"is";
    p.type = L"t";
    p.currentValue = L"cv";
    p.minValue = L"mv";
    p.maxValue = L"xv";
    p.defaultValue = L"dv";
    p.precision = L"p";
    p.maxLength = L"l";
    p.description = L"de";

    CPPUNIT_ASSERT(p.name == L"n");
    CPPUNIT_ASSERT(p.index == L"i");
    CPPUNIT_ASSERT(p.displayName == L"d");
    CPPUNIT_ASSERT(p.category == L"c");
    CPPUNIT_ASSERT(p.isEditable == L"is");
    CPPUNIT_ASSERT(p.type == L"t");
    CPPUNIT_ASSERT(p.currentValue == L"cv");
    CPPUNIT_ASSERT(p.minValue == L"mv");
    CPPUNIT_ASSERT(p.maxValue == L"xv");
    CPPUNIT_ASSERT(p.defaultValue == L"dv");
    CPPUNIT_ASSERT(p.precision == L"p");
    CPPUNIT_ASSERT(p.maxLength == L"l");
    CPPUNIT_ASSERT(p.description == L"de");

    d.Parameters().push_back(p);

    p.name.clear();
    p.index.clear();

    CPPUNIT_ASSERT(d.Parameters().size() == 1);

    d.Clear();

    CPPUNIT_ASSERT(d.ImplementationAttributes().description == L"");
    CPPUNIT_ASSERT(d.ImplementationAttributes().executableFileName == L"");
    CPPUNIT_ASSERT(d.ImplementationAttributes().name == L"");
    CPPUNIT_ASSERT(d.ImplementationAttributes().version == L"");

    CPPUNIT_ASSERT(d.ParameterFileAttributes().company == L"");
    CPPUNIT_ASSERT(d.ParameterFileAttributes().userName == L"");
    CPPUNIT_ASSERT(d.ParameterFileAttributes().contentVersion == L"");

    CPPUNIT_ASSERT(d.Parameters().size() == 0);

}
