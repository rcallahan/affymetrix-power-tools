////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////


#include "mas5-stat/src/IntensityFileType.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_fusion_io;

class IntensityFileTypeTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( IntensityFileTypeTest );

	CPPUNIT_TEST( testFromHP );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
    void testFromHP();
};

CPPUNIT_TEST_SUITE_REGISTRATION( IntensityFileTypeTest );

void IntensityFileTypeTest::setUp()
{
}

void IntensityFileTypeTest::tearDown()
{
}

void IntensityFileTypeTest::testFromHP()
{
    FusionCELData cel;
    cel.SetFileName("../data/T1_r1.CEL");
    cel.Read(false);
    CPPUNIT_ASSERT(FromHP(cel) == true);
    cel.SetFileName("../data/ht.cel");
    cel.Read(false);
    CPPUNIT_ASSERT(FromHP(cel) == false);
    cel.SetFileName("../data/gcos_m10.CEL");
    cel.Read(false);
    CPPUNIT_ASSERT(FromHP(cel) == false);
    cel.SetFileName("../data/hp.cel");
    cel.Read(false);
    CPPUNIT_ASSERT(FromHP(cel) == true);

}
