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

#ifndef __NORMALIZATIONTEST_H_
#define __NORMALIZATIONTEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include <vector>

class NormalizationTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( NormalizationTest );
	//
	CPPUNIT_TEST( test_ExtractSketch       );
	CPPUNIT_TEST( test_InterpolateSketch   );
	CPPUNIT_TEST( test_MedianNormalization );
	CPPUNIT_TEST( test_Normalization       );
        CPPUNIT_TEST( test_BiocCompat          );
        CPPUNIT_TEST( test_InterpolateVsCore   );
        CPPUNIT_TEST( test_DocExample          );
  //
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
  //
	void test_ExtractSketch();
        void test_InterpolateSketch();
	void test_MedianNormalization();
	void test_Normalization();
        void test_BiocCompat();
        void test_InterpolateVsCore();
        bool test_InterpolateVsCoreNorm(vector<vector<float> > inputCoreData,
                                        vector<vector<float> > inputInterpData,
                                        bool doBioc, float maxDiff = .001,
                                        int sketchSize = 0);
        void fillInWRandData(unsigned int seed, int maxVal, int numChips, int numProbes, 
                                   vector<vector<float> > &toFill);
        void test_DocExample();
};

#endif // __NORMALIZATIONTEST_H_
