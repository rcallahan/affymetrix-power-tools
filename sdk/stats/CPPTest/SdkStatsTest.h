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

#ifndef __SDKSTATSTEST_H_
#define __SDKSTATSTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class SdkStatsTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( SdkStatsTest );
	
  CPPUNIT_TEST(test_pnorm);
  CPPUNIT_TEST(test_psignrank);
  CPPUNIT_TEST(test_signedRankTest);
  CPPUNIT_TEST(test_pwilcox);
  CPPUNIT_TEST(test_ranksumTest);
  CPPUNIT_TEST(test_uprob);
  CPPUNIT_TEST(test_chisqrprob);
  CPPUNIT_TEST(test_PearsonCorrelation);
  CPPUNIT_TEST(test_CalcHWEqPValue);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void test_pnorm();
  void test_psignrank();
  void test_signedRankTest();
  void test_pwilcox();
  void test_ranksumTest();
  void test_uprob();
  void test_chisqrprob();
  void test_PearsonCorrelation();
  void test_CalcHWEqPValue();
};

#endif // __SDKSTATSTEST_H_
