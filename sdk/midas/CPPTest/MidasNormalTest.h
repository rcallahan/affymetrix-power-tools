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

/// @file   MidasNormalTest.h
/// @brief  Header for MidasNormalTest.cpp

#ifndef __MIDASTEST_H_
#define __MIDASTEST_H_

#include "util/PgOptions.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
#include <vector> 
using namespace std; //#include <vector> 
//
using namespace std; //#include <vector> 
//
using namespace std; //#include <vector> 
//

class MidasNormalTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MidasNormalTest );
  CPPUNIT_TEST( testSpliceDetector );
  CPPUNIT_TEST( testConfigureRun );
  CPPUNIT_TEST( testKeepPath );
  //CPPUNIT_TEST( testKeepWinPath );
  CPPUNIT_TEST( testKeepPathKeep );
  CPPUNIT_TEST_SUITE_END();

public:
  void testSpliceDetector();
  void testConfigureRun();
  void testKeepPath();
  void testKeepWinPath();
  void testKeepPathKeep();
  bool floatsCloseEnough (const float& f1, const float& f2, int digits = 6);
  bool vectorsCloseEnough (const vector<float>& vf1, const vector<float>& vf2, int digits = 6);
  bool vectorsOfVectorsCloseEnough (const vector<vector<float> >& vvf1,
                                    const vector<vector<float> >& vvf2, int digits = 6);
  void resetOptionValues (PgOpt *options[]);
};

#endif // __MIDASTEST_H_
