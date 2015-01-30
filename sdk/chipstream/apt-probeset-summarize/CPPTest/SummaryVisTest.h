////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/// @file   SummaryVisTest.h
/// @brief  Header for SummaryVisTest.cpp

#ifndef SUMMARY_VIS_TEST_H
#define SUMMARY_VIS_TEST_H

#include <cppunit/extensions/HelperMacros.h>

class SummaryVisTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SummaryVisTest );
  CPPUNIT_TEST( testSummaryVis );
  CPPUNIT_TEST( testWiggleColName );
  CPPUNIT_TEST( testWiggleColIndex );
  CPPUNIT_TEST_SUITE_END();

public:
  void testSummaryVis();
  void testWiggleColName();
  void testWiggleColIndex();
};

#endif /* SUMMARY_VIS_TEST_H */
