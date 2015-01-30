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

/// @file   CelExtractTest.h
/// @brief  Header for CelExtractTest.cpp

#ifndef CEL_EXTRACT_TEST_H
#define CEL_EXTRACT_TEST_H

#include <cppunit/extensions/HelperMacros.h>

class CelExtractTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( CelExtractTest );
  CPPUNIT_TEST( testExtractEntirePgf );
  CPPUNIT_TEST( testExtractProbesetIds );
  CPPUNIT_TEST( testExtractProbeIds );
  CPPUNIT_TEST( testExtractPmGcbg );
  CPPUNIT_TEST( testExtractPmMm );
  CPPUNIT_TEST( testExtractCdf );
  CPPUNIT_TEST_SUITE_END();

public:
  void testExtractEntirePgf();
  void testExtractProbesetIds();
  void testExtractProbeIds();
  void testExtractPmGcbg();
  void testExtractPmMm();
  void testExtractCdf();
};

#endif /* CEL_EXTRACT_TEST_H */
