////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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


/// @file   DumpPgfTest.h
/// @brief  Header for DumpPgfTest.cpp

#ifndef DUMP_PGF_TEST_H
#define DUMP_PGF_TEST_H

#include <cppunit/extensions/HelperMacros.h>

class DumpPgfTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DumpPgfTest );
  CPPUNIT_TEST( testDumpEntirePgf );
  CPPUNIT_TEST( testDumpProbesetIds );
  CPPUNIT_TEST( testDumpProbeIds );
  CPPUNIT_TEST( testDumpSequentialClf );
  CPPUNIT_TEST( testDumpType );
  CPPUNIT_TEST( testDumpTypeOr );
  CPPUNIT_TEST_SUITE_END();

public:
  void testDumpEntirePgf();
  void testDumpProbesetIds();
  void testDumpProbeIds();
  void testDumpSequentialClf();
  void testDumpType();
  void testDumpTypeOr();
};

#endif /* DUMP_PGF_TEST_H */
