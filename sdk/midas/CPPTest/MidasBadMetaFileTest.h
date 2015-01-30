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

/// @file   MidasBadMetaFileTest.h
/// @brief  Header for MidasBadMetaFileTest.cpp

#ifndef __MIDASTESTBADMETAFILE_H_
#define __MIDASTESTBADMETAFILE_H_

#include "util/PgOptions.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class MidasBadMetaFileTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MidasBadMetaFileTest );
  CPPUNIT_TEST( testBadMeta );
  CPPUNIT_TEST_SUITE_END();

public:
  void testBadMeta();
};

#endif // __MIDASTESTBADMETAFILE_H_
