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

//
#include "calvin_files/parsers/test/FileExceptionTest.h"
//
#include "calvin_files/parsers/src/FileException.h"
//

using namespace affymetrix_calvin_exceptions;

CPPUNIT_TEST_SUITE_REGISTRATION( FileExceptionTest );

void FileExceptionTest::setUp()
{
}

void FileExceptionTest::tearDown()
{
}

void FileExceptionTest::testCreation_FileNotFoundException()
{
	FileNotFoundException ex(L"Calvin",L"Default Description, Please Update!",affymetrix_calvin_utilities::DateTime::GetCurrentDateTime().ToString(),std::string(__FILE__),(u_int16_t)__LINE__,0);
	CPPUNIT_ASSERT(1);
}

void FileExceptionTest::testCreation_InvalidVersionException()
{
	InvalidVersionException ex(L"Calvin",L"Default Description, Please Update!",affymetrix_calvin_utilities::DateTime::GetCurrentDateTime().ToString(),std::string(__FILE__),(u_int16_t)__LINE__,0);
	CPPUNIT_ASSERT(1);
}

void FileExceptionTest::testCreation_InvalidFileTypeException()
{
	InvalidFileTypeException ex(L"Calvin",L"Default Description, Please Update!",affymetrix_calvin_utilities::DateTime::GetCurrentDateTime().ToString(),std::string(__FILE__),(u_int16_t)__LINE__,0);
	CPPUNIT_ASSERT(1);
}
