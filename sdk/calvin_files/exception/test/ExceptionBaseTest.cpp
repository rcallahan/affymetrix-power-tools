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

#include "calvin_files/exception/test/ExceptionBaseTest.h"
//
#include "calvin_files/exception/src/ExceptionBase.h"
//

using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ExceptionTest );

void ExceptionTest::setUp()
{
}

void ExceptionTest::tearDown()
{
}

void ExceptionTest::testCreate()
{
	affymetrix_calvin_exceptions::CalvinException* ex = new CalvinException();

	CPPUNIT_ASSERT( ex != NULL );

	delete ex;
	ex = NULL;

	ex = new CalvinException(L"Calvin",L"This is a CPPUnitTest",DateTime::GetCurrentDateTime().ToString(),std::string(__FILE__),(u_int16_t)__LINE__,1234);
	CPPUNIT_ASSERT( ex != NULL );
	delete ex;

	std::wstring datetime = DateTime::GetCurrentDateTime().ToString();
	std::string file = __FILE__;
	u_int16_t line = (u_int16_t)__LINE__;
	u_int64_t code = 1234;
	std::wstring source = L"Calvin";
	std::wstring desc = L"This tests the throw";

	try
	{
		CalvinException exception(source,desc,datetime,file,line,code);
		throw exception;
	}
	catch(CalvinException& cex)
	{
		CPPUNIT_ASSERT(cex.Source().compare(source) == 0);
		CPPUNIT_ASSERT(cex.Description().compare(desc) == 0);
		CPPUNIT_ASSERT(cex.TimeStamp().compare(datetime) == 0);
		CPPUNIT_ASSERT(cex.ErrorCode() == code);
		CPPUNIT_ASSERT(cex.LineNumber() == line);
		CPPUNIT_ASSERT(cex.SourceFile().compare(file) == 0);
	}
}
void ExceptionTest::testToString()
{
	affymetrix_calvin_exceptions::CalvinException ex;
	CPPUNIT_ASSERT( ex.ToString().compare(L"Not implemented yet.") == 0 );
}
