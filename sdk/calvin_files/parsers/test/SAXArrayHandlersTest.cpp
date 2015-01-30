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
#include "calvin_files/parsers/test/SAXArrayHandlersTest.h"
//
#include "calvin_files/parsers/src/SAXArrayHandlers.h"
//
#include <cstring>
#include <string.h>
//

using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( SAXArrayHandlersTest );

void SAXArrayHandlersTest::setUp()
{
}

void SAXArrayHandlersTest::tearDown()
{
}

void SAXArrayHandlersTest::testCreation_startDocument()
{
	ArrayData array;
	SAXArrayHandlers handler(&array);
	handler.startDocument();
	CPPUNIT_ASSERT(1);
}

void SAXArrayHandlersTest::testCreation_endDocument()
{
	ArrayData array;
	SAXArrayHandlers handler(&array);
	handler.endDocument();
	CPPUNIT_ASSERT(1);
}


/*! Converts a XML character string to wide string.
 * @param c1 The XML string to convert
 * @return The wide character string
 */
static const XMLCh * StringToXMLCh_ASCIIONLY(const char * const c1)
{
	int n = (int)strlen(c1);
	XMLCh *s = new XMLCh[n+1];
	s[n] = 0;
	for (int i=0; i<n; i++)
	{
		s[i] = c1[i];
	}
	return s;
}


void SAXArrayHandlersTest::testCreation_startElement_with_bad_data()
{
	class TestAttributes : public XERCES_CPP_NAMESPACE::AttributeList
	{
	public:
		XMLSize_t getLength() const { return 2; }
        const XMLCh* getName(const XMLSize_t) const { return NULL; }
        const XMLCh *getType(const XMLCh *const ) const { return NULL; }
        const XMLCh *getType(const XMLSize_t) const { return NULL; }
        const XMLCh *getValue(const char *const ) const { return NULL; }
        const XMLCh *getValue(const XMLCh *const ) const { return NULL; }
        const XMLCh *getValue(const XMLSize_t) const { return NULL; }
	};
	ArrayData array;
	SAXArrayHandlers handler(&array);
	TestAttributes attributes;
	handler.startElement(NULL, attributes);
	const XMLCh *c = StringToXMLCh_ASCIIONLY("invalid_name");
	handler.startElement(c, attributes);
	delete [] c;
	CPPUNIT_ASSERT(1);
}
