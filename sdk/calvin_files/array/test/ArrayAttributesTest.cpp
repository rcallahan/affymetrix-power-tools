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
#include "calvin_files/array/test/ArrayAttributesTest.h"
//
#include "calvin_files/array/src/ArrayAttributes.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ArrayAttributesTest );

void ArrayAttributesTest::setUp()
{
}

void ArrayAttributesTest::tearDown()
{
}

void ArrayAttributesTest::testCreation()
{
	ArrayAttributes atts;
	CPPUNIT_ASSERT(1);
}

void ArrayAttributesTest::testproperty_Identifier()
{
	string id = "file_id";
	ArrayAttributes atts;
	atts.Identifier() = id;
	CPPUNIT_ASSERT(atts.Identifier() == id);
}

void ArrayAttributesTest::testproperty_ArrayBarcode()
{
	ArrayAttributes atts;
	atts.ArrayBarcode() = "123";
	CPPUNIT_ASSERT(atts.ArrayBarcode() == "123");
}

void ArrayAttributesTest::testproperty_ArrayName()
{
	ArrayAttributes atts;
	atts.ArrayName() = "name";
	CPPUNIT_ASSERT(atts.ArrayName() == "name");
}

void ArrayAttributesTest::testproperty_Comment()
{
	ArrayAttributes atts;
	atts.Comment() = L"comment";
	CPPUNIT_ASSERT(atts.Comment() == L"comment");
}

void ArrayAttributesTest::testproperty_Media()
{
	ArrayAttributes atts;
	atts.Media() = PlateOrStripMedia;
	CPPUNIT_ASSERT(atts.Media() == PlateOrStripMedia);
}

void ArrayAttributesTest::testproperty_CreatedBy()
{
	ArrayAttributes atts;
	atts.CreatedBy() = L"me";
	CPPUNIT_ASSERT(atts.CreatedBy() == L"me");
}

void ArrayAttributesTest::testproperty_CreationDateTime()
{
	ArrayAttributes atts;
	atts.CreationDateTime() = L"now";
	CPPUNIT_ASSERT(atts.CreationDateTime() == L"now");
}

void ArrayAttributesTest::testproperty_MediaFileName()
{
	ArrayAttributes atts;
	atts.MediaFileName() = "MediaFileName";
	CPPUNIT_ASSERT(atts.MediaFileName() == "MediaFileName");
}

void ArrayAttributesTest::testproperty_MediaFileGUID()
{
	ArrayAttributes atts;
	atts.MediaFileGUID() = "MediaFileGUID";
	CPPUNIT_ASSERT(atts.MediaFileGUID() == "MediaFileGUID");
}

void ArrayAttributesTest::testproperty_LibraryPackageName()
{
	ArrayAttributes atts;
	atts.LibraryPackageName() = "LibraryPackageName";
	CPPUNIT_ASSERT(atts.LibraryPackageName() == "LibraryPackageName");
}

void ArrayAttributesTest::testproperty_MasterFile()
{
	ArrayAttributes atts;
	atts.MasterFile() = "master";
	CPPUNIT_ASSERT(atts.MasterFile() == "master");
}

void ArrayAttributesTest::testproperty_MediaCol()
{
	ArrayAttributes atts;
	atts.MediaCol() = 1;
	CPPUNIT_ASSERT(atts.MediaCol() == 1);
}

void ArrayAttributesTest::testproperty_MediaRow()
{
	ArrayAttributes atts;
	atts.MediaRow() = 2;
	CPPUNIT_ASSERT(atts.MediaRow() == 2);
}

void ArrayAttributesTest::testproperty_PatAssignment()
{
	ArrayAttributes atts;
	atts.PatAssignment() = AffyBarcodeAssignment;
	CPPUNIT_ASSERT(atts.PatAssignment() == AffyBarcodeAssignment);
}

void ArrayAttributesTest::testproperty_MasterFileId()
{
	ArrayAttributes atts;
	atts.MasterFileId() = "master_id";
	CPPUNIT_ASSERT(atts.MasterFileId() == "master_id");
}

void ArrayAttributesTest::testmethod_Clear()
{
	string id = "file_id";
	ArrayAttributes atts;
	atts.Identifier() = id;
	atts.Attributes().resize(2);
	atts.ArrayBarcode() = "123";
	atts.ArrayName() = "name";
	atts.Comment() = L"comment";
	atts.Media() = PlateOrStripMedia;
	atts.CreatedBy() = L"me";
	atts.CreationDateTime() = L"now";
	atts.MediaFileName() = "MediaFileName";
	atts.MediaFileGUID() = "MediaFileGUID";
	atts.LibraryPackageName() = "LibraryPackageName";
	atts.MasterFile() = "master";
	atts.MediaCol() = 1;
	atts.MediaRow() = 2;
	atts.PatAssignment() = AffyBarcodeAssignment;
	atts.MasterFileId() = "master_id";
	atts.Clear();

	CPPUNIT_ASSERT(atts.Identifier() == "");
	CPPUNIT_ASSERT(atts.Attributes().size() == 0);
	CPPUNIT_ASSERT(atts.ArrayBarcode() == "");
	CPPUNIT_ASSERT(atts.ArrayName() == "");
	CPPUNIT_ASSERT(atts.Comment() == L"");
	CPPUNIT_ASSERT(atts.Media() == CartridgeMedia);
	CPPUNIT_ASSERT(atts.CreatedBy() == L"");
	CPPUNIT_ASSERT(atts.CreationDateTime() == L"");
	CPPUNIT_ASSERT(atts.MasterFile() == "");
	CPPUNIT_ASSERT(atts.MediaCol() == 0);
	CPPUNIT_ASSERT(atts.MediaRow() == 0);
	CPPUNIT_ASSERT(atts.PatAssignment() == NoAssignment);
	CPPUNIT_ASSERT(atts.MasterFileId() == "");
	CPPUNIT_ASSERT(atts.MediaFileName() == "");
	CPPUNIT_ASSERT(atts.MediaFileGUID() == "");
	CPPUNIT_ASSERT(atts.LibraryPackageName() == "");

}

void ArrayAttributesTest::testproperty_Attributes()
{
	ArrayAttributes atts;
	ParameterNameValuePairVector &params = atts.Attributes();
	ParameterNameValuePair param;

	param.Name = L"array-att-name-1";
	param.Value = L"array-att-value-1";
	params.push_back(param);

	param.Name = L"array-att-name-2";
	param.Value = L"array-att-value-2";
	params.push_back(param);

	CPPUNIT_ASSERT( params.size() == 2 );

	ParameterNameValuePairVector::iterator it = atts.Attributes().begin();
	CPPUNIT_ASSERT( (*it).Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( (*it).Value == L"array-att-value-1" );
	++it;
	CPPUNIT_ASSERT( (*it).Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( (*it).Value == L"array-att-value-2" );
	++it;
	CPPUNIT_ASSERT (it == atts.Attributes().end() );

	CPPUNIT_ASSERT( atts.Attributes()[0].Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( atts.Attributes()[0].Value == L"array-att-value-1" );
	CPPUNIT_ASSERT( atts.Attributes()[1].Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( atts.Attributes()[1].Value == L"array-att-value-2" );
}
