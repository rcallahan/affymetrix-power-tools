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
#include "file/CPPTest/DttArrayFileReaderTest.h"
//
#include "file/DttArrayFileReader.h"
//

using namespace std;
using namespace affymetrix_dttarray;

CPPUNIT_TEST_SUITE_REGISTRATION( DttArrayFileReaderTest );


#define WITH_TEMPLATES "./data/sample_and_exp_template.xml"
#define NO_FILE "no_file.xml"
#define INVALID_FILE "./data/test.1lq"
#define SAMPLE_FILE "./data/gdac-sample.xml"
#define EXP_FILE "./data/gdac-exp.xml"
#define SCAN_FILE "./data/gdac-scan.xml"

void DttArrayFileReaderTest::setUp()
{
}

void DttArrayFileReaderTest::tearDown()
{
}

void DttArrayFileReaderTest::testproperty_ArrayType()
{
	const char *type = "array type";
	DttArrayData data;
	data.SetArrayType(type);
	CPPUNIT_ASSERT( data.GetArrayType() == type );
}

void DttArrayFileReaderTest::testproperty_ExpName()
{
	const char *exp = "experiment name";
	DttArrayData data;
	data.SetExperimentName(exp);
	CPPUNIT_ASSERT( data.GetExperimentName() == exp);
}

void DttArrayFileReaderTest::testproperty_Attributes()
{
	DttArrayData data;
	AttributeNameValueType att;
	att.name = "name1";
	att.value = "value1";
	att.type = "type1";
	data.Attributes().push_back(att);
	att.name = "name2";
	att.value = "value2";
	att.type = "type2";
	data.Attributes().push_back(att);
	CPPUNIT_ASSERT( data.Attributes().size() == 2);

	AttributeNameValueTypeList::iterator it = data.Attributes().begin();

	att = *it;
	CPPUNIT_ASSERT( att.name == "name1");
	CPPUNIT_ASSERT( att.value == "value1");
	CPPUNIT_ASSERT( att.type == "type1");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "name2");
	CPPUNIT_ASSERT( att.value == "value2");
	CPPUNIT_ASSERT( att.type == "type2");
}

void DttArrayFileReaderTest::testmethod_Clear()
{
	const char *exp = "experiment name";
	const char *type = "array type";
	DttArrayData data;
	data.SetArrayType(type);
	data.SetExperimentName(exp);
	AttributeNameValueType att;
	att.name = "name1";
	att.value = "value1";
	att.type = "type1";
	data.Attributes().push_back(att);
	data.Clear();
	CPPUNIT_ASSERT( data.GetArrayType() == "");
	CPPUNIT_ASSERT( data.GetExperimentName() == "");
	CPPUNIT_ASSERT( data.Attributes().size() == 0);
}

void DttArrayFileReaderTest::testmethod_Read()
{
	DttArrayData data;
	DttArrayFileReader reader;
	reader.SetFileName(WITH_TEMPLATES);
	CPPUNIT_ASSERT(reader.Read(data) == true);
	CPPUNIT_ASSERT(data.GetArrayType() == "Test3");
	CPPUNIT_ASSERT(data.GetExperimentName() == "Example2_No_Wash_Stain_Uses_Sample_And_Exp_Template");
	CPPUNIT_ASSERT(data.Attributes().size() == 33);

	AttributeNameValueType att;
	AttributeNameValueTypeList::iterator it = data.Attributes().begin();

	att = *it;
	CPPUNIT_ASSERT( att.name == "GCOS Sample Name");
	CPPUNIT_ASSERT( att.value == "Sample_with_Sample_Template");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Sample Template Name");
	CPPUNIT_ASSERT( att.value == "MIAME Sample Information");
	CPPUNIT_ASSERT( att.type == "string");

	it = data.Attributes().end();
	--it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Grid Alignment Passed");
	CPPUNIT_ASSERT( att.value == "");
	CPPUNIT_ASSERT( att.type == "controlled");
}

void DttArrayFileReaderTest::testmethod_Read_when_file_does_not_exist()
{
	DttArrayData data;
	DttArrayFileReader reader;
	reader.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(reader.Read(data) == false);
}

void DttArrayFileReaderTest::testmethod_Read_when_file_is_not_valid()
{
	DttArrayData data;
	DttArrayFileReader reader;
	reader.SetFileName(INVALID_FILE);
	CPPUNIT_ASSERT(reader.Read(data) == false);
}

void DttArrayFileReaderTest::testmethod_Exists_pass()
{
	DttArrayFileReader reader;
	reader.SetFileName(WITH_TEMPLATES);
	CPPUNIT_ASSERT(reader.Exists() == true);
}

void DttArrayFileReaderTest::testmethod_Exists_fail()
{
	DttArrayFileReader reader;
	reader.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(reader.Exists() == false);
}

void DttArrayFileReaderTest::testmethod_Read_from_gdac_exporter_sample()
{
	DttArrayData data;
	DttArrayFileReader reader;
	reader.SetFileName(SAMPLE_FILE);
	CPPUNIT_ASSERT(reader.Read(data) == true);
	CPPUNIT_ASSERT( data.GetArrayType() == "");
	CPPUNIT_ASSERT( data.GetExperimentName() == "");
	CPPUNIT_ASSERT( data.Attributes().size() == 8);

	AttributeNameValueType att;
	AttributeNameValueTypeList::iterator it = data.Attributes().begin();

	att = *it;
	CPPUNIT_ASSERT( att.name == "GCOS Sample Name");
	CPPUNIT_ASSERT( att.value == "s_with_temp");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Sample Template Name");
	CPPUNIT_ASSERT( att.value == "samp_temp");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Tissue Type");
	CPPUNIT_ASSERT( att.value == "Skin");
	CPPUNIT_ASSERT( att.type == "controlled");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Date");
	CPPUNIT_ASSERT( att.value == "8-1-2002");
	CPPUNIT_ASSERT( att.type == "date");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Time");
	CPPUNIT_ASSERT( att.value == "12.30 PM");
	CPPUNIT_ASSERT( att.type == "time");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Barcode");
	CPPUNIT_ASSERT( att.value == "123123123");
	CPPUNIT_ASSERT( att.type == "integer");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "GCOS Sample Project");
	CPPUNIT_ASSERT( att.value == "s_proj");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "GCOS Sample Type");
	CPPUNIT_ASSERT( att.value == "s_type");
	CPPUNIT_ASSERT( att.type == "string");

}


void DttArrayFileReaderTest::testmethod_Read_from_gdac_exporter_exp()
{
	DttArrayData data;
	DttArrayFileReader reader;
	reader.SetFileName(EXP_FILE);
	CPPUNIT_ASSERT(reader.Read(data) == true);
	CPPUNIT_ASSERT( data.GetArrayType() == "Test3");
	CPPUNIT_ASSERT( data.GetExperimentName() == "e_with_temp");
	CPPUNIT_ASSERT( data.Attributes().size() == 5);

	AttributeNameValueType att;
	AttributeNameValueTypeList::iterator it = data.Attributes().begin();

	att = *it;
	CPPUNIT_ASSERT( att.name == "Experiment Template Name");
	CPPUNIT_ASSERT( att.value == "exp_temp");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "PI");
	CPPUNIT_ASSERT( att.value == "LJ");
	CPPUNIT_ASSERT( att.type == "controlled");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Desc");
	CPPUNIT_ASSERT( att.value == "GDAC Exporter test case");
	CPPUNIT_ASSERT( att.type == "string");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Date");
	CPPUNIT_ASSERT( att.value == "8-01-2002");
	CPPUNIT_ASSERT( att.type == "date");

	++it;
	att = *it;
	CPPUNIT_ASSERT( att.name == "Time");
	CPPUNIT_ASSERT( att.value == "1.45 AM");
	CPPUNIT_ASSERT( att.type == "time");
}
