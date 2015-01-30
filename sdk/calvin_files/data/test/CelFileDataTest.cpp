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
#include "calvin_files/data/test/CelFileDataTest.h"
//
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/utils/src/AffyStlCollectionTypes.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CelFileDataTest );

void CelFileDataTest::setUp()
{
    cfMeta = new CelFileData("CelFileMetaData");
}

void CelFileDataTest::tearDown()
{
    delete cfMeta;
}

void CelFileDataTest::testCreation()
{
	CelFileData hdr("CelFileMetaData");
	CPPUNIT_ASSERT(1);
}

void CelFileDataTest::FilenameTest()
{
    std::string p1 = "CelFileMetaData";
	cfMeta->SetFilename(p1);
    std::string p2 = cfMeta->GetFilename();
    CPPUNIT_ASSERT(p1 == p2);
}

void CelFileDataTest::SetIntensityCountTest()
{
    u_int32_t p1 = 100;
	cfMeta->SetIntensityCount(p1);
    CPPUNIT_ASSERT(1);
}

void CelFileDataTest::SetStdDevCountTest()
{
    u_int32_t p1 = 100;
	cfMeta->SetStdDevCount(p1);
    CPPUNIT_ASSERT(1);
}

void CelFileDataTest::SetPixelCountTest()
{
    u_int32_t p1 = 100;
	cfMeta->SetPixelCount(p1);
    CPPUNIT_ASSERT(1);
}

void CelFileDataTest::SetOutlierCountTest()
{
    u_int32_t p1 = 100;
	cfMeta->SetOutlierCount(p1);
    CPPUNIT_ASSERT(1);
}

void CelFileDataTest::SetMaskCountTest()
{
    u_int32_t p1 = 100;
	cfMeta->SetMaskCount(p1);
    CPPUNIT_ASSERT(1);
}

void CelFileDataTest::RowsTest()
{
	cfMeta->SetRows(45);
	CPPUNIT_ASSERT(cfMeta->GetRows() == 45);
	cfMeta->SetRows(66);
	CPPUNIT_ASSERT(cfMeta->GetRows() == 66);

}

void CelFileDataTest::ColsTest()
{
	cfMeta->SetCols(54);
	CPPUNIT_ASSERT(cfMeta->GetCols() == 54);
	cfMeta->SetCols(112);
	CPPUNIT_ASSERT(cfMeta->GetCols() == 112);
}

void CelFileDataTest::ArrayTypeTest()
{
	cfMeta->SetArrayType(L"MyArray");
	CPPUNIT_ASSERT(cfMeta->GetArrayType() == L"MyArray");
}

void CelFileDataTest::MasterFileTest()
{
    cfMeta->SetMasterFileName(L"masterfile");
    CPPUNIT_ASSERT(cfMeta->GetMasterFileName() == L"masterfile");
}

void CelFileDataTest::LibraryPackageTest()
{
    cfMeta->SetLibraryPackageName(L"libpackage");
    CPPUNIT_ASSERT(cfMeta->GetLibraryPackageName() == L"libpackage");
}

void CelFileDataTest::AlgorithmNameTest()
{
	cfMeta->SetAlgorithmName(L"minimize error");
	CPPUNIT_ASSERT(cfMeta->GetAlgorithmName() == L"minimize error");
}

void CelFileDataTest::AddAlgorithmParametersTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"alpha1");
	nvt.SetValueFloat(0.004f);
	CPPUNIT_ASSERT_NO_THROW(cfMeta->AddAlgorithmParameter(nvt));
	
	// Check that the parameter
	ParameterNameValueTypeIt begin, end;
	cfMeta->GetFileHeader()->GetGenericDataHdr()->GetNameValIterators(begin, end);
	std::wstring parameterName = ALGORITHM_PARAM_NAME_PREFIX;
	parameterName += L"alpha1";
	CPPUNIT_ASSERT(begin->GetName() == parameterName);
}

void CelFileDataTest::GetAlgorithmParametersTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"alpha1");
	nvt.SetValueFloat(0.004f);
	CPPUNIT_ASSERT_NO_THROW(cfMeta->AddAlgorithmParameter(nvt));
	nvt.SetName(L"alpha2");
	nvt.SetValueFloat(0.06f);
	CPPUNIT_ASSERT_NO_THROW(cfMeta->AddAlgorithmParameter(nvt));

	ParameterNameValueTypeVector algParams;
	CPPUNIT_ASSERT_NO_THROW(cfMeta->GetAlgorithmParameters(algParams));
	CPPUNIT_ASSERT(algParams.size() == 2);
	CPPUNIT_ASSERT(algParams.at(0).GetName() == L"alpha1");
	CPPUNIT_ASSERT(algParams.at(1).GetName() == L"alpha2");
}

void CelFileDataTest::FindAlgorithmParameterTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"alpha1");
	nvt.SetValueFloat(0.004f);
	CPPUNIT_ASSERT_NO_THROW(cfMeta->AddAlgorithmParameter(nvt));
	nvt.SetName(L"alpha2");
	nvt.SetValueFloat(0.06f);
	CPPUNIT_ASSERT_NO_THROW(cfMeta->AddAlgorithmParameter(nvt));

	ParameterNameValueType param;
	CPPUNIT_ASSERT(cfMeta->FindAlgorithmParameter(L"alpha2", param));
	CPPUNIT_ASSERT(param.GetValueFloat() == 0.06f);
	CPPUNIT_ASSERT(cfMeta->FindAlgorithmParameter(L"bogus", param) == false);
}
