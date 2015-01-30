////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include "file/CPPTest/GRDFileDataTest.h"
//
#include "file/GRDFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CGRDFileDataTest );

#define TEST_GRD_FILE "GRD/test.grd"
extern std::string externalDataPath; // The path to the data.

static std::string GetDataPath(const char *relPath)
{
//	extern std::string externalDataPath; // The path to the data.
//	std::string path = externalDataPath;
//#ifdef _MSC_VER
//	path += "\\";
//#else
//	path += "/";
//#endif
//	path += relPath;
//	return path;
  return Fs::join(externalDataPath,relPath);
}


void CGRDFileDataTest::setUp()
{
}

void CGRDFileDataTest::tearDown()
{
}

void CGRDFileDataTest::testCreation()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT( 1 );
}

void CGRDFileDataTest::testproperty_FileName()
{
	affxgrd::CGRDFileData GRD;
	std::string path = GetDataPath(TEST_GRD_FILE);
	GRD.SetFileName(path.c_str());
	CPPUNIT_ASSERT( GRD.GetFileName() == path );
}

void CGRDFileDataTest::testmethod_Exists()
{
	affxgrd::CGRDFileData GRD;
	GRD.SetFileName(GetDataPath(TEST_GRD_FILE).c_str());
	CPPUNIT_ASSERT( GRD.Exists() == true );
	GRD.SetFileName("test");
	CPPUNIT_ASSERT( GRD.Exists() == false );
}

void CGRDFileDataTest::testmethod_Read()
{
	affxgrd::CGRDFileData GRD;

	GRD.SetFileName("test");
	CPPUNIT_ASSERT( GRD.Read() == false);

	GRD.SetFileName(GetDataPath(TEST_GRD_FILE).c_str());
	CPPUNIT_ASSERT( GRD.Read() == true);

	double eps = 1e-10;
	double eps3 = 0.001;

	CPPUNIT_ASSERT_DOUBLES_EQUAL(GRD.GetHeader().GetVersion(), 1.0, eps);
	CPPUNIT_ASSERT(GRD.GetHeader().GetCols() == 2560);
	CPPUNIT_ASSERT(GRD.GetHeader().GetRows() == 2560);
	CPPUNIT_ASSERT(GRD.GetHeader().GetRows()*GRD.GetHeader().GetCols() == GRD.GetHeader().GetNumCells());

	CPPUNIT_ASSERT_DOUBLES_EQUAL(GRD.GetHeader().GetFeaturePitchX(), 5.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(GRD.GetHeader().GetFeaturePitchY(), 5.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(GRD.GetHeader().GetFeatureSetbackX(), 2.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(GRD.GetHeader().GetFeatureSetbackY(), 2.0, eps);

	affxgrd::StrStrMap &params = GRD.GetHeader().GetParameters();
	CPPUNIT_ASSERT( params.size() == 3);
	std::string val;

	val = params[affxgrd::SZ_PARENT_DAT_PROP_NAME];
	CPPUNIT_ASSERT( val == "C:\\GeneChip\\Affy_Data\\Data\\Test1.DAT" );

	val = params[affxgrd::SZ_SCAN_DATE_TIME_PROP_NAME];
	CPPUNIT_ASSERT( val == "0 0 0 0 0 0 0 0 0 " );

	val = params[affxgrd::SZ_SCANNER_ID_PROP_NAME];
	CPPUNIT_ASSERT( val == "7G_Feas_3");


	CPPUNIT_ASSERT( GRD.GetHeader().GetNumSubgrids() == 225);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ul.fx, 1030.2562, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ul.fy, 653.97498, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ur.fx, 2742.8, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ur.fy, 676.90002, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ll.fx, 1006.2781, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).ll.fy, 2336.3977, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).lr.fx, 2719.7251, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(0).lr.fy, 2360.2500, eps3);

	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ul.fx, 24237.250, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ul.fy, 24460.137, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ur.fx, 25920.801, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ur.fy, 24481.350, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ll.fx, 24213.0, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).ll.fy, 26171.926, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).lr.fx, 25896.574, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( GRD.GetHeader().GetOptSubgrid(224).lr.fy, 26194.676, eps3);


	affxgrd::FCOORD center;

	center = GRD.GetCenter(0, 0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fx, 1035.1843, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fy, 658.60431, eps3);

	center = GRD.GetCenter(200, 100);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fx, 2990.6855, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fy, 1669.6373, eps3);

	center = GRD.GetCenter(GRD.GetHeader().GetCols()-1, GRD.GetHeader().GetRows()-1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fx, 25892.199, eps3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(center.fy, 26190.176, eps3);
}

void CGRDFileDataTest::VersionTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0f, GRD.GetHeader().GetVersion(), 0.000001);
	GRD.GetHeader().SetVersion(1.1f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.1f, GRD.GetHeader().GetVersion(), 0.000001);
}

void CGRDFileDataTest::RowsTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT(GRD.GetHeader().GetRows() == 0);
	GRD.GetHeader().SetRows(5);
	CPPUNIT_ASSERT(GRD.GetHeader().GetRows() == 5);
}

void CGRDFileDataTest::ColsTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT(GRD.GetHeader().GetCols() == 0);
	GRD.GetHeader().SetCols(4);
	CPPUNIT_ASSERT(GRD.GetHeader().GetCols() == 4);
}

void CGRDFileDataTest::NumCellsTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT(GRD.GetHeader().GetNumCells() == 0);
	GRD.GetHeader().SetRows(3);
	GRD.GetHeader().SetCols(4);
	CPPUNIT_ASSERT(GRD.GetHeader().GetNumCells() == 12);
}

void CGRDFileDataTest::FeaturePitchTest()
{
	double eps = 1e-5;
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, GRD.GetHeader().GetFeaturePitchX(), eps);
	GRD.GetHeader().SetFeaturePitchX(3.4f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4f, GRD.GetHeader().GetFeaturePitchX(), eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, GRD.GetHeader().GetFeaturePitchY(), eps);
	GRD.GetHeader().SetFeaturePitchY(1.3f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3f, GRD.GetHeader().GetFeaturePitchY(), eps);
}

void CGRDFileDataTest::FeatureSetbackTest()
{
	double eps = 1e-5;
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, GRD.GetHeader().GetFeatureSetbackX(), eps);
	GRD.GetHeader().SetFeatureSetbackX(3.4f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4f, GRD.GetHeader().GetFeatureSetbackX(), eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, GRD.GetHeader().GetFeatureSetbackY(), eps);
	GRD.GetHeader().SetFeatureSetbackY(1.3f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3f, GRD.GetHeader().GetFeatureSetbackY(), eps);
}

void CGRDFileDataTest::OptSubgridTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT(GRD.GetHeader().GetNumSubgrids() == 0);
	affxgrd::FRECT r;
	r.ul.fx = 1.1f;
	r.ul.fy = 1.2f;
	r.ur.fx = 10.1f;
	r.ur.fy = 1.3f;
	r.lr.fx = 10.2f;
	r.lr.fy = 10.3f;
	r.ll.fx = 1.4f;
	r.ll.fx = 10.4f;

	affxgrd::FRECT r2 = r;
	r2.ll.fx = 2.3f;

	GRD.GetHeader().AddOptSubgrid(&r);
	GRD.GetHeader().AddOptSubgrid(&r2);

	CPPUNIT_ASSERT(GRD.GetHeader().GetNumSubgrids() == 2);

	affxgrd::FRECT r3 = GRD.GetHeader().GetOptSubgrid(0);
	CPPUNIT_ASSERT(r3.ll.fx == r.ll.fx);
	affxgrd::FRECT r4 = GRD.GetHeader().GetOptSubgrid(1);
	CPPUNIT_ASSERT(r4.ll.fx == r2.ll.fx);
}

void CGRDFileDataTest::ParametersTest()
{
	affxgrd::CGRDFileData GRD;
	CPPUNIT_ASSERT(GRD.GetHeader().GetParameters().size() == 0);

	affxgrd::StrStrMap& params = GRD.GetHeader().GetParameters();
	params["hello"] = "bon jour";
	params["adios"] = "goodbye";

	CPPUNIT_ASSERT(GRD.GetHeader().GetParameters().size() == 2);

	CPPUNIT_ASSERT(params.find("hello")->second == "bon jour");

}
