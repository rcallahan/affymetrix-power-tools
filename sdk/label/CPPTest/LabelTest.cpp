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

//
#include "label/snp.label.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;

double hardTest[] = {-0.297,-0.284,0.7475,0.764,0.787,-0.409,-0.300,-0.810,
		-0.373,-0.245,-0.429,-0.780,-0.278,0.149,-0.803,0.770,
		-0.793,0.752,0.728,-0.375,0.745,0.199,-0.408,0.756};


/**
 * @class LabelTest
 * @brief cppunit class for testing conversion functions.
 */
class LabelTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( LabelTest );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( HardTests );
  CPPUNIT_TEST_SUITE_END();

public:
// three pretend data sets
	vector<double> onecluster;
	vector<double> twocluster;
	vector<double> threecluster;
	vector<int> Hints;
	vector<double> SubsetInbred;

// reference data
	vector<int> oneref;
	vector<int> tworef;
	vector<int> threeref;

  void setUp();
	void testLabel();
	void HardTests();
	static bool closeEnough(double d1, double d2, int digits=6);

};

bool LabelTest::closeEnough(double d1, double d2, int digits) {
  double diff = fabs(d1 - d2);
  if(diff < 1 / (pow((double)10, digits)))
    return true;
  return false;
}

/** Fill in our data matrices. */
void LabelTest::setUp() {
	int i;
	onecluster.resize(20);
	twocluster.resize(20);
	threecluster.resize(20);
	for (i=0; i<20; i++)
	{
		onecluster[i]= .66+i/200;
		twocluster[i]= 0 + i/200;
		twocluster[i/2] = .66+i/200;
		threecluster[i] = -.66+i/200;
		threecluster[2*i/3] = 0+i/200;
		threecluster[i/3] = .66+i/200;
	}
	oneref.resize(20);
	oneref[0]=0;
	oneref[1]=0;
	oneref[2]=0;
	oneref[3]=0;
	oneref[4]=0;
	oneref[5]=0;
	oneref[6]=0;
	oneref[7]=0;
	oneref[8]=0;
	oneref[9]=0;
	oneref[10]=0;
	oneref[11]=0;
	oneref[12]=0;
	oneref[13]=0;
	oneref[14]=0;
	oneref[15]=0;
	oneref[16]=0;
	oneref[17]=0;
	oneref[18]=0;
	oneref[19]=0;
	tworef.resize(20);
	tworef[0]=0;
	tworef[1]=0;
	tworef[2]=0;
	tworef[3]=0;
	tworef[4]=0;
	tworef[5]=0;
	tworef[6]=0;
	tworef[7]=0;
	tworef[8]=0;
	tworef[9]=0;
	tworef[10]=1;
	tworef[11]=1;
	tworef[12]=1;
	tworef[13]=1;
	tworef[14]=1;
	tworef[15]=1;
	tworef[16]=1;
	tworef[17]=1;
	tworef[18]=1;
	tworef[19]=1;

	threeref.resize(20);
	threeref[0]=0;
	threeref[1]=0;
	threeref[2]=0;
	threeref[3]=0;
	threeref[4]=0;
	threeref[5]=0;
	threeref[6]=0;
	threeref[7]=1;
	threeref[8]=1;
	threeref[9]=1;
	threeref[10]=1;
	threeref[11]=1;
	threeref[12]=1;
	threeref[13]=2;
	threeref[14]=2;
	threeref[15]=2;
	threeref[16]=2;
	threeref[17]=2;
	threeref[18]=2;
	threeref[19]=2;

}

void dump_fun(vector<affx::GType> &testcall, vector<double> &testconf, int length)
{
	int i;
	cout << "{";
	for (i=0; i<length; i++)
		cout << (int)testcall[i] << ",";
	cout << "};\n";
	cout << "{";
	for (i=0; i<length; i++)
		cout << testconf[i] << ",";
	cout << "};\n";
}

void LabelTest::HardTests(){
	// do test from some real data
	vector<affx::GType> testcall;
	vector<double> testconf;
	vector<double> vhardTest;
	vector<double> yvals;
	int i;

	vhardTest.resize(24);
	for (i=0; i<24; i++)
		vhardTest[i]=hardTest[i];
	yvals.resize(24);
	for (i=0; i<24; i++)
		yvals[i] = 0; // dummy values for strength

	SubsetInbred.assign(24,0);

	snp_param sp;
	sp.Initialize();
	testcall.resize(24);
	testconf.resize(24);
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	vector<int> trefcall;
	vector<double> trefconf;
  // encoded AA,AB,BB,NN
	int trefla[]= {1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffa[] = {8.9407e-06,2.08616e-06,1.19209e-07,1.19209e-07,1.19209e-07,0.0542988,7.48634e-05,0,0.00129962,5.96046e-07,0.219659,8.78572e-05,1.07288e-06,0.0223923,0,1.19209e-07,4.76837e-07,1.19209e-07,1.19209e-07,0.0043909,1.19209e-07,0.27395,0.0184641,1.19209e-07};

	for (i=0; i<24; i++)
	{
    // the calls are coded correctly?
    CPPUNIT_ASSERT(affx::GType_valid(trefla[i]));
    CPPUNIT_ASSERT(affx::GType_valid(testcall[i]));
    // and they do match.
		CPPUNIT_ASSERT(testcall[i]==trefla[i]);
    // are the confs close enough?
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffa[i],4));
	}

	// systematically test the parameters of sp to verify
	// that they do what they should.
	// comvar = 0 (1 is default)
	sp.Initialize();
	sp.comvar=0;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflb[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffb[]={7.15256e-07,7.15256e-07,-1.19209e-07,-1.19209e-07,-1.19209e-07,0.000139713,1.07288e-06,5.96046e-07,4.41074e-06,7.15256e-07,0.00105572,0.00340241,7.15256e-07,0.00046581,4.41074e-06,-1.19209e-07,6.7234e-05,-1.19209e-07,2.14577e-06,1.24574e-05,-1.19209e-07,0.000467122,4.43459e-05,-1.19209e-07};

	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflb[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffb[i],4));
	}

	// lambda = 0.5 (1 is default), comvar=1
	sp.Initialize();
	sp.lambda = 0.5;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflc[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffc[]={1.13249e-06,5.96046e-07,0,0,0,0.00361925,5.72205e-06,0,7.50422e-05,3.57628e-07,0.025741,0.000107646,4.76837e-07,0.0117684,0,0,2.38419e-07,0,0,0.000248611,0,0.0129521,0.00106764,0};

	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflc[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffc[i],4));
	}


	// callmethod = 1 (0 is default in code)
	sp.Initialize();
	sp.callmethod=1;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int trefld[] = {1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffd[] = {0.00215191,0.00141519,0,0,0,0.0743099,0.0023703,2.96235e-05,0.0244855,0.000402033,0.132802,7.79629e-05,0.00116622,0.00489968,3.71337e-05,0,5.126e-05,0,0,0.0260765,0,0.128166,0.0721188,0};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],trefld[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffd[i],4));
	}


	// hardshell = 0, 1, 2, 3, SB = 0,0.05, 0.05, 0.5 (2, 0.05 is default)
	// check isotonic regression, hardshell = 3
	sp.Initialize();
	sp.hardshell=3;
	sp.shellbarrier = 0.5;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int trefle[] = {1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffe[] = {2.83122e-05,4.94719e-06,5.96046e-08,5.96046e-08,5.96046e-08,0.0736424,0.000262558,0,0.00395501,1.3113e-06,0.222974,2.92063e-06,2.44379e-06,0.000485122,0,5.96046e-08,0,5.96046e-08,5.96046e-08,0.0110618,5.96046e-08,0.0619845,0.0339516,5.96046e-08};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],trefle[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffe[i],4));
	}


	// bins = 5 (0 is default)
	sp.Initialize();
	sp.bins = 5;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflf[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double trefff[]={8.34465e-07,8.34465e-07,0,0,0,8.34465e-07,8.34465e-07,0,8.34465e-07,8.34465e-07,8.34465e-07,0,8.34465e-07,0.0340099,0,0,0,0,0,8.34465e-07,0,0.0340099,8.34465e-07,0};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflf[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],trefff[i],4));
	}


	// copynumber = 1 (x chromosome) (2 is default)
	sp.Initialize();
	sp.copynumber=1;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflg[] = {2,2,0,0,0,2,2,2,2,2,2,2,2,0,2,0,2,0,0,2,0,0,2,0,};
	double treffg[]={0,0,0,0,0,0,0,0,0,4.41074e-06,0,0,0,0.000879705,0,0,0,0,0,0,0,1.19209e-07,0,0,};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflg[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffg[i],4));
	}

	// hints=1,Hints=<....>, CP = 4, CP=16 (0 is default/0 is default)
	sp.Initialize();
	sp.hints=1;
	Hints.resize(24);
	for(i=0; i<24; i++)
		Hints[i] = treflg[i]; // badly set hints (x chromosome) just to get some action!
	sp.contradictionpenalty=4;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflh[]={2,2,0,0,0,2,2,2,2,2,2,2,2,1,2,0,2,0,0,2,0,1,2,0};
	double treffh[]={0.230311,0.255853,5.96046e-08,5.96046e-08,5.96046e-08,0.0193979,0.194542,0,0.100341,0.36749,0.00874978,6.55651e-07,0.283264,0.020083,0,5.96046e-08,0,5.96046e-08,3.1352e-05,0.0616282,1.78814e-07,0.0933905,0.0327486,5.96046e-08};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflh[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffh[i],4));
	}

	// mix = 0,1,2 (0 is default)
	sp.Initialize();
	sp.mix=1;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int trefli[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffi[]={9.23872e-06,4.29153e-06,-1.19209e-07,-1.19209e-07,-1.19209e-07,0.0266182,4.40478e-05,0,0.000579476,2.08616e-06,0.143353,0.00015825,3.09944e-06,0.0201183,1.78814e-07,-1.19209e-07,1.66893e-06,-1.19209e-07,-1.19209e-07,0.00184906,-1.19209e-07,0.253318,0.00797242,-1.19209e-07};
	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],trefli[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffi[i],4));
	}

	// bic = 2; (0 is default)
	sp.Initialize();
	sp.bic=15; // crank it up to get some action going
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflj[]={1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0};
	double treffj[]={0.00335956,0.00335771,-1.19209e-07,-1.19209e-07,-1.19209e-07,0.0179384,0.00337726,0.271897,0.0037061,0.00335729,0.0623444,0.271874,0.00335747,0.00937039,0.271897,-1.19209e-07,0.271897,-1.19209e-07,-1.19209e-07,0.00453621,-1.19209e-07,0.223942,0.00831538,-1.19209e-07};

	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflj[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffj[i],4));
	}

	// wobble = 5 (tiny is default)
	sp.Initialize();
	sp.wobble = 5;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflk[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffk[]={0.000961304,0.000900745,1.19209e-07,1.19209e-07,1.19209e-07,0.00967038,0.00110471,0,0.00187355,0.000755131,0.0467294,2.39611e-05,0.000856459,0.0047546,0,1.19209e-07,0,1.19209e-07,2.98023e-07,0.00272113,1.19209e-07,0.0689227,0.00490648,1.19209e-07};

	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflk[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffk[i],4));
	}

	// cross-terms xah = 0.3, xhb = 0.3, xhb = -0.3, aa.k=1, ab.k=1, bb.k=1
	sp.Initialize();
	sp.prior.xah=0.3;
	sp.prior.xhb=0.3;
	sp.prior.xab=-0.3;
	sp.prior.aa.k=1;
	sp.prior.ab.k=1;
	sp.prior.bb.k=1;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int trefll[]={1,1,0,0,0,1,1,2,1,1,1,2,1,1,2,0,2,0,0,1,0,1,1,0};
	double treffl[]={0.00184649,0.00149536,1.19209e-07,1.19209e-07,1.19209e-07,0.0497121,0.00280714,0,0.00868309,0.000947833,0.149507,1.10269e-05,0.00128591,0.00152892,0,1.19209e-07,0,1.19209e-07,1.19209e-07,0.0145928,1.19209e-07,0.0519889,0.0279695,1.19209e-07};


	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],trefll[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffl[i],4));
	}

	// prior location
	sp.Initialize();
	sp.prior.aa.m=0;
	sp.prior.bb.m=0;
	bayes_label(testcall,testconf,vhardTest, yvals, Hints,SubsetInbred,sp);
	//dump_fun(testcall, testconf, 24);
	int treflm[]={2,2,1,1,0,2,2,2,2,2,2,2,2,2,2,0,2,1,1,2,1,2,2,1};
	double treffm[]={1.07288e-06,1.07288e-06,1.07288e-06,1.07288e-06,0,1.07288e-06,1.07288e-06,0,1.07288e-06,1.07288e-06,1.07288e-06,1.07288e-06,1.07288e-06,3.51071e-05,1.07288e-06,0,1.07288e-06,1.07288e-06,1.07288e-06,1.07288e-06,1.07288e-06,0.00857729,1.07288e-06,1.07288e-06};

	for (i=0; i<24; i++)
	{
		CPPUNIT_ASSERT(closeEnough (testcall[i],treflm[i],4));
		CPPUNIT_ASSERT(closeEnough (testconf[i],treffm[i],4));
	}

}

void LabelTest::testLabel() {
// test variables
	vector<affx::GType> testcall;
	vector<double> testconf;
	vector<double> yvals;
	int i;

	snp_param sp;
	sp.Initialize();

	testconf.resize(onecluster.size());
	testcall.resize(onecluster.size());
	yvals.resize(onecluster.size());
	SubsetInbred.assign(onecluster.size(),0);
	bayes_label(testcall,testconf,onecluster,yvals,Hints,SubsetInbred,sp);
	for (i=0; i<onecluster.size(); i++)
	{
		CPPUNIT_ASSERT( closeEnough( testcall[i], oneref[i] ) );
		CPPUNIT_ASSERT(closeEnough(testconf[i],0));
	}
	yvals.resize(twocluster.size());
	SubsetInbred.assign(twocluster.size(),0);
	bayes_label(testcall,testconf,twocluster,yvals,Hints,SubsetInbred,sp);
	for (i=0; i<twocluster.size(); i++)
	{
		CPPUNIT_ASSERT(closeEnough(testcall[i],tworef[i]));
		CPPUNIT_ASSERT(closeEnough(testconf[i],0));
	}
	yvals.resize(threecluster.size());
	SubsetInbred.assign(threecluster.size(),0);
	bayes_label(testcall,testconf,threecluster,yvals,Hints,SubsetInbred,sp);
	for (i=0; i<threecluster.size(); i++)
	{
		//cout << testconf[i] << "\t" << threeref[i] << endl;
		CPPUNIT_ASSERT(closeEnough(testcall[i],threeref[i]));
		CPPUNIT_ASSERT(closeEnough(testconf[i],0));
	}

}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LabelTest );

