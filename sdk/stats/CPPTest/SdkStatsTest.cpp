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
#include "stats/CPPTest/SdkStatsTest.h"
//
#include "stats/statfun.h"
#include "stats/stats-distributions.h"
#include "stats/stats.h"
//
#include "util/AffxConv.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
//


// ----------
CPPUNIT_TEST_SUITE_REGISTRATION( SdkStatsTest );

void SdkStatsTest::setUp()
{
}

void SdkStatsTest::tearDown()
{
}


// ==============================

// currently this is the limit of accuracy for pnorm(), 
// it could be made more accurate.
#define TOLERANCE 1e-6 
//
#define MAX_LINE_LEN 1024
#define PNORM_FILE "../test_data/pnorm.tab"
#define PWILCOX_FILE "../test_data/pwilcox.tab"
#define PSIGNRANK_FILE "../test_data/psignrank.tab"
#define SIGNED_RANK_TEST_FILE "../test_data/signrank_test.tab"
#define RANK_SUM_TEST_FILE "../test_data/ranksum_test.tab"

int test_pnorm(void);
int test_psignrank(void);
int test_pwilcox(void);
int test_signedRankTest(void);
int test_ranksumTest(void);

// Test pnorm
void SdkStatsTest::test_pnorm(void) {
  int exit_status = EXIT_SUCCESS;

  float x,val,ref,error;
  ifstream pnormFile(PNORM_FILE);
  int lineNo=0;
  float max_error=0;
  char line[MAX_LINE_LEN];
  while(pnormFile.getline(line,MAX_LINE_LEN)) {
    lineNo++;
    if(2 != sscanf(line,"%f\t%f",&x,&ref)) {
      fprintf(stderr,"bad input in line %d of %s\n",lineNo,PNORM_FILE);
      exit(EXIT_FAILURE);
    }
    val = affxstat::pnorm(x,0,1,1,0);
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for pnorm(%f), got %.10f, should be %.10f\n",fabs(ref-val),x,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;
  }
  pnormFile.close();
  printf("Maximum observed error for pnorm() is %.10f\n",max_error);

  CPPUNIT_ASSERT(exit_status==0);
  //return(exit_status);
}


// Test psignrank
void SdkStatsTest::test_psignrank(void) {
  int exit_status = EXIT_SUCCESS;

  float val,ref,error;
  unsigned int n,t;
  ifstream psignrankFile(PSIGNRANK_FILE);
  int lineNo=0;
  float max_error=0;
  char line[MAX_LINE_LEN];
  while(psignrankFile.getline(line,MAX_LINE_LEN)) {
    lineNo++;
    if(3 != sscanf(line,"%u\t%u\t%f",&n,&t,&ref)) {
      fprintf(stderr,"bad input in line %d of %s\n",lineNo,PSIGNRANK_FILE);
      exit(EXIT_FAILURE);
    }
    val = affxstat::psignrank(t,n,true,false);
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for psignrank(%u,%u), got %.10f, should be %.10f\n",error,t,n,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;
  }
  psignrankFile.close();
  printf("Maximum observed error for psignrank() is %.10f\n",max_error);

  CPPUNIT_ASSERT(exit_status==0);
  //return(exit_status);
}

// Test pwilcox
void SdkStatsTest::test_pwilcox(void) {
  int exit_status = EXIT_SUCCESS;

  float val,ref,error;
  unsigned int m,n,x;
  ifstream pwilcoxFile(PWILCOX_FILE);
  int lineNo=0;
  float max_error=0;
  char line[MAX_LINE_LEN];
  while(pwilcoxFile.getline(line,MAX_LINE_LEN)) {
    lineNo++;
    if(4 != sscanf(line,"%u\t%u\t%u\t%f",&m,&n,&x,&ref)) {
      fprintf(stderr,"bad input in line %d of %s\n",lineNo,PWILCOX_FILE);
      exit(EXIT_FAILURE);
    }
    val = affxstat::pwilcox(x,m,n,true,false);
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for pwilcox(%u,%u,%u), got %.10f, should be %.10f\n",fabs(ref-val),x,m,n,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;
  }
  pwilcoxFile.close();
  printf("Maximum observed error for pwilcox() is %.10f\n",max_error);

  CPPUNIT_ASSERT(exit_status==0);
  //return(exit_status);
}

// Test signedRankTest
void SdkStatsTest::test_signedRankTest(void) {
  int exit_status = EXIT_SUCCESS;

  float val,ref,error;
  unsigned int n;
  ifstream inFile(SIGNED_RANK_TEST_FILE);
  int lineNo=0;
  int tailType;
  float max_error=0;
  while(!inFile.eof()) {
    lineNo++;
    inFile >> ref >> tailType >> n;
    if(inFile.fail())
      break;
    double *x = new double[n];
    for(unsigned int i=0; i<n; i++)
      inFile >> x[i];

    // Check p-value on linear scale
    val = affxstat::signedRankTest(x,n,tailType,false);
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for signedRankTest on line %d of %s, got %.10f, should be %.10f\n",error,lineNo,SIGNED_RANK_TEST_FILE,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;

    // Check p-value on log scale
    val = exp(affxstat::signedRankTest(x,n,tailType,true));
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for signedRankTest on line %d of %s (logged pvalue), got %.10f, should be %.10f\n",error,lineNo,SIGNED_RANK_TEST_FILE,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;

    delete[] x;
  }
  inFile.close();
  printf("Maximum observed error for signedRankTest() is %.10f\n",max_error);

  CPPUNIT_ASSERT(exit_status==0);
  //return(exit_status);
}

// Test rank sum test
void SdkStatsTest::test_ranksumTest(void) {
  int exit_status = EXIT_SUCCESS;

  float val,ref,error;
  unsigned int n,m;
  ifstream inFile(RANK_SUM_TEST_FILE);
  int lineNo=0;
  int tailType;
  float max_error=0;
  while(!inFile.eof()) {
    lineNo++;
    inFile >> ref >> tailType >> n >> m;
    if(inFile.fail())
      break;
    double *x = new double[n];
    double *y = new double[m];
    for(unsigned int i=0; i<n; i++)
      inFile >> x[i];
    for(unsigned int i=0; i<m; i++)
      inFile >> y[i];

    // test on linear pvalue
    val = affxstat::ranksumTest(x,n,y,m,tailType,false);
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for ranksumTest on line %d of %s, got %.10f, should be %.10f\n",error,lineNo,RANK_SUM_TEST_FILE,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;

    // test on logged pvalue
    val = exp(affxstat::ranksumTest(x,n,y,m,tailType,true));
    error = fabs(val-ref);
    if(error > TOLERANCE) {
      printf("error of size %.10f for ranksumTest on line %d of %s (logged version), got %.10f, should be %.10f\n",error,lineNo,RANK_SUM_TEST_FILE,val,ref);
      exit_status = EXIT_FAILURE;
    }
    if(error > max_error)
      max_error = error;

    delete[] x;
    delete[] y;
  }
  inFile.close();
  printf("Maximum observed error for ranksumTest() is %.10f\n",max_error);

  CPPUNIT_ASSERT(exit_status==0);
  //return(exit_status);
}


//void SdkStatsTest::test_SdkStats()
//{
//  test_pnorm
//
//  if(test_pnorm()==EXIT_FAILURE)
//    exit_status = EXIT_FAILURE;
//
//  if(test_psignrank()==EXIT_FAILURE)
//    exit_status = EXIT_FAILURE;
//
//  if(test_signedRankTest()==EXIT_FAILURE)
//    exit_status = EXIT_FAILURE;
//
//  if(test_pwilcox()==EXIT_FAILURE)
//    exit_status = EXIT_FAILURE;
//
//  if(test_ranksumTest()==EXIT_FAILURE)
//    exit_status = EXIT_FAILURE;
//
//  exit(exit_status);
//}
//
//}

// ====================

// rounded to 5 places
float test_uprob_data[]={
  0.50000f,
  0.46017f,
  0.42074f,
  0.38209f,
  0.34458f,
  0.30854f,
  0.27425f,
  0.24196f,
  0.21186f,
  0.18406f };

void 
SdkStatsTest::test_uprob()
{
  //printf("test_uprob:\n");
  for (int n=0;n<10;n++) {
    float p=affxstat::uprob(n/10.0);
    //printf("uprob(%.2f)=%.5f\n",n/10.0,p);
    CPPUNIT_ASSERT(fabs(p-test_uprob_data[n])<0.00001);
  }
}

float test_chisqrprob_data[]={
  1.00000f,
  0.95123f,
  0.90484f,
  0.86071f,
  0.81873f,
  0.77880f,
  0.74082f,
  0.70469f,
  0.67032f,
  0.63763f };

void 
SdkStatsTest::test_chisqrprob()
{
  //printf("test_chisqrprob:\n");
  for (int n=2;n<3;n++) {
    for (int x=0;x<10;x++) {
      float p=affxstat::chisqrprob(n,x/10.0);
      //printf("chisqrprob(%d,%.2f)=%.5f\n",n,x/10.0,p);
      CPPUNIT_ASSERT(fabs(p-test_chisqrprob_data[x])<0.00001);
    }
  }
}

float log10transform(float x)
{
	return log10(x);
}

void
SdkStatsTest::test_PearsonCorrelation()
{
	float x1[] = {
		0.552628642f,
		0.791982055f,
		0.730823791f,
		0.790206725f,
		0.962450187f,
		0.825606768f,
		0.180500524f,
		0.56035321f,
		0.624721932f,
		0.590874686f,
		0.858872657f,
		0.004050031f,
		0.155189576f,
		0.743305177f,
		0.676323844f,
		0.881118343f,
		0.29733674f,
		0.37398256f,
		0.43675202f,
		0.762386693f
	};
	float x2[] = {
		0.725074466f,
		0.030247519f,
		0.899631008f,
		0.02357856f,
		0.717061428f,
		0.492287055f,
		0.060332963f,
		0.672098711f,
		0.557733722f,
		0.605800393f,
		0.615076848f,
		0.627646562f,
		0.029739631f,
		0.527833176f,
		0.841452178f,
		0.812548337f,
		0.72512321f,
		0.826932845f,
		0.189030606f,
		0.396031268f
	};
	int size = sizeof(x1) / sizeof(float);
	float expected_result = 0.193280778f;
	const float eps = 1e-6f;
	double result = affxstat::PearsonCorrelation(x1, x2, size);		
	CPPUNIT_ASSERT_DOUBLES_EQUAL(result, expected_result, eps);

	// Now with LOG10 transform.
	expected_result = 0.051642286f;
	result = affxstat::PearsonCorrelation(x1, x2, size, log10transform);		
	CPPUNIT_ASSERT_DOUBLES_EQUAL(result, expected_result, eps);

	// Now with no data.
	expected_result = NAN;
	result = affxstat::PearsonCorrelation(NULL, NULL, 0);		
//  CPPUNIT_ASSERT(isnan(result));	
    CPPUNIT_ASSERT(result != result); // This is an alternative test for NaN.


}
void SdkStatsTest::test_CalcHWEqPValue()
{
    double expected[] = {
        1.0,
        0.63735188823393707,
        1.0000000000000000,
        0.72903448953880390,
        0.15729920705028488,
        0.56370286165077310,
        0.083264516663550503,
        0.35064788970445893
    };
    int index = 0;
    for (int ia=0; ia<2; ia++)
    {
        for (int ib=1; ib<3; ib++)
        {
            for (int iab=0; iab<2; iab++)
            {
                double hw = affxstat::CalcHWEqPValue(ia, iab, ib);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(hw, expected[index++], 0.000001);
            }
        }
    }
}
