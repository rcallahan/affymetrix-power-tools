////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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
#include "stats/statfun.h"
//
#include <cstdlib>
#include <fstream>
#include <iostream>
//

using namespace std;

#define TOLERANCE 1e-6  // currently this is the limit of accuracy for pnorm(), it could be made more accurate.
#define MAX_LINE_LEN 1024
#define PNORM_FILE "test_data/pnorm.tab"
#define PWILCOX_FILE "test_data/pwilcox.tab"
#define PSIGNRANK_FILE "test_data/psignrank.tab"
#define SIGNED_RANK_TEST_FILE "test_data/signrank_test.tab"
#define RANK_SUM_TEST_FILE "test_data/ranksum_test.tab"

int test_pnorm(void);
int test_psignrank(void);
int test_pwilcox(void);
int test_signedRankTest(void);
int test_ranksumTest(void);

int main(int argc, char *argv[]) {
  int exit_status = EXIT_SUCCESS;

  if(test_pnorm()==EXIT_FAILURE)
    exit_status = EXIT_FAILURE;

  if(test_psignrank()==EXIT_FAILURE)
    exit_status = EXIT_FAILURE;

  if(test_signedRankTest()==EXIT_FAILURE)
    exit_status = EXIT_FAILURE;

  if(test_pwilcox()==EXIT_FAILURE)
    exit_status = EXIT_FAILURE;

  if(test_ranksumTest()==EXIT_FAILURE)
    exit_status = EXIT_FAILURE;

  exit(exit_status);
}

// Test pnorm
int test_pnorm(void) {
  int exit_status = EXIT_SUCCESS;

  float x,val,ref,error;
  ifstream pnormFile(PNORM_FILE,ios_base::in);
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

  return(exit_status);
}


// Test psignrank
int test_psignrank(void) {
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

  return(exit_status);
}

// Test pwilcox
int test_pwilcox(void) {
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

  return(exit_status);
}

// Test signedRankTest
int test_signedRankTest(void) {
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

  return(exit_status);
}

// Test rank sum test
int test_ranksumTest(void) {
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

  return(exit_status);
}
