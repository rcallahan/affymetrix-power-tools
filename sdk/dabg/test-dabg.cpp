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

/// @file   test-dabg.cpp
/// @author harley
/// @brief  Collection of tests of the AffyDabg SDK

//
#include "dabg/Dabg.h"
//
#include "stats/stats-distributions.h"
//
#include <cassert>
#include <cmath>

// =======

#define TEST_DATA_FILENAME "check-data-1.txt"

// ==========

void
dabg_uniform_fill(affx::Dabg& dabg)
{
  // fill with a uniform set of values.
  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    for (int i=0;i<100;i++) {
      dabg.add_intensity(gcbin,i);
    }
    assert(dabg.intensity_count(gcbin)==100);
  }
}


/// @brief          Test adding values to a dabg object.
void
run_test_1()
{
  affx::Dabg dabg;
  int gcbin=0;
  //DABG_T intensity;
  double pval;

  //
  dabg.intensity_reserve(2000);

  //
  pval=dabg.intensity_to_pvalue( 0,  0);
  pval=dabg.intensity_to_pvalue( 0,100);
  pval=dabg.intensity_to_pvalue(25,  0);
  pval=dabg.intensity_to_pvalue(25,100);

  //
  for (int i=0;i<=AFFYDABG_GC_COUNT_MAX;i++) { // <= for 0-25
    dabg.add_intensity(i,(affx::Dabg::intensity_t)i);
  }

  // These two datafile should be the same
  dabg.dabgfile_write("dabg-1.data");
  dabg.dabgfile_read("dabg-1.data");
  dabg.dabgfile_write("dabg-2.data");

  //
#define I_MAX   99
#define I_PNTS  20
  for (int i=1;i<I_MAX;i++) {
    dabg.add_intensity(gcbin,(affx::Dabg::intensity_t)(I_MAX-i));
  }
  // dump the data before and after the prepare.
  // dabg.dump();
  // dabg.prepare();
  // dabg.dump(3);

  for (float iv=-(I_MAX/I_PNTS);iv<(I_MAX+(I_MAX/I_PNTS));iv+=(int)(I_MAX/I_PNTS)) {
    pval=dabg.intensity_to_pvalue(gcbin,iv);
  }
}


// ==========

void
print_dabg_results(int probe_cnt,double* target_out,double* pvalue_out)
{
  printf("target=%8.6f ",*target_out);
  for (int p=0;p<probe_cnt;p++) {
    printf("pval[%2d]=%6.8f ",p,pvalue_out[p]);
  }
  printf("\n");
}



// how many probes in this test
#define T2_PROBE_CNT   5
// the PM value
#define T2_PM       10.0

/// @brief          Another test of AffyDabg
void
run_test_2()
{
  affx::Dabg dabg;
  //dabg_uniform_fill(dabg);

  int probe_cnt=T2_PROBE_CNT; // sizeof(gc_arr)/sizeof(gc_arr[0]);
  int gc_arr[T2_PROBE_CNT]={10,11,12,13,14};
  double pm_arr[T2_PROBE_CNT]={T2_PM,T2_PM,T2_PM,T2_PM,T2_PM};
  //
  double target_out=0.0;
  double pvalue_out[T2_PROBE_CNT]={0.0,0.0,0.0,0.0,0.0};

  // test a single chip
  dabg.compute_target_fisher(probe_cnt,
                             (int*)gc_arr,
                             (double*)pm_arr, 0,
                             (double*)&target_out,
                             (double*)pvalue_out);

  print_dabg_results(T2_PROBE_CNT,&target_out,pvalue_out);
}


// ==========

// Data from a test run of dabg_classic
//
// DB<3> x @atomID
// 0  2598546
// 1  2598547
// 2  2598548
// 3  2598549
// # x @int
//   DB<4> x @int
// 0  66
// 1  86
// 2  84
// 3  142
//    # x @pvalues
//  DB<5> x @pvalues
//0  0.750762970498474
//1  0.750251762336354
//2  0.766364551863041
//3  0.696266397578204
//    # x @dist

// LOWER: target=0.967344 pval[ 0]=0.75584942 pval[ 1]=0.75125879 pval[ 2]=0.77139980 pval[ 3]=0.69657254
// PERL:                           0.75076297          0.75025176          0.76636455          0.69626639
// AVG:   target=0.964861 pval[ 0]=0.74771106 pval[ 1]=0.74723065 pval[ 2]=0.76636457 pval[ 3]=0.69254029
// UPPER: target=0.962130 pval[ 0]=0.73957276 pval[ 1]=0.74320245 pval[ 2]=0.76032227 pval[ 3]=0.68850803

//    # x scalar($iholder->tIntensities())
//    # x $testStatistic

// DB<19> x $self->{'target'}{$celname}
// 0  0.96605

#define T3_PROBE_CNT 4

/// @brief          Yet another DABG test
void
run_test_3(const std::string& fname)
{
  affx::Dabg dabg;
  int probe_cnt=T3_PROBE_CNT; // sizeof(gc_arr)/sizeof(gc_arr[0]);
  int gc_arr[T3_PROBE_CNT]={0,1,1,3};
  double pm_arr[T3_PROBE_CNT]={66,86,84,142};
  //
  double target_out=0.0;
  double pvalue_out[T3_PROBE_CNT]={0.0,0.0,0.0,0.0};

  // fill with a uniform set of values.
  dabg.dabgfile_read(fname);
  // dabg.dump(1); //

  // test a single chip
  dabg.compute_target_fisher(probe_cnt,
                             (int*)gc_arr,
                             (double*)pm_arr, 0,
                             (double*)&target_out,
                             (double*)pvalue_out);

  print_dabg_results(T3_PROBE_CNT,&target_out,pvalue_out);
}

// ==========
// test the chisqrprob function with the perl generated
// pvals. our pvals are off which makes the chisqrprob off
// when using our data.

// the perl output
// classic: ii= 0 pm=66  dist=983 pval=  0.75076297 pval_product=  0.75076297
// classic: ii= 1 pm=86       993 pval=  0.75025176 pval_product=  0.56326124
// classic: ii= 2 pm=84       993 pval=  0.76636455 pval_product=  0.43166345
// classic: ii= 3 pm=142      991 pval=  0.69626640 pval_product=  0.30055275
// classic: ========== target:   0.96605000

// c++ output
// test4: i= 0 pval=  0.75076297 pval_product=  0.75076297
// test4: i= 1 pval=  0.75025176 pval_product=  0.56326124
// test4: i= 2 pval=  0.76636455 pval_product=  0.43166345
// test4: i= 3 pval=  0.69626640 pval_product=  0.30055275
// test4: target=  0.96604574


/// @brief          Test of affxstat::chisqrprob



void
test_target(int pval_cnt,double* pval_arr,double target_ref)
{
  double pval_product=1.0;
  double pval;
  double target;

  for (int i=0;i<pval_cnt;i++) {
    pval=pval_arr[i];
    pval_product=pval_product*pval;
    printf("target: %02d: pval=%12.8f pval_product=%12.8f\n",i,pval,pval_product);
  }

  //
  double test_statistic=-2*log(pval_product);
  target=affxstat::chisqrprob(2*pval_cnt,test_statistic);
  printf("target: test_statistic=%12.8f\n",test_statistic);
  printf("target: target=%12.8f (perl=%12.8f) (diff=%12.8f)\n",target,target_ref,(target-target_ref));
  //
  // assert(fabs(target-target_ref)<0.001);
}

void
run_test_4()
{
  // these values came from the classic
  double pval_arr[]={
    0.750762970498474,
    0.750251762336354,
    0.766364551863041,
    0.696266397578204};
  int pval_cnt=sizeof(pval_arr)/sizeof(pval_arr[0]);

  test_target(pval_cnt,pval_arr,0.96605);
}


void
run_test_5()
{
  double pval_arr[]={0.464466393640087 ,
                     0.0921228304405874,
                     0.122945430637738 ,
                     0.4543063773833   ,
                     0.415962837837838 ,
                     0.52548435171386  ,
                     0.256191102344948 ,
                     0.200674536256324 ,
                     0.106044252563411 ,
                     0.834927140255009 ,
                     0.279932546374368 };
  int pval_cnt=sizeof(pval_arr)/sizeof(pval_arr[0]);

  // teststat= 28.4450193272307
  test_target(pval_cnt,pval_arr,0.16140);
}

//////////

#define T6_PROBE_CNT  6
#define T6_PM        50
#define HBAR(x) { printf("%-10s==============================\n",x); }

void
run_test_6()
{
  affx::Dabg dabg;

  // make a reference
  dabg_uniform_fill(dabg);

  int probe_cnt=T6_PROBE_CNT;
  int gc_arr[T6_PROBE_CNT]   ={10   ,   11,   12,   13,   14,   15};
  double pm_arr[T6_PROBE_CNT]={T6_PM,T6_PM,T6_PM,T6_PM,T6_PM,T6_PM};
  double out_target=0.0;
  double out_pvalue[T6_PROBE_CNT]={0.0,0.0,0.0,0.0,0.0,0.0};

  HBAR("1");

  // flat
  dabg.compute_target_percentile(probe_cnt,
                                 (int*)gc_arr,
                                 (double*)pm_arr,
                                 0.5, // per_cut
                                 (double*)&out_target,
                                 (double*)out_pvalue);

  // test a single point
  HBAR("2");
#define T6_AMAX (10.0)
  for (int i=0;i<=T6_AMAX;i++) {
    dabg.compute_target_percentile(1,
                                   (int*)gc_arr,
                                   (double*)pm_arr,
                                   (i/T6_AMAX), // per_cut
                                   (double*)&out_target,
                                   (double*)out_pvalue);
    assert(out_target==0.5);
  }

  // test interpolation between two points over the range 0.01 to 1.00
  // It should be linear
  HBAR("3");
#define T6_BMAX (10.0)
  pm_arr[0]=  0.0;
  pm_arr[1]=100.0;
  // calc by hand
  double out_b[]={0.010,0.109,0.208,0.307,0.406,0.505,0.604,0.703,0.802,0.901,1.000};
  for (int i=0;i<=T6_BMAX;i++) {
    dabg.compute_target_percentile(2,
                                   (int*)gc_arr,
                                   (double*)pm_arr,
                                   (i/T6_BMAX), // per_cut
                                   (double*)&out_target,
                                   (double*)out_pvalue);
    // double diff=fabs(out_target-out_b[i]);
    // if (diff<
    // assert(<0.0001);
    assert(fabs(out_target-out_b[i])<0.00001);
  }

  // Test that we get a step as the probes are marked down.
  // As the cross the middle it will interpolate.
  HBAR("4"); //                 V middle value
  double out_c[]={0.50,0.50,0.50,0.75,1.00,1.00,1.00};
  for (int i=0;i<T6_PROBE_CNT;i++) {
    pm_arr[i]=50.0;
  }
  for (int i=0;i<T6_PROBE_CNT+1;i++) {
    dabg.compute_target_percentile(T6_PROBE_CNT,
                                   (int*)gc_arr,
                                   (double*)pm_arr,
                                   0.50, // per_cut at the middle
                                   (double*)&out_target,
                                   (double*)out_pvalue);
    assert(out_target==out_c[i]);
    if (i<T6_PROBE_CNT) {
      pm_arr[i]=00.0;
    }
  }
}

//////////

/// @brief          Entry point for running the test functions.
int
main(int argc,const char* argv[])
{
  //convert_file();
  run_test_1();
  run_test_2();

  if (argc<2) {
    //printf("### %s: skipping test3 -- no input file given\n",argv[0]);
    run_test_3(TEST_DATA_FILENAME);
  } else {
    run_test_3(argv[1]);
  }

  // test combining pvalues to targets
  run_test_4();
  run_test_5();

  run_test_6();

  return 0;
}
