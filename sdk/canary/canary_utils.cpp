////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#include "canary/canary_utils.h"
//
#include "canary/CanaryModel.h"
#include "canary/CanaryPrior.h"
//
#include "file/TsvFile/TsvFile.h"
#include "stats/stats-distributions.h"
//
#include <cfloat> // for FLT_MIN
//
using namespace affx;


pair<list<string>,map<string,valarray<double> > >
load_broad_intensity_map(string ifname)
  {

  // List of column names which are sample names
  list<string>sample_name;

  // Associate cnv region names with intensity vectors
  map<string,valarray<double> >intensity_map;

  // Open the ascii tab delimited input file
  TsvFile infile;
  if(infile.openTable(ifname) != TSV_OK)
      Err::errAbort("Unable to open intensity map file '" + ifname + "'");


  // Retrieve the set of column names - no error checking for now
  infile.nextLevel(0);
  for (int k=1; k<infile.getColumnCount(infile.lineLevel()); k++)
    {
    string column_str;
    infile.get(infile.lineLevel(),k,column_str);
    sample_name.push_back(column_str);
    }

  // Error checking will be easier if I know the number of columns
  unsigned long ncol=sample_name.size();

  
  while(infile.nextLevel(0)==TSV_OK)
    {
    string cnv_region;
    double intensity;

    // load the set of cnv region names
    infile.get(infile.lineLevel(),0,cnv_region);

    // a place to put the intensities for a cnv region
    vector<double>intensity_vec;

    // go fetch them.
    for (int k=1; k<=ncol; k++)
      {
      infile.get(infile.lineLevel(),k,intensity);
      intensity_vec.push_back(intensity);
      }

    // store the intensities
    valarray<double>vals(intensity_vec.size());
    for (int k=0; k<intensity_vec.size(); k++) vals[k] = intensity_vec[k];
    intensity_map[cnv_region].resize(vals.size());
    intensity_map[cnv_region] = vals;
    }

  infile.close();

  return pair<list<string>,map<string,valarray<double> > >
         (sample_name,intensity_map);
  }


map<string,CanaryPrior>
load_broad_prior_map(string ifname)
  {
  map<string,CanaryPrior>prior_map;

  // Open the ascii tab delimited input file
//cnv_region	P0	P1	P2	P3	P4	M0	M1	M2	M3	M4	V0	V1	V2	V3	V4
//CNP10	0.0076	0.0608	0.9316	0	0	0.3162	0.6032	0.8655	1.1052	1.3243	0.0015	0.001	0.0015	0.0032	0.0047
  TsvFile infile;

  string cnv_region;
  double p0,p1,p2,p3,p4,m0,m1,m2,m3,m4,v0,v1,v2,v3,v4;
  infile.bind(0,"cnv_region",&cnv_region, TSV_BIND_REQUIRED);
  infile.bind(0,"P0",&p0, TSV_BIND_REQUIRED);
  infile.bind(0,"P1",&p1, TSV_BIND_REQUIRED);
  infile.bind(0,"P2",&p2, TSV_BIND_REQUIRED);
  infile.bind(0,"P3",&p3, TSV_BIND_REQUIRED);
  infile.bind(0,"P4",&p4, TSV_BIND_REQUIRED);
  infile.bind(0,"M0",&m0, TSV_BIND_REQUIRED);
  infile.bind(0,"M1",&m1, TSV_BIND_REQUIRED);
  infile.bind(0,"M2",&m2, TSV_BIND_REQUIRED);
  infile.bind(0,"M3",&m3, TSV_BIND_REQUIRED);
  infile.bind(0,"M4",&m4, TSV_BIND_REQUIRED);
  infile.bind(0,"V0",&v0, TSV_BIND_REQUIRED);
  infile.bind(0,"V1",&v1, TSV_BIND_REQUIRED);
  infile.bind(0,"V2",&v2, TSV_BIND_REQUIRED);
  infile.bind(0,"V3",&v3, TSV_BIND_REQUIRED);
  infile.bind(0,"V4",&v4, TSV_BIND_REQUIRED);

  if(infile.open(ifname) != TSV_OK)
      Err::errAbort("Couldn't open file: " + ifname + " to read.");

  // Just go through and pick up the region, probabilities, means and variances.
  while(infile.nextLevel(0)==TSV_OK)
    {
    valarray<double>prop(5);
    valarray<double>mean(5);
    valarray<double>var(5);

    prop[0]=p0;
    prop[1]=p1;
    prop[2]=p2;
    prop[3]=p3;
    prop[4]=p4;
    mean[0]=m0;
    mean[1]=m1;
    mean[2]=m2;
    mean[3]=m3;
    mean[4]=m4;
    var[0]=v0;
    var[1]=v1;
    var[2]=v2;
    var[3]=v3;
    var[4]=v4;

    prior_map[cnv_region] = CanaryPrior(cnv_region,prop,mean,var);
    }

  return prior_map;
  }


double hardy_weinberg(vector<int>cluster_vec,valarray<double>props)
  {
  CanaryModel CM;
  if (cluster_vec == CM["model0"]) return 1.0;

  if (cluster_vec == CM["model2"]) return 1.0;

  if (cluster_vec == CM["model01"])
    return hardy_weinberg2(props[0],props[1],0.0);

  if (cluster_vec == CM["model12"])
    return hardy_weinberg2(0.0,props[0],props[1]);

  if (cluster_vec == CM["model23"])
    return hardy_weinberg2(props[0],props[1],0.0);

  if (cluster_vec == CM["model34"])
    return hardy_weinberg2(0.0,props[0],props[1]);

  if (cluster_vec == CM["model012"])
    return hardy_weinberg2(props[0],props[1],props[2]);

  if (cluster_vec == CM["model234"])
    return hardy_weinberg2(props[0],props[1],props[2]);

  if (cluster_vec == CM["model123"])
//	return hardy_weinberg3(0.0,props[0],props[1],props[2],0.0);
	return hardy_weinberg3_bad(0.0,props[0],props[1],props[2],0.0);

  if (cluster_vec == CM["model0123"])
//	return hardy_weinberg3(props[0],props[1],props[2],props[3],0.0);
	return hardy_weinberg3_bad(props[0],props[1],props[2],props[3],0.0);

  if (cluster_vec == CM["model1234"])
//	return hardy_weinberg3(0.0,props[0],props[1],props[2],props[3]);
	return hardy_weinberg3_bad(0.0,props[0],props[1],props[2],props[3]);

  if (cluster_vec == CM["model01234"])
//	return hardy_weinberg3(props[0],props[1],props[2],props[3],props[4]);
	return hardy_weinberg3_bad(props[0],props[1],props[2],props[3],props[4]);

  else 
      Err::errAbort("!!!found no model!!!\n");

  return 1.0;
  }


double hardy_weinberg2(double Naa, double Nab, double Nbb)
  {
  // Nxx is the proportion multiplied by the count
  // Exx is the expected hw value 
  // prob_a is the probability of an a allele.
  double count = (Naa + Nab + Nbb);

  double prob_a = (2*Naa + Nab)/(2*count);

  double Eaa = count*prob_a*prob_a;
  double Eab = count*2.0*prob_a*(1.0 - prob_a);
  double Ebb = count*(1.0 - prob_a)*(1.0 - prob_a);

  double statistic = 0.0;
  if (Eaa > FLT_MIN) statistic += (Eaa-Naa)*(Eaa-Naa)/Eaa;
  if (Eab > FLT_MIN) statistic += (Eab-Nab)*(Eab-Nab)/Eab;
  if (Ebb > FLT_MIN) statistic += (Ebb-Nbb)*(Ebb-Nbb)/Ebb;

  return 1.0 - affxstat::chisqrprob(1,statistic);
  }

// N0 - N4 are the counts of observed copy numbers
double hardy_weinberg3(double N0,double N1,double N2,double N3,double N4)
  {
  double count = N0 + N1 + N2 + N3 + N4;

  double f0 = N0/count; // observed proportion with copy number of 0
  double f1 = N1/count; // and so on for f1 - f4
  double f2 = N2/count;
  double f3 = N3/count;
  double f4 = N4/count;

  double p0 = 1.0/3.0; // starting point for iterations
  double p1 = 1.0/3.0;	
  double p2 = 1.0/3.0;

  double was_p0 = p0; // old values for checking convergence
  double was_p1 = p1;
  double was_p2 = p2;

  double converge = 0.0000001; // did it converge

  // There probably is a closed form solution for this but instead let
  // the computer do the work.
  //
  // The observed counts with a copy number of 2 can come either from
  // a pair of strands with single copies or one strand having a copy
  // number of 0 and the other a copy number of 2.  This deconvolution
  // causes a loss of 3 degrees of freedom because observed values are
  // being imputed from expected values which would not happen with 
  // three alleles producing distinct phenotypes for each combination.
  //
  // The plan is to start with initial values and substitute the unobserved
  // proportions with imputed proportions for copy number of 2.  This can 
  // probably be written out as an expectation step in an EM algorithm if
  // fitting the H-W model to the data.
  
  double diff;
  do
    {
    p0 = (2.0*f0 + f1 + f2*2.0*p0*p2/(2.0*p0*p2 + p1*p1))/2.0;
    p1 = (f2*2.0*p1*p1/(2.0*p0*p2 + p1*p1) + f1 + f3)/2.0;
    p2 = (2*f4 + f2*2.0*p0*p2/(2.0*p0*p2 + p1*p1) + f3)/2.0;

	double diff_max = 0.0;

	diff = abs(p0 - was_p0);
    was_p0 = p0;
    if (diff > diff_max) diff_max = diff;

	diff = abs(p1 - was_p1);
    was_p1 = p1;
    if (diff > diff_max) diff_max = diff;

	diff = abs(p2 - was_p2);
    was_p2 = p2;
    if (diff > diff_max) diff_max = diff;

    } while (diff > converge);

  double statistic = 0.0;
  double E0 = count*p0*p0;
  statistic += (E0-N0)*(E0-N0)/E0;

  double E1 = count*2.0*p0*p1;
  statistic += (E1-N1)*(E1-N1)/E1;

  double E2 = count*(2.0*p0*p2 +  p1*p1);
  statistic += (E2-N2)*(E2-N2)/E2;

  double E3 = count*2.0*p1*p2;
  statistic += (E3-N3)*(E3-N3)/E3;

  double E4 = count*p2*p2;
  statistic += (E4-N4)*(E4-N4)/E4;

  return 1.0 - affxstat::chisqrprob(2,statistic);
  }


// Use this defective version to repeat broad results.
double hardy_weinberg3_bad(double N0,double N1,double N2,double N3,double N4)
  {
  double count = N0 + N1 + N2 + N3 + N4;

  double f0 = N0/count; // observed proportion with copy number of 0
  double f1 = N1/count; // and so on for f1 - f4
//  double f2 = N2/count;
  double f3 = N3/count;
  double f4 = N4/count;

  double freq0 = (2*f0 + f1)/2;
  double freq2 = (2*f4 + f3)/2;
  double freq1 = 1.0 - freq0 - freq2;

  double statistic=0.0;


  double E0 = count*freq0*freq0;
  if (E0 > FLT_MIN) statistic += (E0-N0)*(E0-N0)/E0;

  double E1 = count*2.0*freq0*freq1;
  if (E1 > FLT_MIN) statistic += (E1-N1)*(E1-N1)/E1;

  double E2 = count*(freq1*freq1 + 2*freq0*freq2);
  if (E2 > FLT_MIN) statistic += (E2-N2)*(E2-N2)/E2;

  double E3 = count*2.0*freq1*freq2;
  if (E3 > FLT_MIN) statistic += (E3-N3)*(E3-N3)/E3;

  double E4 = count*freq2*freq2;
  if (E4 > FLT_MIN) statistic += (E4-N4)*(E4-N4)/E4;

  return 1.0 - affxstat::chisqrprob(2,statistic);
  }


// Repeats the final fix-it-up step in the Broad Implementation.
BroadEstepper1 directional_impute(CanaryOptions& opts,
                                  pair<string,BroadEstepper1> model_pair,
                                  CanaryPrior CP)
  {
  valarray<double> freqs = CP.prop();
  string model = model_pair.first;
  BroadEstepper1 BE = model_pair.second;

  if ((freqs[3]+freqs[4])==0.0)
    {
    if (model=="model012") return BE;

    vector<int>cvec_forced;
    cvec_forced.push_back(0);
    cvec_forced.push_back(1);
    cvec_forced.push_back(2);

    if (model=="model0")
      {
      double delta = BE.mean(0) - CP.mean(0);
      valarray<double>mean(3);
      mean[0] = BE.mean(0);
      mean[1] = CP.mean(1) + delta;
      mean[2] = CP.mean(2) + delta;
      
      valarray<double>var(3);
      var[0] = BE.var(0);
      var[1] = CP.var(1);
      var[2] = CP.var(2);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }

    if (model=="model2")
      {
      double delta = BE.mean(0) - CP.mean(2);
      valarray<double>mean(3);
      mean[0] = CP.mean(0) + delta;
      mean[1] = CP.mean(1) + delta;
      mean[2] = BE.mean(0);
      
      valarray<double>var(3);
      var[0] = CP.var(0);
      var[1] = CP.var(1);
      var[2] = BE.var(0);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }


    if (model=="model01")
      {
      double delta = BE.mean(0) - CP.mean(0);
      delta += BE.mean(1) - CP.mean(1);
      delta /= 2;

      valarray<double>mean(3);
      mean[0] = BE.mean(0);
      mean[1] = BE.mean(1);
      mean[2] = CP.mean(2) + delta;
      
      valarray<double>var(3);
      var[0] = BE.var(0);
      var[1] = BE.var(1);
      var[2] = CP.var(2);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }

    if (model=="model12")
      {
      double delta = BE.mean(0) - CP.mean(1);
      delta += BE.mean(1) - CP.mean(2);
      delta /= 2;

      valarray<double>mean(3);
      mean[0] = CP.mean(0) + delta;
      mean[1] = BE.mean(0);
      mean[2] = BE.mean(1);
      
      valarray<double>var(3);
      var[0] = CP.var(0);
      var[1] = BE.var(0);
      var[2] = BE.var(1);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }
    }
  

  if ((freqs[0]+freqs[1])==0.0)
    {
    if (model=="model234") return BE;

    vector<int>cvec_forced;
    cvec_forced.push_back(2);
    cvec_forced.push_back(3);
    cvec_forced.push_back(4);

    if (model=="model2")
      {
      double delta = BE.mean(0) - CP.mean(2);
      valarray<double>mean(3);
      mean[0] = BE.mean(0);
      mean[1] = CP.mean(3) + delta;
      mean[2] = CP.mean(4) + delta;
      
      valarray<double>var(3);
      var[0] = BE.var(0);
      var[1] = CP.var(3);
      var[2] = CP.var(4);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }

    if (model=="model23")
      {
      double delta = BE.mean(0) - CP.mean(2);
      delta += BE.mean(1) - CP.mean(3);
      delta /= 2;

      valarray<double>mean(3);
      mean[0] = BE.mean(0);
      mean[1] = BE.mean(1);
      mean[2] = CP.mean(4) + delta;
      
      valarray<double>var(3);
      var[0] = BE.var(0);
      var[1] = BE.var(1);
      var[2] = CP.var(4);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }
  
    if (model=="model34")
      {
      double delta = BE.mean(0) - CP.mean(3);
      delta += BE.mean(1) - CP.mean(4);
      delta /= 2;

      valarray<double>mean(3);
      mean[0] = CP.mean(2) + delta;
      mean[1] = BE.mean(0);
      mean[2] = BE.mean(1);
      
      valarray<double>var(3);
      var[0] = CP.var(2);
      var[1] = BE.var(0);
      var[2] = BE.var(1);

      valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

      BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
      BE1.prop(prop);
      BE1.mean(mean);
      BE1.var(var);
      return BE1;
      }
    }

  // Otherwise
  return directional_impute_01234(opts,model_pair,CP);
  }

// Repeats the final fix-it-up step in the Broad Implementation for 
// only the full model

BroadEstepper1 directional_impute_01234(CanaryOptions& opts,
                                        pair<string,BroadEstepper1> model_pair,
    CanaryPrior CP)
  {
//  valarray<double> freqs = CP.prop();
  string model = model_pair.first;
  BroadEstepper1 BE = model_pair.second;
  
  if (model=="model01234") return BE;

  vector<int>cvec_forced;
  cvec_forced.push_back(0);
  cvec_forced.push_back(1);
  cvec_forced.push_back(2);
  cvec_forced.push_back(3);
  cvec_forced.push_back(4);

  if (model=="model0")
    {
    double delta = BE.mean(0) - CP.mean(0);

    valarray<double> mean(5);
    mean[0] = BE.mean(0);
    mean[1] = CP.mean(1) + delta;
    mean[2] = CP.mean(2) + delta;
    mean[3] = CP.mean(3) + delta;
    mean[4] = CP.mean(4) + delta;
      
    valarray<double> var(5);
    var[0] = BE.var(0);
    var[1] = CP.var(1);
    var[2] = CP.var(2);
    var[3] = CP.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model2")
    {
    double delta = BE.mean(0) - CP.mean(2);

    valarray<double> mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = CP.mean(1) + delta;
    mean[2] = BE.mean(0);
    mean[3] = CP.mean(3) + delta;
    mean[4] = CP.mean(4) + delta;
      
    valarray<double> var(5);
    var[0] = CP.var(0);
    var[1] = CP.var(1);
    var[2] = BE.var(0);
    var[3] = CP.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model01")
    {
    double delta = BE.mean(0) - CP.mean(0);
    delta += BE.mean(1) - CP.mean(1);
    delta /= 2;

    valarray<double>mean(5);
    mean[0] = BE.mean(0);
    mean[1] = BE.mean(1);
    mean[2] = CP.mean(2) + delta;
    mean[3] = CP.mean(3) + delta;
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = BE.var(0);
    var[1] = BE.var(1);
    var[2] = CP.var(2);
    var[3] = CP.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model12")
    {
    double delta = BE.mean(0) - CP.mean(1);
    delta += BE.mean(1) - CP.mean(2);
    delta /= 2;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = BE.mean(0);
    mean[2] = BE.mean(1);
    mean[3] = CP.mean(3) + delta;
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = BE.var(0);
    var[2] = BE.var(1);
    var[3] = CP.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model23")
    {
    double delta = BE.mean(0) - CP.mean(2);
    delta += BE.mean(1) - CP.mean(3);
    delta /= 2;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = CP.mean(1) + delta;
    mean[2] = BE.mean(0);
    mean[3] = BE.mean(1);
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = CP.var(1);
    var[2] = BE.var(0);
    var[3] = BE.var(1);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model34")
    {
    double delta = BE.mean(0) - CP.mean(3);
    delta += BE.mean(1) - CP.mean(4);
    delta /= 2;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = CP.mean(1) + delta;
    mean[2] = CP.mean(2) + delta;
    mean[3] = BE.mean(0);
    mean[4] = BE.mean(1);
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = CP.var(1);
    var[2] = CP.var(2);
    var[3] = BE.var(0);
    var[4] = BE.var(1);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model012")
    {
    double delta = BE.mean(0) - CP.mean(0);
    delta += BE.mean(1) - CP.mean(1);
    delta += BE.mean(2) - CP.mean(2);
    delta /= 3;

    valarray<double>mean(5);
    mean[0] = BE.mean(0);
    mean[1] = BE.mean(1);
    mean[2] = BE.mean(2);
    mean[3] = CP.mean(3) + delta;
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = BE.var(0);
    var[1] = BE.var(1);
    var[2] = BE.var(2);
    var[3] = CP.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model123")
    {
    double delta = BE.mean(0) - CP.mean(1);
    delta += BE.mean(1) - CP.mean(2);
    delta += BE.mean(2) - CP.mean(3);
    delta /= 3;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = BE.mean(0);
    mean[2] = BE.mean(1);
    mean[3] = BE.mean(2);
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = BE.var(0);
    var[2] = BE.var(1);
    var[3] = BE.var(2);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model234")
    {
    double delta = BE.mean(0) - CP.mean(2);
    delta += BE.mean(1) - CP.mean(3);
    delta += BE.mean(2) - CP.mean(4);
    delta /= 3;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = CP.mean(1) + delta;
    mean[2] = BE.mean(0);
    mean[3] = BE.mean(1);
    mean[4] = BE.mean(2);
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = CP.var(1);
    var[2] = BE.var(0);
    var[3] = BE.var(1);
    var[4] = BE.var(2);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model0123")
    {
    double delta = BE.mean(0) - CP.mean(0);
    delta += BE.mean(1) - CP.mean(1);
    delta += BE.mean(2) - CP.mean(2);
    delta += BE.mean(3) - CP.mean(3);
    delta /= 4;

    valarray<double>mean(5);
    mean[0] = BE.mean(0);
    mean[1] = BE.mean(1);
    mean[2] = BE.mean(2);
    mean[3] = BE.mean(3);
    mean[4] = CP.mean(4) + delta;
      
    valarray<double>var(5);
    var[0] = BE.var(0);
    var[1] = BE.var(1);
    var[2] = BE.var(2);
    var[3] = BE.var(3);
    var[4] = CP.var(4);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  if (model=="model1234")
    {
    double delta = BE.mean(0) - CP.mean(1);
    delta += BE.mean(1) - CP.mean(2);
    delta += BE.mean(2) - CP.mean(3);
    delta += BE.mean(3) - CP.mean(4);
    delta /= 4;

    valarray<double>mean(5);
    mean[0] = CP.mean(0) + delta;
    mean[1] = BE.mean(0);
    mean[2] = BE.mean(1);
    mean[3] = BE.mean(2);
    mean[4] = BE.mean(3);
      
    valarray<double>var(5);
    var[0] = CP.var(0);
    var[1] = BE.var(0);
    var[2] = BE.var(1);
    var[3] = BE.var(2);
    var[4] = BE.var(3);

    valarray<double>prop = fill_prop(opts,BE.prop(),BE.cvec(),cvec_forced);

    BroadEstepper1 BE1(opts,BE.vals(),CP,cvec_forced);
    BE1.prop(prop);
    BE1.mean(mean);
    BE1.var(var);
    return BE1;
    }

  Err::errAbort("never should have got here\n");
  return BE;
  }

//  Following Broad algorithm which is not documented
valarray<double> fill_prop(CanaryOptions& opts,
                           valarray<double>prop_in,vector<int>cvec_in,
                           vector<int>cvec_forced)
  {
  double min_prop = opts.tune_min_fill_prop;
  double delta = 0.0;
  unsigned int max_pos = 0;
  double max_prop = 0.0;
  for (unsigned int i=0; i<prop_in.size(); i++)
    {
    if (prop_in[i] > max_prop) { max_pos = i; max_prop = prop_in[i]; }
    if (prop_in[i] < min_prop)
      {
      delta += min_prop - prop_in[i]; 
      prop_in[i] = min_prop; 
      }
    }
  prop_in[max_pos] -= delta;

   
  unsigned int nmissing = cvec_forced.size() - cvec_in.size();
  double fraction_giveaway = opts.tune_fraction_giveaway_0;

  if (nmissing==1) fraction_giveaway = opts.tune_fraction_giveaway_1;
  if (nmissing==2) fraction_giveaway = opts.tune_fraction_giveaway_2;
  if (nmissing==3) fraction_giveaway = opts.tune_fraction_giveaway_3;
  if (nmissing==4) fraction_giveaway = opts.tune_fraction_giveaway_4;

  double per_cluster_need = fraction_giveaway/nmissing; // divide by zero?

  valarray<double>ret_arr(per_cluster_need,cvec_forced.size());
  for (unsigned int i=0; i<cvec_in.size(); i++)
    for (unsigned int k=0; k<cvec_forced.size(); k++)
      if (cvec_in[i] == cvec_forced[k]) 
        {
        ret_arr[k] = (1.0 - fraction_giveaway)*prop_in[i];
        break;
        }

  return ret_arr;
  }

valarray<double> median_polish(Matrix & X, int niter)
	{
	int N = X.Nrows();
	int P = X.Ncols();

	valarray<double>row_med(0.0,N);
	valarray<double>col_med(0.0,P);

	double median;
	double converge;

	int count=0;
	// Probably more samples than probes, so sample(columns) medians first;
	do {
		valarray<double>old_row = row_med;
		valarray<double>old_col = col_med;
		// run through the columns(probes)
	 	for (int j=0; j<P; j++)
			{
			valarray<double> resids(N);
			for (int i=0; i<N; i++)
				resids[i] = X.element(i,j) - row_med[i] - col_med[j];
			col_med[j] += compute_median(resids);
			}

		// The median of the probe(column) effects should be zero so
		// transfer it over to the row(sample) effects.
		median = compute_median(col_med);
		row_med += median;
		col_med -= median;

		// run through the rows(samples)
		for (int i=0; i<N; i++)
			{
			valarray<double> resids(P);
			for (int j=0; j<P; j++)
				resids[j] = X.element(i,j) - row_med[i] - col_med[j];
			row_med[i] += compute_median(resids);
			}

		count++;
		converge=0.0;
		for (int i=0; i<N; i++)
			for (int j=0; j<P; j++)
				{
				converge += abs(old_row[i] - row_med[i]);
				converge += abs(old_col[j] - col_med[j]);
				}

		if (converge < 0.01) break;

		} while (count < niter);

	return row_med;
	}


double compute_median(valarray<double> valvec)
	{
	unsigned int N = valvec.size();
	vector<double>vec(N);
	for (unsigned int i=0; i<N; i++) vec[i] = valvec[i];
	sort(vec.begin(),vec.end());

	if (N % 2) return vec[(N - 1)/2];
	else return 0.5*(vec[N/2 - 1] + vec[N/2]);
	}


string swap_suffix(string name1, string new_suffix)
	{
	name1.replace(name1.rfind("."),name1.size(),new_suffix);
	return name1;
	}


