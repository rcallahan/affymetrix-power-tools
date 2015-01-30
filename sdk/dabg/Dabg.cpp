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

/// @file   Dabg.cpp
/// @author harley
/// @brief  A c++ version of the Detect Above BackGround code.

//
#include "dabg/Dabg.h"
//
#include "portability/affy-base-types.h"
#include "stats/stats-distributions.h"
//
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>

/// @brief Creates a new Dabg object and inits it
affx::Dabg::Dabg()
{
  clear();
}

/// @brief               Creates a new Dabg object with space reserved.
/// @param chip_cap      reserve space for this many chips
/// @param intensity_cap reserve space for this many intensities per chip
affx::Dabg::Dabg(int capacity)
{
  clear();
  intensity_reserve(capacity);
}

/// @brief Initialize the object to empty.
void
affx::Dabg::clear()
{
  is_prepared=0;

  data.resize(AFFYDABG_GC_COUNT_MAX+1);
  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    data[gcbin].clear();
  }
}

/// @brief          reserve space for intensity values
/// @param chip_idx the chip to reserve space for.
void
affx::Dabg::intensity_reserve(int capacity)
{
  // capacity is the total number of probes
  // each gc_vector will have a fraction of the probes.
  int gcbin_capacity_guess;
  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    // estimate capacitys for the bins.
    if ((gcbin<5)||(gcbin>19)) {
      gcbin_capacity_guess=capacity*0.0005;
    }
    else if ((gcbin<10)||(gcbin>16)) {
      gcbin_capacity_guess=capacity*0.0800;
    }
    else {
      gcbin_capacity_guess=capacity*0.1600;
    }
    //
    data[gcbin].reserve(gcbin_capacity_guess);
  }
}

size_t
affx::Dabg::intensity_count(int gc_idx)
{
  return data.at(gc_idx).size();
}

/// @brief     Copy the reservations from another dabg
/// @param     dabg_from     Dabg to copy from
void
affx::Dabg::intensity_reserve_copy(affx::Dabg& dabg_from)
{
  data.resize(AFFYDABG_GC_COUNT_MAX+1);
  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    data[gcbin].reserve(dabg_from.data[gcbin].size());
  }
}

/// @brief add a vector of background intensities
/// @param chip_idx      the index of the chip to add to
/// @param gc_vec        a vector of gc content of the intensity values
/// @param intensity_vec a vector of the intensity values 
void
affx::Dabg::add_intensitiy_vec(std::vector<int> gc_vec,std::vector<float> intensity_vec)
{
  assert(gc_vec.size()==intensity_vec.size());

  is_prepared=0; // not ready

  uint32_t gc_vec_size=gc_vec.size();
  for (uint32_t i=0;i<gc_vec_size;i++) {
    data[gc_vec[i]].push_back(intensity_vec[i]);
  }
}

/// @brief             Add a single intensity
/// @param chip_idx    the chip to add it to.
/// @param gcbin       the gc value of this probe
/// @param intensity   the intensity value of this probe
void
affx::Dabg::add_intensity(int gcbin,affx::Dabg::intensity_t intensity)
{
  is_prepared=0; // not ready
  data[gcbin].push_back(intensity);
}

//! prepare this object for use.
//! At the moment, this just means sorting the vectors.
//! also see prepare_maybe in "Dabg.h"
/// @brief   prepare the object for use (sort all the intensities)
void 
affx::Dabg::prepare()
{
  // printf("PREPARE!!!\n"); // 
  for (uint32_t gcbin=0;gcbin<data.size();gcbin++) {
    sort(data[gcbin].begin(),data[gcbin].end());
  }
  is_prepared=1; // ready
}

/// @brief            Convert an intensity to a p-value
/// @param chip       which chip to use
/// @param gcbin      which gc bin the probe is in
/// @param intensity  the intensity of the probe
/// @return           the probablity that this p-value is above the background.
/// @remark
//! this is the left-tail probablity.  
//! Prob=value is larger than all items to the left.
//! vectors are sorted in increasing order.
double
affx::Dabg::intensity_to_pvalue(int gcbin,Dabg::intensity_t intensity)
{
  double theval,pval;
  //Dabg::intensity_t val;
  // insure data is prepared (ie sorted)
  prepare_maybe();

  int l_idx;
  //int u_idx;
  l_idx=0;

  int icnt=data[gcbin].size();
  //printf("icnt=%d\n",icnt);
  // no values?
  if (icnt==0) {
    pval=0; // this will be inverted.
  } 
  // one value?
  else if (icnt==1) {
    theval=data[gcbin][0];
    if (intensity<theval) {
      pval=0.0; // perl dabg returns 1; we invert
    }
    else if (intensity>theval) {
      pval=1.0; // perl dabg returns 0; we invert
    }
    else {
      pval=0.5;
    }
  } 
  else {
    // l_idx is the index of the value ">=" intensity
    l_idx=std::lower_bound(data[gcbin].begin(),data[gcbin].end(),intensity)
      - data[gcbin].begin();
    // fix it up to be less than intensity
    if (l_idx==icnt) {
      l_idx=icnt-1;
    }
    pval=((double)l_idx)/((double)icnt);
  }

  // now that we have our pval, invert it for our chances of being greater.
  pval=(1.0-pval);
  
  //printf("dabg: intensity_to_pvalue(%2d,%2d,%8.4f)=%8.4f\n",chip,gcbin,intensity,pval);
  //printf("dabg: n=%3d : in=%3d : %10.6f : pos=%3d : val=%3d\n",
  //       (int)data[gcbin].size(),(int)intensity,pval,l_idx,(int)data[gcbin][l_idx]);
  return pval;
}

/// @brief       Dump the contents of the object
/// @param level The level of verbosity
/// @remark
//! level=0: summary
//! level=1: gcbins with sizes>0
//! level=2: all gcbins
//! level=3: all data
void
affx::Dabg::dump(int level)
{
  printf("Dabg dump (%d) ==========\n",level);

  if (level==0) { return; }

  for (uint32_t gcbin=0;gcbin<data.size();gcbin++) {
    int v_size=(int)data[gcbin].size();
    int v_cap=(int)data[gcbin].capacity();
    if ((level>=2)||((level>=1)&&(v_size>0))) {
      printf("# data[gcbin=%2d]==%d/%d\n",gcbin,v_size,v_cap);
      if (level>=3) {
        for (uint32_t i=0;i<data[gcbin].size();i++) {
          printf("%3d: %f\n",i,(double)data[gcbin][i]);
        }
      }
    }
  }
}

//////////////////////////////

/// @brief       Compute the targets and pvals for all the chips.
/// @param probe_cnt   number of probes
/// @param gc_arr      array of probe gc values
/// @param pm_arr      array of probe pm intensities
/// @param target_arr  output array of target probabilities
/// @param pvalues_arr output array of probe probabilities
//void
//Dabg::rundabg(int probe_cnt,
//              int* gc_arr,double** pm_arr,
//              double* target_arr,double** pvalues_arr)
//{
//  for (int idx=0;idx<cel_cnt;idx++) {
//    rundabg_idx(idx,probe_cnt,gc_arr,pm_arr[idx],&target_arr[idx],pvalues_arr[idx]);
//  }
//}

/// @brief          Compute the target and pvals for one chip.
/// @param chip_idx    Index of the chip
/// @param probe_cnt   number of probes in set
/// @param gc_arr      gc values of the probles
/// @param pm_arr      pm intensity values of the probles
/// @param target      Output: probability of the probset above background.
/// @param pval_arr    Output: probability of the probe above background
void
affx::Dabg::compute_target_fisher(int probe_cnt,int* gc_arr,double* pm_arr,int dof_adjustment,
                                  double* out_target,double* out_pval_arr)
{
  double pval_product;

  // fill pvalues_arr and compute product
  pval_product=0.0;
  for (int pi=0;pi<probe_cnt;pi++) {
    out_pval_arr[pi]=intensity_to_pvalue(gc_arr[pi],(affx::Dabg::intensity_t)pm_arr[pi]);
    pval_product=pval_product+log(out_pval_arr[pi]);
    //printf("dabg2: pi=%2d pm=%12.8f pval=%12.8f  product=%12.8f\n",pi,pm_arr[pi],pval[pi],pval_product);
  }

  /// @todo This should be replaced with a constant bias rather than a step function.
  if (pval_product==0.0) {
    pval_product=0.000001; // copied from perl
  }
  
  /// @todo  Why '-2*'?
  double test_statistic=-2*pval_product;
  *out_target=affxstat::chisqrprob(2*probe_cnt - dof_adjustment,(float)test_statistic);
  //printf("dabg2: test_statistic=%12.8f\n",test_statistic);
  //printf("dabg2: target=%12.8f\n",*target);
}

void
affx::Dabg::compute_target_percentile(int probe_cnt,int* gc_arr,double* pm_arr,
                                      double per_cut,
                                      double* out_target,double* out_pval_arr)
{
    // in bounds
  assert(per_cut>=0.0);
  assert(per_cut<=1.0);
  //
  std::vector<double> pval_vec;
  pval_vec.reserve(probe_cnt);

  // fill pvalues_arr and compute product
  for (int pi=0;pi<probe_cnt;pi++) {
    // save one to return
    out_pval_arr[pi]=intensity_to_pvalue(gc_arr[pi],(affx::Dabg::intensity_t)pm_arr[pi]);
    // save one to work with.
    pval_vec.push_back(out_pval_arr[pi]);
    //printf("pval_arr[%d]=%f\n",pi,out_pval_arr[pi]);
  }

  // only one value, so return it as the answer.
  if (pval_vec.size()==1) {
    *out_target=pval_vec[0];
  }
  // 
  else {
    // the order is ascending...
    sort(pval_vec.begin(),pval_vec.end());
    // so we need to invert per_cutoff
    //per_cut=1.0-per_cut;
    // now lookup the pvalue (and interpolate)
    double pos=per_cut*(probe_cnt-1);
    int    idx=(int)pos;
    *out_target=pval_vec[idx]+((pos-idx)*(pval_vec[idx+1]-pval_vec[idx]));
    //printf("target=%f\n",*out_target);
  }
}

/// @todo Turn these into ">>" and "<<" operators.

void
affx::Dabg::dabgfile_write(const std::string& filename)
{
  std::fstream file;

  file.open(filename.c_str(),std::fstream::out);
  file << "# This is a dabg data file. Format: (bin[0-25] count \\n (value \\n)* \\n)*\n";

  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    file << gcbin << "    " << data[gcbin].size() << "\n";
    int i_max=data[gcbin].size();
    for (int i=0;i<i_max;i++) {
      file << data[gcbin][i] << "\n";
    }
    file << "\n";
  }
    
  file.close();
}


void
affx::Dabg::dabgfile_read(const std::string& filename)
{
  std::fstream file;
  char line[1000];

  // zero out the contents
  clear(); 

  file.open(filename.c_str(),std::fstream::in);
  // discard first line of the file
  file.getline(line,sizeof(line));

  //
  int gcbin_in;
  int gcbin_size;
  float val;
  for (int gcbin=0;gcbin<=AFFYDABG_GC_COUNT_MAX;gcbin++) {
    file >> gcbin_in;
    assert(gcbin_in==gcbin);
    file >> gcbin_size;
    data[gcbin].reserve(gcbin_size);
    for (int i=0;i<gcbin_size;i++) {
      file >> val;
      data[gcbin].push_back(val);
    }
  }

  file.close();
}
  
