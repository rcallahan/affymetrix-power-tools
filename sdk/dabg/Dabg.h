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

/// @file   Dabg.h
/// @author harley
/// @brief  The headers for the Detect Above BackGround SDK.
#ifndef DABG_H
#define DABG_H

#include <cstring>
#include <string>
#include <vector>
//

//! AFFYDABG_GC_COUNT_MAX the max gc count
#define AFFYDABG_GC_COUNT_MAX (25)

namespace affx {
  class Dabg;
  class DabgGroup;
}

/// @todo template <typename Intensity_t>
class affx::Dabg
{
public:
  typedef float intensity_t;

  //
  int chip_capacity;
  int intensity_capacity;
  int is_prepared;

  //! where the data is held. dabg[gc_count][idx]=intensity
  std::vector< std::vector<intensity_t> > data;

  //
  Dabg();
  Dabg(int intensity_cap);

  //
  void clear();

  //
  void intensity_reserve(int approx_probe_cnt);
  void intensity_reserve_copy(affx::Dabg& dabg2);

  // The number of intensities in this gc bin.
  size_t intensity_count(int gc_bin);

  //
  void add_intensitiy_vec(std::vector<int> gc_vec,std::vector<intensity_t> intensity_vec);
  void add_intensity(int gc_idx,intensity_t intensity);
  //
  double intensity_to_pvalue(int gc_idx,intensity_t intensity);
  void dump(int level);
  //
  void prepare();
  void prepare_maybe() { if (is_prepared!=1) { prepare(); } };

  // two ways of computing the target value
  void compute_target_fisher(int probe_cnt,int* gc_arr,double* pm_arr, int dof_adjustment,
                             double* out_target,double* out_pval_arr);
  void compute_target_percentile(int probe_cnt,int* gc_arr,double* pm_arr,
                                 double per_cut,
                                 double* out_target,double* out_pval_arr);
  //
  void dabgfile_write(const std::string& filename);
  void dabgfile_read(const std::string& filename);
};

#endif /* DABG_H */
