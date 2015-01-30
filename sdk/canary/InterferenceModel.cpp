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

#include "canary/InterferenceModel.h"
//

InterferenceModel::InterferenceModel(CanaryOptions& opts,
                                     std::valarray<double> means_in,
                                     std::valarray<double> vars_in)
  : _opts(opts)
  {
  unsigned int sz = means_in.size();
  _means.resize(sz); _means = means_in;
  _vars.resize(sz); _vars = vars_in;
  }

double InterferenceModel::interference()
  {
  if (_means.size() == 1) return 0.0;

  
  double ci_span = _opts.tune_conf_interval_half_width;

  double total_interval=0.0;
  for (unsigned int k=0; k<_means.size(); k++)
		total_interval += 2.0*ci_span*sqrt(_vars[k]);

  double overlap = 0.0;

  for (unsigned int k=1; k<_means.size(); k++)
    {
    double right_pos = _means[k-1] + ci_span*sqrt(_vars[k-1]);
    double left_pos = _means[k] - ci_span*sqrt(_vars[k]);
	double diff = right_pos - left_pos;
    if (diff > 0) overlap += diff;
    }

  return _opts.tune_inflation*(overlap/total_interval);
  }

