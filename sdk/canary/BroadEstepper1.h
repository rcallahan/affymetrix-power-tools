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

#ifndef _BROAD_ESTEPPER1_H
#define _BROAD_ESTEPPER1_H

#include "canary/CanaryOptions.h"
#include "canary/CanaryPrior.h"
//
#include "newmat.h"
//
#include <valarray>
#include <vector>
//

/// Translation from R to C++ of the Broad's expectation step
/** Expectation step from the R mclust package translated to C++ for the
    case of spherical clusters.  Updating of all parameters according to
    Canary from Birdsuite version 1.3.  For a set of points, each model
    in the model selection scheme will require a BroadEstepper1 instance.
    General idea is to use an instance to update parameter estimates
    by executing the member update_broad().
 */

class BroadEstepper1 {
public:
/// Constructor
/** @param opts stores program execution options, defaults and derived options.
    @param vals_in holds the data points to cluster
    @param prior_in contains prior cluster information
    @clusters_in numeric cluster labels for the model
 */

  BroadEstepper1(CanaryOptions& opts,
                 valarray<double> vals_in,
                 CanaryPrior & prior_in,
                 vector<int> & clusters_in);

  /// Default constructor
  BroadEstepper1(CanaryOptions& opts);

  /// Copy constructor
  BroadEstepper1(const BroadEstepper1 & BE);

  /// Assignment
  BroadEstepper1 & operator=(const BroadEstepper1 & BE);

  /// @return copy of the data
  valarray<double> vals() { return _vals; }

  /// @return copy of the prior
  CanaryPrior prior() { return _prior; }

  /// @return numeric cluster labels
  vector<int> cvec() { return _cvec; }

  /// @return copy of the proportions of cluster memberships
  valarray<double> prop() { return _prop; }

  /// @return value of kth member of cluster memberships
  double prop(int k) { return _prop[k]; }

  /// Breaks encapsulation to assign proportions of cluster memberships
  void prop(valarray<double> prop_in) { _prop = prop_in; }

  /// @return copy of cluster means
  valarray<double> mean() { return _mean; }

  /// @return value of kth cluster mean
  double mean(int k) { return _mean[k]; }

  /// Breaks encapsulation to assign cluster means
  void mean(valarray<double> mean_in) { _mean = mean_in; }

  /// @return copy of cluster variances
  valarray<double> var() { return _var; }

  /// @return value of kth cluster variance
  double var(int k) { return _var[k]; }

  /// Breaks encapsulation to assign cluster variances
  void var(valarray<double> var_in) { _var = var_in; }

  /// @return matrix of probabilities of cluster memberships
  /** (row,col) <=> (i,j) <=> (sample,cluster)
      probability that sample i belongs to cluster j
   */
  Matrix prob() { return _prob; }

//herehere
  void update_broad();
  void update_prop_broad();
  void update_mean_broad();
  void update_var_broad();
  void update_prob();
  double log_likelihood();
  vector<double> confidences();
  vector<int> assignments();

private:
  unsigned long N;
  unsigned long G;
  valarray<double>_vals;
  CanaryPrior _prior;
  vector<int> _cvec;
  valarray<double> _prop;
  valarray<double> _mean;
  valarray<double> _var;
  Matrix _prob;
  CanaryOptions& _opts;
};

#endif
