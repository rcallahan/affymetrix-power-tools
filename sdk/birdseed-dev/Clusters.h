////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

#ifndef _CLUSTERS_H
#define _CLUSTERS_H

//
#include "birdseed-dev/Matrix.h"
#include "birdseed-dev/birdseeddefs.h"
//
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
//

#ifdef _MSC_VER
#include <limits>
#define INFINITY std::numeric_limits<double>::infinity()
#endif

#ifdef __sun__
#define INFINITY std::numeric_limits<double>::infinity()
#endif

namespace birdseed
{
namespace dev {
// This class represents the best clustering returned by FitSNPGaussiansPriors3.
class Clusters
{
public:
  // Number of clusters found
  size_t k;
  // Goodness of this clustering
  double log_likelihood;
  // Ultimately these have MAX_NUM_CLUSTERS rows, but not initially
  // 0th element is for AA, 1st is for AB, 2nd is for BB.
  VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> means;
  VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> vars;
  VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> weights;
  double covar;
  // matrix of prob of a sample being in a particular cluster NxK
  PXI_ZJMatrix pxi_zj;

  // Which iteration of algorithm produced this result.
  size_t iteration;

  Clusters(size_t numSamples, size_t k, size_t iteration);
  std::string tostring() const;

  std::string standardClusterString() const;

};

};
};


#endif /* _CLUSTERS_H */

/******************************************************************/
/**************************[END OF Clusters.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
