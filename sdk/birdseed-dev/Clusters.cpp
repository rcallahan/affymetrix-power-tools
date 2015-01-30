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

//
#include "birdseed-dev/Clusters.h"
//
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
//

namespace birdseed
{
namespace dev {
// This class represents the best clustering returned by FitSNPGaussiansPriors3.

Clusters::Clusters(size_t numSamples, size_t k, size_t iteration):
    k(k),
    log_likelihood(-INFINITY),
    means(k),
    vars(k),
    weights(k),
    pxi_zj(numSamples, k, 0.0),
    iteration(iteration)
{
  /*
    std::cout << "k: " << k << "; iteration: " << iteration << "; log_likelihood: " << log_likelihood << "; means: " << matrixToString(means) << "; vars: " << matrixToString(vars) << "; covar: " << covar << "; weights: " << vectorToString(weights) << "\n";
    std::cout.flush();
  */
}

std::string Clusters::tostring() const
{
  std::stringstream strm;
  strm << "k: " << k << "; iteration: " << iteration << "; log_likelihood: " << log_likelihood << "; means: " << matrixToString(means) << "; vars: " << matrixToString(vars) << "; covar: " << covar << "; weights: " << vectorToString(weights);
  return strm.str();
}

std::string Clusters::standardClusterString() const
{
  // Produce a string representation of the cluster in Finny-order
  std::stringstream strm;
  size_t numClusters = (this->means.numRows());
  for (int i = numClusters - 1; i >= 0; --i) {
    double covar = this->covar * sqrt(this->vars[i][0] * this->vars[i][1]);
    strm << ';' <<  this->means[i][0] << ' ' << this->means[i][1] << ' '
    << this->vars[i][0] << ' ' << covar << ' ' << this->vars[i][1] << ' '
    << this->weights[i];
  }
  return strm.str();
}

};
};


