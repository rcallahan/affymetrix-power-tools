////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#ifndef GAUSSIAN_KERNEL_H
#define GAUSSIAN_KERNEL_H

// After this many standard deviations set kernel weight to 0
// This will affect the number of observations the kernel spans
#define G_KERNEL_SPAN 2.0

#include <valarray>

using namespace std;

class GaussianKernel
{
public:
  // Ready to configure after declaration. But no way to do it yet.
  GaussianKernel();

  // Construct the weights inside the constructor.  This could be moved
  // outside of the constructor.
  // Preconfigured instance
  GaussianKernel(const double sigma_in);

  // How many observations lead in front of the kernel center
  unsigned int lead();

  // How may observations lag behind the kernel center
  unsigned int lag();

  // The weights the kernel uses
  valarray<double>weights();

  // The standard deviation the kernel uses for weights
  double sigma();

  // The number of probes that the kernel spans
  unsigned int size();
  // What happens when a datum is fed into the leading edge of the kernel
  double step(const double val);

private:
  unsigned int m_lead;
  unsigned int m_lag;
  unsigned int m_size;
  double m_sigma;
  valarray<double>m_weights;
  valarray<double>m_data;
};


#endif
