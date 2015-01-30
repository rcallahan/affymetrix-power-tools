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

// After this many standard deviations set kernel weight to 0
// This will affect the number of observations the kernel spans

#include "copynumber/GKernel.h"

GaussianKernel::GaussianKernel()
{
};
// How many observations lead in front of the kernel center
unsigned int GaussianKernel::lead()
{
  return m_lead;
}

// How may observations lag behind the kernel center
unsigned int GaussianKernel::lag()
{
  return m_lag;
}

// The weights the kernel uses
valarray<double> GaussianKernel::weights()
{
  return m_weights;
}

// The standard deviation the kernel uses for weights
double GaussianKernel::sigma()
{
  return m_sigma;
}

// The number of probes that the kernel spans
unsigned int GaussianKernel::size()
{
  return m_size;
}

GaussianKernel::GaussianKernel(const double sigma_in)
{
  m_sigma = sigma_in;
  unsigned int probe_span = (unsigned int)(G_KERNEL_SPAN * m_sigma);
  m_lead = probe_span;
  m_lag = probe_span;
  m_size = m_lead + m_lag + 1;
  m_weights.resize(m_size);
  m_data.resize(m_size);
  m_data = 0.0;

  double precision = 1.0 / (m_sigma * m_sigma);

  m_weights[m_lag] = 1.0;
  for (unsigned int i = 0; i <= m_lag;i++) {
    double weight = exp(-0.5 * precision * (double)(i * i));
    m_weights[m_lag - i] = weight;
    m_weights[m_lag + i] = weight;
  }

  // If the sum of weights is zero, it makes it quicker to step
  // because the innerproduct will be normalized.
  m_weights /= m_weights.sum();
}

// What happens when a datum is fed into the leading edge of the kernel
double GaussianKernel::step(const double val)
{
  // Accumulate the projection.
  double the_sum = 0.0;

  // Go over all the weights shifting them and taking the inner product
  // at the same time.
  for (unsigned int i = 1; i < m_size; i++) {
    double tmp_val = m_data[i];         // temporary save a read
    m_data[i-1] = tmp_val;              // shift the cached datum
    the_sum += tmp_val * m_weights[i-1];  // accumulate the inner product
  }

  m_data[m_size - 1] = val;            // finish the shift
  the_sum += val * m_weights[m_size-1];  // finish the inner product

  return the_sum;
}
