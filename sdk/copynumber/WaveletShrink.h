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
#ifndef WAVELETSHRINK_H
#define WAVELETSHRINK_H

#include "../../external/pywavelets/common.h"
#include "../../external/pywavelets/convolution.h"
#include "../../external/pywavelets/wavelets.h"
#include "../../external/pywavelets/wt.h"

#include <utility>
#include <valarray>
#include <vector>

void wavelet_sandbox(std::valarray<double> & observed, Wavelet * W,
	std::vector<double> likelihood_prec, double prior_prec,
	MODE mode, double converge, bool outlier, double nu, int maxiter);

#endif
