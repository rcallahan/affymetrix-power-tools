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

#include "canary/CanaryPenalty.h"

// store the penalty here instead of in a data frame really a fancy struct

CanaryPenalty::CanaryPenalty()
  {
  _loglikelihood = 0.0;
  _bic = 0.0;
  _closeness = 0.0;
  _af = 0.0;
  _fit = 0.0;
  _hwe = 0.0;
  _overlap = 0.0;
  }


CanaryPenalty::CanaryPenalty(const CanaryPenalty & CP)
  {
  _loglikelihood = CP._loglikelihood;
  _bic = CP._bic;
  _closeness = CP._closeness;
  _af = CP._af;
  _fit = CP._fit;
  _hwe = CP._hwe;
  _overlap = CP._overlap;
  }


CanaryPenalty::CanaryPenalty(double loglik_in, double bic_in,
    double closeness_in, double af_in, double fit_in, double hwe_in,
    double overlap_in)
  {
  _loglikelihood = loglik_in;
  _bic = bic_in;
  _closeness = closeness_in;
  _af = af_in;
  _fit = fit_in;
  _hwe = hwe_in;
  _overlap = overlap_in;
  }
