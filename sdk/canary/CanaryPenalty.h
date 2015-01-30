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

#ifndef _CANARYPENALTY_H
#define _CANARYPENALTY_H

// store the penalty here instead of in a data frame found in the R code
class CanaryPenalty {
public:
  CanaryPenalty(); // make the compiler happy (the map<T>)

  CanaryPenalty(const CanaryPenalty & CP);

  CanaryPenalty(double loglik_in, double bic_in, double closeness_in,
      double af_in, double fit_in, double hwe_in, double overlap_in);

  double loglikelihood() { return _loglikelihood; }
  double bic() { return _bic; }
  double closeness() { return  _closeness; }
  double af() { return _af; }
  double fit() { return _fit; }
  double hwe() { return _hwe; }
  double overlap() { return _overlap; }

private:
  double _loglikelihood;
  double _bic;
  double _closeness;
  double _af;
  double _fit;
  double _hwe;
  double _overlap;
};

#endif
