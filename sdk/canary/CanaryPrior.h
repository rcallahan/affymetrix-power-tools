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

#ifndef _CANARYPRIOR_H
#define _CANARYPRIOR_H

//
#include <cstring>
#include <string>
#include <valarray>
//

class CanaryPrior {

public:
  // Check the default constructors to see if the empty fields might be
  // filled in other than empty.  Otherwise let the compiler do it.
  CanaryPrior() { }


  CanaryPrior(std::string region, std::valarray<double> props,
      std::valarray<double> means, std::valarray<double> vars);


  // Copy constructor just to make sure that references to vector<T> 
  // are not held.  Want actual copies to be safe.
  CanaryPrior(const CanaryPrior & CP);

  // Assignment 
  CanaryPrior & operator=(const CanaryPrior & CP);

  std::string region() { return _region; }
  std::valarray<double>prop() { return  _prop; }
  double prop(int k) { return  _prop[k]; }

  std::valarray<double>mean() { return  _mean; }
  double mean(int k) { return  _mean[k]; }

  std::valarray<double>var() { return _var; }
  double var(int k) { return  _var[k]; }
  void var(std::valarray<double>var_in) { _var = var_in; }

private:
  std::string _region;
  std::valarray<double> _prop;
  std::valarray<double> _mean;
  std::valarray<double> _var;
};

#endif
