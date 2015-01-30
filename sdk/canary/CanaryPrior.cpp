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

#include "canary/CanaryPrior.h"

CanaryPrior::CanaryPrior(std::string region, std::valarray<double> props,
  std::valarray<double> means, std::valarray<double> vars)
  {
  unsigned int sz = props.size();
  _region = region;
  _prop.resize(sz); _prop = props;
  _mean.resize(sz); _mean = means;
  _var.resize(sz); _var = vars;
  }


// Copy constructor just to make sure that references to vector<T> 
// are not held.  Want actual copies to be safe.
CanaryPrior::CanaryPrior(const CanaryPrior & CP)
  {
  unsigned int sz = CP._prop.size();
  _region = CP._region;
  _prop.resize(sz); _prop = CP._prop;
  _mean.resize(sz); _mean = CP._mean;
  _var.resize(sz); _var = CP._var;
  }


// Assignment 
CanaryPrior & CanaryPrior::operator=(const CanaryPrior & CP)
  {
  if (this == &CP) return *this;
  unsigned int sz = CP._prop.size();
  this->_region = CP._region;
  this->_prop.resize(sz); _prop = CP._prop;
  this->_mean.resize(sz); _mean = CP._mean;
  this->_var.resize(sz); _var = CP._var;
  return *this;
  }
