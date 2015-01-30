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

#ifndef _INTERFERENCEMODEL_H
#define _INTERFERENCEMODEL_H

//
#include "canary/CanaryOptions.h"
//
#include <valarray>

class InterferenceModel {
public:
  InterferenceModel(CanaryOptions& opts,
                    std::valarray<double> means_in,
                    std::valarray<double> vars_in);

  double interference();

private:
  CanaryOptions& _opts;
  std::valarray<double> _means;
  std::valarray<double> _vars;
};

#endif
