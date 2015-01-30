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

#ifndef _CANARYWITHPRIOR_H
#define _CANARYWITHPRIOR_H

//
#include "canary/BroadEstepper1.h"
#include "canary/CanaryModel.h"
#include "canary/CanaryOptions.h"
#include "canary/CanaryPenalty.h"
#include "canary/CanaryPrior.h"
//
#include <map>
#include <valarray>
#include <vector>
//

class CanaryWithPrior {
public:
  // Constructed in the cpp file
  CanaryWithPrior(CanaryOptions& opts, 
                  std::valarray<double>& vals_in,
                  CanaryPrior& prior_in);

  valarray<double> vals() { return _vals; }
  CanaryPrior prior() { return _prior; }

  CanaryPenalty penalty(std::string model_name) {
    return CP_map[model_name]; 
    }

  BroadEstepper1 fitted_values(std::string model_name) {
    // cant say this as it creates a BroadEstepper1 if the model name isnt found.
    //return BE_map[model_name];
    std::map<std::string,BroadEstepper1>::iterator i;
    i=BE_map.find(model_name);
    if (i==BE_map.end()) {
      Err::errAbort("Didnt find entry for model: '" + model_name + "'");
    }
    return i->second;
  }

  BroadEstepper1 gmm_withprior(valarray<double> vals,
      std::vector<int> cluster_vec);

  CanaryPenalty assay_fit(BroadEstepper1 BE,std::vector<int> cluster_vec);

  std::string best_fit();

private:
  valarray<double>_vals;
  CanaryPrior _prior;
  std::map<std::string,CanaryPenalty>CP_map;
  std::map<std::string,BroadEstepper1>BE_map;
  CanaryModel _canary_model;
  //
  CanaryOptions& _opts;
};

#endif
