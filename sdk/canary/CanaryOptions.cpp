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

// 
// affy/sdk/canary/CanaryOptions.cpp ---
// 
// $Id: CanaryOptions.cpp,v 1.13 2009-04-13 19:24:05 awilli Exp $
//

//
#include "canary/CanaryOptions.h"
//
#include "util/Convert.h"
#include "util/Util.h"

/////

CanaryOptions::CanaryOptions()
{
  verbosity = 1;
  force = false;
  blockSize = 0;
  precision = 6;
  ccMDChpOutput = true;
  tableOutput = false;
  aptSummarizeAnalysis = "quant-norm.target=1000,pm-only,plier.optmethod=1,expr.genotype=true";

  // set all the tuneables to defaults.
  aptCanaryAnalysis = "af-weight=1.0,TOL=1e-3,hwe_tol=1e-4,hwe_tol2=1e-11,fraction-giveaway-0=0.20,fraction-giveaway-1=0.10,fraction-giveaway-2=0.15,fraction-giveaway-3=0.20,fraction-giveaway-4=0.30,min-fill-prop=0.10,conf-interval-half-width=1.959963984540054,inflation=1.3,min-cluster-variance=0.0001,pseudopoint-factor=100,regularize_variance_factor=0.4";
  // ensure that the defaults are set from this string
  setTuneableStr();
}

CanaryOptions::~CanaryOptions()
{
  // empty
}

//////////

#define SET_TUNE_INT(ARG_STR,ARG_VAR)                                   \
  {                                                                     \
    if (var_name==ARG_STR) {                                            \
      bool ok;                                                          \
      ARG_VAR=Convert::toIntCheck(val,&ok);                             \
      if (!ok) {                                                        \
        Err::errAbort("Cant convert canary tunable value: '"+val+"'."); \
      }                                                                 \
      return;                                                           \
    }                                                                   \
  }
#define SET_TUNE_DOUBLE(ARG_STR,ARG_VAR)                                \
  {                                                                     \
    if (var_name==ARG_STR) {                                            \
      bool ok=true;                                                     \
      ARG_VAR=Convert::toDoubleCheck(val,&ok);                          \
      if (!ok) {                                                        \
        Err::errAbort("Cant convert canary tunable value: '"+val+"'."); \
      }                                                                 \
      return;                                                           \
    }                                                                   \
  }

void
CanaryOptions::setTuneable(const std::string& var_name,const std::string& val)
{
  // for debugging.
  // printf("### CanaryOptions::setTuneable('%s','%s')\n",var_name.c_str(),val.c_str());
  //
  SET_TUNE_DOUBLE("af-weight",tune_af_weight);
  SET_TUNE_DOUBLE("TOL",tune_TOL);
  SET_TUNE_DOUBLE("hwe_tol",tune_hwe_tol);
  SET_TUNE_DOUBLE("hwe_tol2",tune_hwe_tol2);
  //
  SET_TUNE_DOUBLE("fraction-giveaway-0",tune_fraction_giveaway_0);
  SET_TUNE_DOUBLE("fraction-giveaway-1",tune_fraction_giveaway_1);
  SET_TUNE_DOUBLE("fraction-giveaway-2",tune_fraction_giveaway_2);
  SET_TUNE_DOUBLE("fraction-giveaway-3",tune_fraction_giveaway_3);
  SET_TUNE_DOUBLE("fraction-giveaway-4",tune_fraction_giveaway_4);
  //
  SET_TUNE_DOUBLE("min-fill-prop",tune_min_fill_prop);
  //
  SET_TUNE_DOUBLE("conf-interval-half-width",tune_conf_interval_half_width);
  SET_TUNE_DOUBLE("inflation",tune_inflation);
  //
  SET_TUNE_DOUBLE("min-cluster-variance",tune_min_cluster_variance);
  //
  SET_TUNE_DOUBLE("regularize_variance_factor",tune_regularize_variance_factor);
  //
  SET_TUNE_INT("pseudopoint-factor",tune_pseudopoint_factor);
  //
  Err::errAbort("Bad tunable name: '"+var_name+"'"+
                "Please see the help for a list of valid canary tunables." );
}
// dont need anymore.
#undef SET_TUNE_INT
#undef SET_TUNE_DOUBLE

void
//CanaryOptions::setTuneableStr(const std::string& keyval_str)
//{
CanaryOptions::setTuneableStr()
{
//  if (keyval_str=="") {
  if (aptCanaryAnalysis=="") {
    return;
  }

  // make a copy just in case
  std::string keyval_str = aptCanaryAnalysis;
  
  std::vector<std::string> keyval_vec;
  // split on ","
  Util::chopString(keyval_str,',',keyval_vec);
  for (int i=0;i<keyval_vec.size();i++) {
    std::string kv_str=keyval_vec[i];
    Util::trimString(kv_str);
    std::vector<std::string> kv_vec;
    Util::chopString(kv_str,'=',kv_vec);
    //
    if (kv_vec.size()!=2) {
      Err::errAbort("CanaryOptions::setTuneable(): Cant parse: '"+kv_str+"' of '"+keyval_str+"'");
    }
    setTuneable(kv_vec[0],kv_vec[1]);
  }
}
