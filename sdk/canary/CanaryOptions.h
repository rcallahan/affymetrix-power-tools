////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#ifndef CANARYOPTIONS_H
#define CANARYOPTIONS_H

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/AnalysisStream.h"
#include "chipstream/QuantExprMethod.h"
#include "util/Util.h"
//
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//

/** Place to group program options. */
class CanaryOptions {

public:
  // constructor
  CanaryOptions();
  ~CanaryOptions(); 

  //
  void setTuneable(const std::string& key,const std::string& val);
  void setTuneableStr();

  std::string help;
  int verbosity;
  std::string aptSummarizeAnalysis;
  std::string aptCanaryAnalysis;
  bool force;
  std::string setAnalysisName;
  int precision;
  std::string outDir;
  bool tableOutput;
  bool ccMDChpOutput;
  std::string progName;
  int blockSize;
  std::string version;
  std::string cvsId;
  std::string commandLine;
  uint64_t freeMemAtStart;

  unsigned int numRows;
  unsigned int numCols;
  unsigned int probeCount;
  std::string chipType;
  std::string mapName;
  std::string mapVersion;

  std::map<std::string,std::string> stdMethods;

  std::string timeStr;           // Date and time of execution.  
  std::string execGuid;          // Execution id.
  std::string analysisGuid;      // Execution id.
  std::string cdfFile;           // Name of cdf-file.
  std::string spfFile;           // Name of cdf-file.
  std::string genomeVersion;     // verion of genome from regions file

  std::string cnvRegionFile;     // CNV Region.
  std::string cnvMapFile;        // CNV Map.
  std::string cnvNormFile;       // Normalization probes.
  std::string canaryPriorFile;

  std::vector<std::string> celFiles;       ///< Cel files to analyze.
  std::vector<std::string> celGuids;       ///< Cel files to analyze.
  std::vector<std::string> chpFiles;       ///< CHP files to analyze.

//  std::string missingProbeFile;

  // These "tunables" were magic numbers in the R code.
  // they can be adjusted at runtime with setTuneable(name,val);
  double tune_af_weight;
  double tune_TOL;
  double tune_hwe_tol;
  double tune_hwe_tol2;
  //
  double tune_fraction_giveaway_0;
  double tune_fraction_giveaway_1;
  double tune_fraction_giveaway_2;
  double tune_fraction_giveaway_3;
  double tune_fraction_giveaway_4;
  //
  double tune_min_fill_prop;
  //
  double tune_conf_interval_half_width;
  double tune_inflation;
  //
  double tune_min_cluster_variance;
  double tune_regularize_variance_factor;
  //
  int    tune_pseudopoint_factor;
};

#endif
