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

/**
 * @file   ClusterZ.h
 * @author Earl Hubbell
 * @date   Wed Aug 2 2006
 * 
* Strip out the clustering classes from the genotyping routines
* to avoid potential redundancy
* i.e. pretend we have a genotyping SDK
 * * 
 */

#ifndef CLUSTERZ_H 
#define CLUSTERZ_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/ProbeSet.h"
//
#include "stats/stats.h"
//
#include "newmat.h"
//
#include <cfloat>
#include <cstring>
#include <string>
#include <vector>
//

/*
 * this class contains cluster prior data for BRLMM-style clusters
 */
class ClusterPrior {

public: 
  std::string id;
  unsigned int snpCount;
  ColumnVector centers;
  ColumnVector centerVars;
  Matrix covars;
  ClusterPrior() : snpCount(0) {}
};

/*
 * this class contains cluster "model" data for BRLMM-style clusters
 */
class ClusterModel {

public:
  Matrix clusterMeans;
  Matrix clusterVars;
};

/*
 * this class contains summary data for clusters
 */
class ClusterStats {

public:

  /* Set all our counts to three for the three possible genotypes:
     AA,AB,BB */
  ClusterStats();

  // counts an means should usually have 3 observations each one eache
  // for AA, AB, and BB respectively. Check these to make sure they are
  // at least two, otherwise the variance calculations are invalid.
  std::vector<int> counts;   ///< How many of each type of SNP genotype were observed.
  
  // These next two vectors keep track of values for genotyping equivalent of MA plots in
  // R, essentially A is the difference and M is a measure of joint intensity 
  std::vector<double> mMeans; ///< Mean M (intensity) of each genotype (center of the clusters);
  std::vector<double> aMeans; ///< Mean A (difference) of each genotype (center of the clusters);
  std::vector<double> mVars;
  std::vector<double> aVars;
  std::vector<double> covars;
  /// Simon style vars object, not sure why indexed like this but it
  /// goes:
  /// mVars[0],aVars[0],covars[0],mVars[1],aVars[1],covars[1],etc..
  std::vector<double> vars; 
  /// Simon style means object, not sure why indexed like this but it
  /// goes:
  /// mMeans[0],aMeans[0],mMeans[1],aMeans[1],etc..
  std::vector<double> means;
};

// utility functions for dealing with priors

/** utility function: generate cluster from strings */
ClusterPrior x_clusterPriorFromStrings(const std::string& fileName, std::string &id, std::string &center, 
                                       std::string &var, std::string &covar);

/** utility function: write out a ClusterPrior object */
void x_writePriorOut(ClusterPrior &prior, std::ostream &out, 
                     const std::string& fieldDelim, const std::string& sep);

/** utility function: write out a ClusterModel object */
void x_writeModelOut(ClusterModel &model, ostream &out, const std::string &name,
                     int copies, const std::string& fieldDelim, const std::string& sep);

/** utility function: write out a ClusterPrior object */
void x_writePrior(ClusterPrior &prior, const std::string& fileName);

/** utility function: generate a ClusterModel from strings */
ClusterModel x_clusterModelFromStrings(const std::string& fileName, std::string &id, 
                                       std::string &center, std::string &var);

/** load snp specific "priors" from file */
void x_loadSnpPriors(std::map<std::string, ClusterPrior> &snpPriors, 
                     const std::string& fileName);

/** load global prior */
ClusterPrior x_loadPrior(const std::string& fileName);

/** load snp-specific models from file */
void x_loadSnpModels(std::map<std::string, ClusterModel> &snpModels, 
                     const std::string& fileName);

/** produce a string describing a model */
std::string x_getModelString(string m_ProbesetName, ClusterPrior *m_CurrentPrior);

#endif /* CLUSTERZ_H */
