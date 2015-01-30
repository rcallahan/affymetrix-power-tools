////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   SnpClusterStore.h
 * @author Chuck Sugnet
 * @date   Tue Feb  2 11:05:21 2010
 * 
 * @brief Class for reading the posterior from snp
 * clustering. Currently supports only Axiom but goal is to support
 * all GTC Affymetrix Genotyping algorithms. Goal is to provide
 * relatively quick random access to snp posteriors for display in
 * cluster plots.
 * 
 */
#ifndef SNPCLUSTERSTORE_H

#include "label/snp.label.h"
#include "chipstream/SnpModelDb.h"
#include "chipstream/QuantLabelZ.h"

#include <vector>
#include <string>

/**
 * Class for reading the posterior from snp clustering. Currently
 * supports only Axiom but goal is to support all GTC Affymetrix
 * Genotyping algorithms. Goal is to provide relatively quick random
 * access to snp posteriors for display in cluster plots.
 * 
 */
class SnpClusterStore {

public:

  /** Enum for types supported. */
  enum ClusterModelType {
    UNKNOWN,
    AXIOM
  };

  /** 
   * Constructor that will open up the hdf5 file supplied and index probe sets by position.
   *  
   * @param fileName - Name of the model file of interest (eg AxiomGT1.snp-posteriors.a5)
   * @param tmpFileName - Optional name of temp file to write if necessary.
   */
  SnpClusterStore(const std::string &fileName, const std::string &tmpName="");
  
  /**
   * Clean up
   */
  ~SnpClusterStore();

  /** 
   * Return true if the model exists for this snp probeset for this copynumber.
   * 
   * @param snpName - Name of snp probeset (eg AX-11088371)
   * @param copyNumber - Usually 2, but can be 1 for chr X snps in men or chrY snps in men
   * 
   * @return - true if there is a model available, false otherwise
   */
  bool snpClusterExists(const std::string &snpName, int copyNumber) const;
  
  /** 
   * Given the name of a snp probeset and the copynumber desired
   * return the cluster distributions in a snp_distribution
   * object. For example the Axiom SNP has copynumber 2 in women but
   * copynumber 1 in men. This may be extended to other duplications
   * and deletions in the future.
   * 
   * @param snpName - Name of snp probeset (eg AX-11088371)
   * @param copyNumber - Usually 2, but can be 1 for chr X snps in men or chrY snps in men
   * 
   * @return - object with cluster information filled in.
   */
  snp_distribution getSnpCluster(const std::string &snpName, int copyNumber) const;

  /** 
   * Get the values associated with a particular key in the header information.
   * 
   * @param key - Name of the value to search for.
   * @param values - Values associated with the key specified.
   */
  void getHeaderValue(const std::string &key, std::vector<std::string> &values) const;
  
private:

    SnpModelDb m_SnpDb;
    std::string m_TmpFile;
    mutable snp_labeled_distribution m_DistCache;

};

#endif /* CLUSTER_STORE */
