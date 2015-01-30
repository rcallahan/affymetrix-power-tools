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

/// @file   MidasSpliceDetector.h
/// @brief  Headers for MidasSpliceDetector.

#ifndef MIDASSPLICEDETECTOR_H
#define MIDASSPLICEDETECTOR_H

// This #undef of New is required because of a conflict with the perl include
// file handy.h, which #define's New; New() is also a method in newmat.
// handy.h is included in the swig wrapper for MidasSpliceDetector.
#include "newmat.h"
//
#include <vector>
//
#ifdef New
#undef New
#endif /* New */

/**
 *  midasSpliceDetector
 *  \brief Object for Midas, Microarray Detection of Alternative Splicing.
 *
 *  Performs one-way ANOVA to detect alternative spliced cassette
 *  exons from exon array data.  See the white paper 'Alternative
 *  Analysis Methods for Exon Arrays' for a detailed description
 *  of the method.
 *
 *  The constructor midasSpliceDetector() has input parameter a vector
 *  of integers, e.g., [ 1 1 1 2 2 3 ], one for each microarray scan,
 *  describing the experimental groups.  These can be in any order
 *  and need not be consecutive; [ 2 1 1 5 7 7 3 ] is valid.
 *  The order of group ids matches the order of the exon
 *  and gene data passed to runMidas.
 *
 *  Subsequent data analysis can be performed either on a single
 *  or multiple exon basis.  For processing a single exon,
 *  the caller provides a vector of floats, one entry per scan, of
 *  exon data, and another vector of floats of gene data.
 *  These may be produced by an earlier Plier run and must be
 *  of the same size as the groups vector. RunMidas calculates
 *  the F statistic, normalized signal, and p value for this exon.
 *  The normalized signal, based on a log value, is a normalized measure of
 *  how much exon specific expression differs between groups;
 *  it is returned as a vector of floats, of the same dimension
 *  as the exon and gene input.  The midasSpliceDetector
 *  output is returned by reference to lvalues allocated by the caller.
 *
 *  For processing multiple exons, the caller provides a vector
 *  of vectors of exon data, one vector for each exon,
 *  and one vector of gene data; all exons are compared with
 *  the single gene data vector.  Output is returned as a vector
 *  of F statistics, one entry per exon, a vector of vectors of
 *  normalized exon signal, one vector per exon of length the number
 *  of scans, and a vector of p values, one entry per exon.
 *  Here also, output is returned by reference to lvalues allocated
 *  by the caller.
 *
 *  The caller may choose not to receive one or more of the outputs;
 *  this is signaled by flags passed to the constructor.
 */
class midasSpliceDetector
{
public:
  /** Constructor.
   * @param groups Experimental groups.
   * @param wantPvalues Caller requests p values.
   * @param wantFstats Caller requests F statistics.
   * @param wantNormalized Caller requests normalized exon signal.
   * @param logStabilize Amount to add to data before taking logarithm.
   * @param noLogTransform Do not log transform data.
   */

  /* Note: the parameters wantPvalues, wantFstats, wantNormalized, and noLogTransform,
     user configurable options, are functionally either true or false and could be
     treated as bools.  Here they are defined as ints so that this code could be
     swig wrapped; in Visual Studio builds of swig wrapped code, using bools for these
     parameters produced unresolved external symbol errors at link time. */

  midasSpliceDetector(const std::vector<int>& groups,
                      int wantPvalues,
                      int wantFstats,
                      int wantNormalized,
                      const float& logStabilize,
                      int noLogTransform);

  /** Calculates stats on single exon, gene data.
   * @param exonData         Single exon data.
   * @param geneData         Gene data.
   * @param fstat            F statistic.
   * @param normalizedSignal Vector of normalized exon signal.
   * @param pvalue           P value.
   */
  void runMidasSingle (const std::vector<float>& exonData,
                       const std::vector<float>& geneData,
                       float* fstat, 
                       std::vector<float>* normalizedSignal,
                       float* pvalue);

  /** Calculates stats on multiple exons, gene data.
   * @param exonData         Multiple exon data.
   * @param geneData         Gene data.
   * @param fstat            returned F statistic vector
   * @param normalizedSignal returned vector of vectors of normalized exon signal
   * @param pvalue           returned P values vector
   */
  void runMidasMultiple (const std::vector<std::vector<float> >& exonData,
                         const std::vector<float>& geneData,
                         std::vector<float>* fstat,
                         std::vector<std::vector<float> >* normalizedSignal,
                         std::vector<float>* pvalue);

private:

  /** Sets up sample (group) offsets into unique list of ids.
   */
  void sampleSetup ();

  /** Generates full design matrix.
   */
  void buildDesignMatrix ();

  /** Generates reduced design matrix.
   */
  void buildSampleOnlyDesignMatrix ();

  /** Finds effects for given design matrix, data.
   * @param design        Design matrix.
   * @param adjExonData   Exon data.
   * @param geneData      Gene data.
   * @param fittedData    Column vector containing fitted data.
   * @param residSsq      Residual sum of squares.
   * @param modelSsq      Model sum of squares.
   */
  void getEffects (const Matrix& design,
                   const std::vector<float>& exonData,
                   const std::vector<float>& geneData,
                   ColumnVector& fittedData, 
                   double& residSsq,
                   double& modelSsq);

  /// private data
  /// experimental groups
  const std::vector<int> groups;
  /// caller requests p values
  const int wantPvalues;
  /// caller requests F statistics
  const int wantFstats;
  /// caller requests normalized exon signal
  const int wantNormalized;
  /// Amount to add to data before taking logarithm
  const float logStabilize;
  /// caller requests no log transform of data
  const int noLogTransform;
  /// offset of sample ids in vector of unique vector ids
  std::vector<int> sample_offsets;
  /// total count of samples
  const unsigned int N;
  /// count of unique sample ids
  int T;
  /// full design matrix
  Matrix fullDesign;
  /// rank of full design matrix
  double df1;
  /// reduced design matrix
  Matrix reducedDesign;
  /// rank of reduced design matrix
  double dfModel;
};

#endif /* MIDASSPLICEDETECTOR_H */
