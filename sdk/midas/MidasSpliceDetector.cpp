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

/// @file   MidasSpliceDetector.cpp
/// @brief  Microarray detection of alternative splicing methods.

#include "midas/MidasSpliceDetector.h"
//
#include "stats/statfun.h"
#include "util/Err.h"
//
#include "newmatap.h"
//
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iterator>
//

using namespace std;

/**
 *  \brief Constructor.
 *
 *  Translates group ids into offsets into list of unique ids,
 *  builds full, reduced design matrices.
 *
 *  @param groups Vector of experimental group ids.
 *  @param wantPvalues Caller requests p values.
 *  @param wantFstats Caller requests F statistics.
 *  @param wantNormalized Caller requests normalized exon signal.
 *  @param logStabilize Amount to add to data before taking logarithm.
 *  @param noLogTransform Do not log transform data.
 *
 *  Errors: call errAbort if all requested output flags are false.
*/
midasSpliceDetector::midasSpliceDetector (const std::vector<int>& groups, int wantPvalues,
      int wantFstats, int wantNormalized, const float& logStabilize, int noLogTransform)
    : groups (groups), wantPvalues (wantPvalues), wantFstats (wantFstats),
      wantNormalized (wantNormalized), logStabilize (logStabilize),
      noLogTransform (noLogTransform), N (groups.size())
{
  // if the caller doesn't want any output, (s)he shouldn't bother calling
  if (! (wantPvalues || wantFstats || wantNormalized))
    Err::errAbort("midasSpliceDetector: at least one output is required");

  sampleSetup ();
  buildDesignMatrix ();
  buildSampleOnlyDesignMatrix ();
}

// Debug -- print the input values for comparison to the original
void
dump_vector(std::string vec_label,const std::vector<float>& vec)
{
  int vec_size;
  vec_size=vec.size();
  printf("%s[%d]==",vec_label.c_str(),vec_size);
  for (int i=0;i<vec_size;i++) {
    printf("%f, ",vec[i]);
  }
  printf("\n");
}


/**
 *  \brief Process single exon, gene data.
 *
 *  Matlab: Do a 1-way ANOVA testing for no sample effects.
 *  model tests:
 *  log (exon signal) - log (gene signal) = constant
 *  against
 *  log (exon signal) - log (gene signal) = constant + sample effect
 *
 *  Output:
 *  (1) fstat: for this nexon, the F-statistic on dfModel (reduced design)
 *  numerator degrees of freedom and df1 (full design) denominator degrees of freedom.
 *  (2) normalizedSignal: a vector of nscans normalized exon signal.
 *  Each vector is relative to the nexon, max is 0.
 *  Within groups the values are the same.
 *  (3) pvalue: p-value corresponding to the f statistic.
 *  The output normalizedSignal vector is preallocated by the caller.
 *
 *  Errors: call errAbort if the exon and gene data vectors are not both
 *  equal in size to the count N of group ids (scans).
 *  If normalizedSignal output is requested, call errAbort if the normalizedSignal
 *  vector is not equal in size to the count N of group ids.
 *
 *  @param exonData Exon data.
 *  @param geneData Gene data.
 *  @return fstat F statistic.
 *  @return normalizedSignal Vector of normalized exon signals.
 *  @return pvalue P value.
*/
void midasSpliceDetector::runMidasSingle (const std::vector<float>& exonData,
                                          const std::vector<float>& geneData,
                                          float* fstat,
                                          std::vector<float>* normalizedSignal,
                                          float* pvalue)
{
  // Debug -- print the input values for comparison to the original
  //printf("--------------------\n");
  //dump_vector("exonData",exonData);
  //dump_vector("geneData",geneData);

  // validate vector sizes
  if (! ( exonData.size() == N && geneData.size() == N && (!wantNormalized || normalizedSignal->size() == N) ))
    Err::errAbort("midasSpliceDetector: invalid vector sizes");

  // Matlab: use exon level estimates
  ColumnVector fittedData;
  double residSsq;
  double m1Ssq;
  getEffects (fullDesign, exonData, geneData, fittedData, residSsq, m1Ssq);

  double m2Ssq;
  double modelSsq;
  getEffects (reducedDesign, exonData, geneData, fittedData, m2Ssq, modelSsq);

  // residSsq will be zero if the exonData and geneData vectors are identical
  double myFstat = 0.0;
  if (residSsq != 0.0)
    myFstat = modelSsq/residSsq;

  if (wantFstats)
    *fstat = (float) myFstat;

  if (wantNormalized)
  {
    const double maxFitted = fittedData.Maximum();
    for (unsigned int i = 0; i < N; ++i)
      (*normalizedSignal)[i] = (float) (fittedData (i+1) - maxFitted);
  }

  // convert the F statistic to a p value via the cumulative F distribution function
  if (wantPvalues)
    *pvalue = (float) affxstat::Ftest (myFstat, dfModel, df1);
}

/**
 *  \brief Process multiple exons, single gene data.
 *
 *  Invokes the single exon, gene data runMidas for each
 *  exon in the input vector of vectors.
 *  Produce vector of fstat and pvalue results, vector
 *  of vector of normalizedSignal results.
 *
 *  Errors: call errAbort if the respective output requested flag
 *  is set, and the fstat, normalizedSignal, and pvalue
 *  vectors are not each equal in size to that of the exon vector.
 *
 *  @param exonData Exon data as a vector of vectors.
 *  @param geneData Gene data as a vector.
 *  @return fstat F statistics as a vector.
 *  @return normalizedSignal(s) as a vector of vectors.
 *  @return pvalue P values as a vector.
*/
void midasSpliceDetector::runMidasMultiple (const std::vector<std::vector<float> >& exonData,
                                            const std::vector<float>& geneData,
                                            std::vector<float>* fstat,
                                            std::vector<std::vector<float> >* normalizedSignal,
                                            std::vector<float>* pvalue)
{
  const unsigned int dataSize = exonData.size();
  if (! ( (!wantFstats || fstat->size() == dataSize) && (!wantNormalized || normalizedSignal->size() == dataSize)
      && (!wantPvalues || pvalue->size() == dataSize) ))
    Err::errAbort("midasSpliceDetector: invalid vector sizes");

  // allow caller to pass a null pointer in the case of unwanted output
  for (unsigned int i = 0; i < dataSize; ++i)
  {
    float* pFstat = fstat ? &(*fstat)[i] : 0;
    vector<float>* pNormalizedSignal = normalizedSignal ? &(*normalizedSignal)[i] : 0;
    float* pPvalue = pvalue ? &(*pvalue)[i] : 0;
    runMidasSingle (exonData[i], geneData, pFstat, pNormalizedSignal, pPvalue);
  }
}

/**
 *  \brief Set up sample (group) offsets into unique list of ids.
 *
 *  Matlab: groupIdentifiers is a vector of integers identifying samples.
 *  Biological samples belonging to the same group have the same integer value.

 *  Output: vector containing offset of sample ids in unique sample id list, count of unique sample ids.
 *  Errors: none except allocation - if passed an empty vector, produce an empty vector, zero count.
*/
void midasSpliceDetector::sampleSetup ()
{
  sample_offsets = vector<int> (N);
  vector<int> s_sorted (groups);
  sort (s_sorted.begin(), s_sorted.end());	// sort to identify unique ids

  vector<int> unique_samples;	// unique ids
  unique_copy (s_sorted.begin(), s_sorted.end(), back_inserter(unique_samples));
  T = unique_samples.size();

  // see 'Effective STL' 'Distinguish among count, find, binary_search, lower_bound,
  // upper_bound, and equal_range' for a good discussion on searching containers
  typedef vector<int>::iterator VIIter_t;
  const VIIter_t unique_begin = unique_samples.begin();
  const VIIter_t unique_end = unique_samples.end();
  for (unsigned int i = 0; i < N; ++i)
    sample_offsets[i] = distance (unique_begin, lower_bound (unique_begin, unique_end, groups[i]));
}

/**
 *  \brief Generate full design matrix.
 *
 *  The full design matrix assumes no alternative splicing.
 *  Matlab: Construct the design matrix for M probe effects (for a single exon).
 *      Given N samples with T unique samples, so T sample effects
 *      and all interactions of probes within probeset by tissue.
 *      Samples are an array of integer replicate ids.
 *  Input: vector of sample ids translated into offset into unique id list, count N of samples,
 *         count T of unique ids.
 *  Output: full matrix describing the experimental design.
 *  Note: in the midasSpliceDetector Matlab script, the count M of probe effects = ndataPerPS
 *  is hard-coded to 1 - this simplification is reflected here in the buildDesignMatrix()
 *  and buildSampleOnlyDesignMatrix() methods.
 *  Comment on probe effects count: use an exon-level estimate from Plier not 4 from probe data.
*/
void midasSpliceDetector::buildDesignMatrix ()
{
  // Matlab: design matrix is of full rank since including all effects
  fullDesign = Matrix (N, T);
  fullDesign = 0.0;

  for (unsigned int i = 0; i < N; ++i)
    // Matlab: each interaction matrix for the probe effects x i-th sample effect
    fullDesign (i + 1, (sample_offsets[i] + 1)) = 1.0;
  df1 = fullDesign.Ncols();
}

/**
 *  \brief Generate reduced design matrix.
 *
 *  The reduced design matrix allows for alternative splicing.
 *  Matlab: Given N samples with T unique samples and M probes,
 *      construct the design matrix for T sample effects (for a single exon).
 *      Samples are an array of integer replicate ids.
 *  Input: vector of sample ids translated into offset into unique id list, count N of samples,
 *         count T of unique ids.
 *  Output: Reduced matrix describing the experimental design.
 *  Note: as in buildDesignMatrix(), the count M of probe effects is assumed to be 1.
*/
void midasSpliceDetector::buildSampleOnlyDesignMatrix ()
{
  // Matlab: design matrix is of rank (T-1) for sample effects only
  Matrix sampleEffect (N, T);
  sampleEffect = 0.0;

  // Matlab: i-th sample effect
  for (unsigned int i = 0; i < N; ++i)
    // column is the offset of the sample id in the vector of unique ids
    sampleEffect (i + 1, sample_offsets[i] + 1) = 1.0;

  // Matlab: project onto a subspace orthogonal to the mean and
  // probe effect subspace given the design matrix A for mean and probe effects
  // A = repmat(eye(M), N, 1); design matrix for mean and probe effects
  // B = nulbasis(A); column basis for orthogonal subspace to A
  // P = B*inv(B'*B)*B'is the projection matrix
  // C = colbasis(sampleEffect); C is the column basis of sampleEffect
  // P*C will be the column basis of sampleEffect projected onto B
  // nulbasis relies on Matlab's rref() row-reduced echelon form
  // function; prefer not to implement that, just generate the matrix B
  // eliminate call C = colbasis(sampleEffect) - for a full rank sampleEffect,
  // the case here, colbasis(sampleEffect) equals sampleEffect.
  Matrix B (N, (N - 1));
  for (unsigned int i = 1; i <= N - 1; ++i)	// newmat indices range from 1 to N (here N - 1)
    B (1, i) = -1.0;	// top row all elements -1.0
  IdentityMatrix eyeB (N - 1);
  B.SubMatrix (2, N, 1, (N - 1)) = eyeB;	// the rest of B

  // Matlab: recommended the Matlab way to do this P*C = B*(B\C)
  // i.e. solving Bx = C
  // design = colbasis(B*(B\C));
  // QRZ destroys the input matrices - copy B for later use
  Matrix X = B;
  UpperTriangularMatrix U (N - 1);
  Matrix M (N - 1, T);
  QRZ (X, U);
  QRZ (X, sampleEffect, M);
  Matrix temp = B * U.i() * M;	// design = B*(B\sampleEffect)

  reducedDesign = temp.SubMatrix (1, N, 1, (T - 1));	// what colbasis() does
  dfModel = reducedDesign.Ncols();
}

/**
 *  \brief Find effects for given design matrix, data.
 *
 *  Input: design matrix, vectors of probe and gene data.
 *  Output: fitted data vector (for normalized signal), residual
 *  and model sums of squares, model degrees of freedom.
 *
 *  If so requested by the user, do not log transform the data.
 *
 *  Errors: call errAbort if arguments for log() are negative.
 *  Assume that the vectors are of the appropriate size
 *  (checked in runMidas()).
 *
 *  @param design Design matrix.
 *  @param adjExonData Exon data.
 *  @param geneData Gene data.
 *  @return fittedData Column vector containing fitted data.
 *  @return residSsq Residual sum of squares.
 *  @return modelSsq Model sum of squares.
*/
void midasSpliceDetector::getEffects (const Matrix& design,
                                      const std::vector<float>& exonData,
                                      const std::vector<float>& geneData,
                                      ColumnVector& fittedData,
                                      double& residSsq,
                                      double& modelSsq)
{
  ColumnVector residuals (N);	// residuals are returned here
  for (unsigned int i = 0; i < N; ++i)
  {
    if (noLogTransform)
      residuals (i+1) = exonData[i] - geneData[i];
    else
    {
        if ((exonData[i] + logStabilize) <= 0 || (geneData[i] + logStabilize) <= 0)
            Err::errAbort("Attempt to take the logarithm of a non-positive number");

        // Matlab: log is actually log + small ... logStabilize (default 8)
        residuals (i+1) = log (exonData[i] + logStabilize) - log (geneData[i] + logStabilize);
    }
  }

  // fit = design\(pData-gData);
  Matrix X = design;	// QRZ destroys input matrices
  const int dfModel = design.Ncols();
  UpperTriangularMatrix U (dfModel);
  Matrix M (dfModel, 1);
  QRZ (X, U);
  QRZ (X, residuals, M);
  Matrix fit = U.i() * M;

  // fittedData=design*fit;
  fittedData = design * fit;
  // residuals = pData - gData - fittedData;
  // residuals have already been calculated by QRZ

  // dfModel = rank(design);
  // dfFull = (length(gData) - dfModel);
  // residSsq = residuals' * residuals/dfFull;
  Matrix Ssq = (residuals.t() * residuals) / (N - dfModel);
  residSsq = Ssq (1,1);	// can't directly read product into a double

  // modelSsq = fittedData'*fittedData/dfModel;
  Ssq = (fittedData.t() * fittedData) / dfModel;
  modelSsq = Ssq (1,1);
}
