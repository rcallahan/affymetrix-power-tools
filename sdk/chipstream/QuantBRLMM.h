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
 * @file   QuantBRLMM.h
 * @author Chuck Sugnet
 * @date   Wed Feb 22 09:50:58 2006
 *
 * Implements a SNP genotype calling method motivated by Nusrat Rabbee's
 * RLMM (Robust Linear Model with Mahalanobis Distance Classifier) that utilizes
 * a Baysian prior about the centers of SNP clusters. Idea is to produce a
 * summary value for all probes that measure same genotype (default using RMA to
 * summarize) and then classify the genotype based on the relative intensity of
 * those estimates. For example in genotype AA the A probesets should be
 * relatively bright and B probesets should be relatively dim. For genotype AB
 * both probesets should be lighting up. Exactly how much intensity is expected
 * is learned from a random set of SNPs in the data using high confidence DM
 * calls as reference.
 *
 * For further information about RLMM see:
 * http://www.stat.berkeley.edu/users/nrabbee/RLMM/
 *
 *
 */
/*
   We use floats to be compatible with R.

   Questions:

   - What happens in snp.clusterStats.R::clusterStats() when there are no
   genotypes observe? Looks like a bunch of NA's are passed back, what
   happens to them?
   - Couldn't calculate B summary values for SNP: AFFX-barcodeS
*/

#ifndef QUANTBRLMM_H
#define QUANTBRLMM_H

#define QUANTBRLMM_VERSION "1.0"

//
#include "chipstream/BioTypes.h"
#include "chipstream/ClusterZ.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/TsvReport.h"
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


#if defined (WIN32)
#define asinh(x) (log(x+sqrt(x*x+1)))
#endif

/// String describing quant method
#define QUANTBRLMM "brlmm"

/**
 * Implements a SNP genotype calling method motivated by Nusrat Rabbee's
 * RLMM (Robust Linear Model with Mahalanobis Distance Classifier) that utilizes
 * a Baysian prior about the centers of SNP clusters. Idea is to produce a
 * summary value for all probes that measure same genotype (default using RMA to
 * summarize) and then classify the genotype based on the relative intensity of
 * those estimates. For example in genotype AA the A probesets should be
 * relatively bright and B probesets should be relatively dim. For genotype AB
 * both probesets should be lighting up. Exactly how much intensity is expected
 * is learned from a random set of SNPs in the data using high confidence DM
 * calls as reference.
 *
 * For further information about RLMM see:
 * http://www.stat.berkeley.edu/users/nrabbee/RLMM/
 */
class QuantBRLMM : public QuantGTypeMethod {

public:

  /// What transformation of intensity estimates to feed into classifier.
  enum Transformation {
    MvA, ///< Minus Vs Average
    RvT, ///< R vs Theta
    CES, ///< Contrast Extremes Stretch (previously known as south san francisco).
    CCS  ///< Contrast Centers Stretch (previously known as alternative south san francisco).
  };

  class BrlmmParam {

  public:

    BrlmmParam();

    double m_HetMultiplier;          ///< Number to balance het calls with to balance performance on het/hom calls.
    int m_Iterations;                ///< How many iterations to do of BRLMM, i.e. boostrapping from initial calls.
    double m_IterationThreshold;      ///< What is our threshold (maximum value) for confidences befor
                                     ///  calling 'NN' in iterative mode?
    double m_K;                      ///< Scaling Parameter to use in transformations.
    enum Transformation m_Transform; ///< What transformation of initial data are we feeding into the classifier?
    int m_PriorPseudoCount;          ///< How much weight should the prior get?
    int m_MinCluster;                ///< Minimum examples for each genotype to be included in prior estimation.
    bool m_LowPrecision;             ///< Summary values fed into R are usually truncated, should we do
                                     /// that to be compatible? (used for regression testss)
  };



  /** Constructor, currently creates prior estimates from pre-calculated data from R. */
  QuantBRLMM(double hetMultiplier, int iterations, double iterationThreshold,
             double K, enum Transformation transform,
             int priorWeight, int priorMinCall, bool lowPrecision);

  /** Fill in our self documentation from current state. */
  static void setupSelfDoc(SelfDoc &doc);

  /** Destructor. */
  ~QuantBRLMM();

  /**
   * Tell brlmm to output snp models to file.
   *
   * @param fileName - relative or full path of file to write models to.
   */
  void setSnpModelOut(const std::string& fileName,affx::TsvReport::TsvReportFmt_t format);

  /**
   * Covert string representations into ClusterPrior
   *
   * @param id - identifying string.
   * @param center - centers.
   * @param var - variance of center variances.
   * @param covar - covariance of centers.
   *
   * @return - ClusterPrior.
   */
  static ClusterPrior clusterPriorFromStrings(const std::string& fileName, std::string &id, std::string &center,
                                              std::string &var, std::string &covar);

  /**
   * Write out a text representation of the prior to the file specified.
   * @param prior - Prior to be written.
   * @param out - ostream to write to.
   */
  //  static void writePriorOut(ClusterPrior &prior,affx::TsvFile& tsv,const );
  static void writePriorOut(ClusterPrior &prior, std::ostream &out,
                            const std::string& fieldDelim, const std::string& sep);

  /**
   * Write out a text representation of the model to the file specified.
   * @param model - Model to be written.
   * @param out - ostream to write to.
   */
  void writeModelOut(ClusterModel &model,affx::TsvReport& tsv, const std::string& name,int copies);

  /**
   * Write out a text representation of the prior to the file specified.
   * @param prior
   * @param fileName
   */
  static void writePrior(ClusterPrior &prior, const std::string& fileName);

  /**
   * Covert string representations into ClusterModel
   *
   * @param id - identifying string.
   * @param center - centers.
   * @param var - variance of center variances.
   *
   * @return - ClusterModel
   */
  ClusterModel clusterModelFromStrings(const std::string& fileName, std::string &id,
                                       std::string &center, std::string &var);

  /**
   * Read the prior from the test file specified.
   * @param fileName - Name of file to read from.
   * @return - The prior read in from file.
   */
  static ClusterPrior loadPrior(const std::string& fileName);

  /**
   * Generate a string representation of a snp model.
   * @return string - text description of snp clusters.
   */
  std::string getModelString();

  /**
   * Load in a number of snp specific models.
   * @param snpModels - map of snp probeset names to models.
   * @param fileName - filename to read models from.
   */
  void loadSnpModels(std::map<std::string, ClusterModel> &snpModels,
                     const std::string& fileName);

  /**
   * Load up the snp models from a text file.
   * @param fileName - text file to read models from.
   */
  void readSnpModels(const std::string& fileName);

  /**
   * Load in a number of snp specific priors.
   * @param snpPriors - map of snp probeset names to priors.
   * @param fileName - filename to read priors from.
   */
  static void loadSnpPriors(std::map<std::string, ClusterPrior> &snpPriors, const std::string& fileName);

  /**
   * @brief Set up the quantification method given all the data about the probe
   * set, chip layout and data.
   *
   * @param psGroup - Probes to be used for final estimate.
   * @param layout - Chip layout annotation.
   * @param iMart - Raw data from chips.
   * @param iTrans - Transformations to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

  static void transformData(std::vector<double> &aVals,
                            std::vector<double> &bVals,
                            QuantBRLMM::BrlmmParam &param);

  ClusterPrior *priorForSnp();

  /**
   * @brief Do the heavy lifting of estimation.
   */
  void computeEstimate();

  /**
   * Remove the het cluster from a prior by setting center to be very
   * large negative number and removing covariance with other
   * clusters.
   * @param prior - Prior to remove het cluster from.
   * @param trans - transformation used for clustering.
   */
  void removeHetFromPrior(ClusterPrior &prior, enum Transformation trans);

  /**
   * This is the fix for the snps on chromosome X that don't appear
   * diploid in men (i.e. aren't in the pseudoautosomal region). The
   * concept is to fit and call the model separately for men and
   * women. Implementation is to set all the male initial calls to NN
   * when fitting for women. For fitting and calling men the female
   * samples are set to NN and further the prior is ajusted such that
   * het calls are impossible to get.
   */
  void callHaploidSnp(ClusterPrior &prior);

  static void valueForTransformation(double aVal, double bVal, enum QuantBRLMM::Transformation trans,
                                     double K, double &xVal, double &yVal);

  /**
   * Transform the data in aValues and bValues using the method specified by trans.
   *
   * @param aValues - Summary values for A allele.
   * @param bValues - Summary values for B allele.
   * @param genoTypes - Genotypes for each sample.
   * @param Transformation - Type of transformation to use.
   * @param K - Scale parameter used in CCS and CES transformations.
   * @param xDim - two dimensional matrix for x dimension with rows
   * (index 0) being genotype (0,1,2) and columns being values observed.
   * @param yDim - two dimensional matrix for y dimension with rows
   * (index 0) being genotype (0,1,2) and columns being values observed.
   */
  static void fillInTransformedData(std::vector<double> &aValues, std::vector<double> &bValues,
                                    std::vector<affx::GType> &genoTypes, enum Transformation trans,
                                    double K,
                                    std::vector<std::vector<double> > &xDim, std::vector<std::vector<double> > &yDim);
  /**
   * Specify the prior that will be used for predictions.
   *
   * @param prior - Prior data to use.
   */
  void setPrior(const ClusterPrior &prior) ;

  /**
   * @brief What is the name of the quantification method?
   * @return name of adjuster.
   */
  std::string getType();


  /** Conversion utility for transforming vectors to column vectors. */
  static ColumnVector columnVectorFromVector(std::vector<double> &v);

  /** Conversion utility for transforming arrays to column vectors. */
  static ColumnVector columnVectorFromArray(double *array, int count);

  /** Conversion utility for transforming arrays to matrices. */
  static Matrix matrixFromArray(double *array, int count, int nRow, int nCol);

  /**
   * What version of the BRLMM algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion();

  /** Make banded count matrix. */
  static Matrix makeCountMat(std::vector<int> &counts);

  /** Reformat variance/covariance vector as a banded matrix. */
  static Matrix formatAsBigVarMat(ColumnVector v);

  /** Bayesian update of the mean vector, this is where the voodoo
      happens. Have to query Simon more about the motivation. */
  static Matrix bayesianMeanUpdate(ColumnVector &newVar, std::vector<int> &counts, ColumnVector &means,
                                   Matrix &priorCenterVar, ColumnVector &priorCenter);

  /// @todo should degFree be weight or counts rather than degrees of freedom?
  static ColumnVector shrinkVariance(ClusterStats &empirical, ClusterPrior &prior, int degFree);

  /** Wrapper function for getting bayesian adjusted parameters given data. */
  static ClusterModel fitModel(std::vector<double> &aValues, std::vector<double> &bValues,
                               std::vector<std::vector<double> > &genoDist,
                               std::vector<affx::GType> &genoTypes, ClusterPrior &prior, const std::string &name,
                               BrlmmParam &param);

  /** Count up some elementary statistics about a particular snp over
      multiple experiments. */
  static ClusterStats computeClusterStats(std::vector<double> &aValues, std::vector<double> &bValues,
                                          std::vector<affx::GType> &genoTypes, enum Transformation trans,
                                          double K);
  /**
   * Compute the prior from a collection of snps and fill in the prior.
   *
   * @param snpStats - Stats for a relatively large number of snps
   * @param prior - Fill in the prior.
   */
  void computeClusterPrior(const std::vector<ClusterStats> &snpStats, ClusterPrior &prior);

  /** 
   * Construct a prior from a collection of probesets. Note that to be part of
   * the prior a SNP must have at least 2 examples per AA, AB, and BB otherwise
   * the variance will be zero.
   * 
   * @param probeSets - ProbeSets to be used in computing prior.
   * @param layout - description and chip and probesets.
   * @param iMart - Cel file data for probesets.
   * @param iTrans - Transformers to operate on raw data (i.e. normalization and bg subtraction).
   * @param pmAdjust - Individual probe adjustment like pm-mm or pm-only.
   * 
   * @return - Prior constructed from summary of all the probesets.
   */
  ClusterPrior makePrior(	std::vector<string> &probeSetNames, 
				ChipLayout &layout,
                         	IntensityMart &iMart, 
				std::vector<ChipStream *> &iTrans,
                         	PmAdjuster &pmAdjust);

  /* Make the inverse matrix from the variance and covariances. */
  static Matrix makeVarInv(Matrix &clusterVars, int index);

  /* Take the model derived in fitModel and use it to do predictions. */
  static void brlmmCall(	const std::vector<double> &aValues, 
				const std::vector<double> &bValues,
                        	ClusterModel &model, std::vector<affx::GType> &calls,
                        	std::vector<double> &confidences, 
				std::vector<std::vector<double> > &genoDist,  
				BrlmmParam &param);

  /* Set the quantification method to use for summarizing individual alleles */
  void setQuantExprMethod(QuantExprMethod *qMethod);

  /* Get the quantification method to use for summarizing individual alleles */
  QuantExprMethod *getQuantExprMethod();


  /** 
   * Set the initial guesses of genotypes, usually from DM predictions.
   * 
   * @param genotypes vector of genotypes AA=0, AB=1, AB=2, NN=3,
   */
  void setKnownGenoTypes(std::map<const char *,std::vector<affx::GType>, Util::ltstr> &genotypes);

  /**
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  /**
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();

  /**
   * Convert a (text insensitive) text repesentation of a transformation into
   * the enumration.
   * @param s - string specifying transformation.
   * @return - Enumeration version of transformation.
   */
  static enum Transformation transformationForString(const std::string &s);

  /**
   * Convert a transformation enum to string representation
   * @param t - tranformation enum value.
   * @return - text description of enumeration.
   */
  static std::string stringForTransformation(enum Transformation t);

  /**
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   *
   * @param param - Map of key/value pairs to initialize the object.
   *
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /**
   * Clear out all of our results and values.
   */
  void blankSelf();

  /**
   * How many genotyping calls do we have? Also indicates how many
   * samples there are.
   *
   * @return Number of genotyping calls made.
   */
  size_t getNumCalls();

  /**
   * Get our confidence value for a particular call in a particular sample.
   * @param index - sample of interest.
   * @return - Confidence in the predicted call.
   */
  double getConfidence(unsigned int index);

  /**
   * Get our confidence value for a particular call in a particular sample.
   * @param genoType -  Genotype of interest.
   * @param index - sample of interest.
   * @return - Confidence in the predicted call.
   */
  double getConfGenotype(affx::GType  genoType, unsigned int index);

  /**
   * Get the summary value after transformation for the A and B alleles respectively.
   *
   * @param index - Which chip to get values for.
   * @param aValue - Filled in with transformed a value.
   * @param bValue - Filled in with transformed b value.
   */
  void getAlleleValues(unsigned int index, double &aValue, double &bValue);

  /**
   * Get the summary value names for the A and B alleles respectively.
   *
   * @param aName - Filled in with a value name.
   * @param bName - Filled in with b value name.
   */
  void getAlleleValueNames(std::string &aName, std::string &bName);

  /**
   * Get the name of the probeset that these calls are being made for.
   * @return name of probeset.
   */
  const std::string &getProbeSetName();

  /**
   * Add an expression reporter to output the summary values for each
   * allele, residuals, etc.
   *
   * @param reporter - the object that will output all summary values.
   */
  void addExprReporter(QuantMethodReport *reporter);

  /**
   * Get the parameters used for this brlmm run.
   */
  BrlmmParam getParam();

  void setHaploidSnps(std::map<std::string,bool> &haploidSnps);

  void setSnpSpecificPriors(std::map<std::string,ClusterPrior> &snpPriorMap);

  void setGenders(std::vector<affx::Gender> &genders);

  double getMaxScore();
  void setMaxScore(float v);

  typedef std::map<const char *,std::vector<affx::GType>, Util::ltstr >  KnownGenoTypes_t;
  typedef KnownGenoTypes_t::iterator KnownGenoTypes_iter_t;

  double getMinThresh();
  double getMaxThresh();

  void setDmCallsOut(const std::string& fileName, std::vector<std::string> &col_names);

  virtual bool prepare(const IntensityMart &iMart);

  virtual bool finish();

private:
  std::vector<affx::GType> m_InitialCalls; ///< Genotype calls were are supplied with for learning.
  std::vector<double> m_Confidences; ///< Our resulting confidences in those calls.
  std::vector<std::vector<double> > m_Distances; ///< Our resulting distances for AA, AB, and BB genotypes.
  QuantExprMethod *m_QuantMethod; ///< Quantification method used for summarizing snps.
  std::vector<QuantMethodReport *> m_Reporters;   ///< Report object for outputting results
  ClusterPrior *m_Prior; ///< Prior information about typical snp cluster layouts.
  ClusterPrior *m_CurrentPrior; ///< Prior information about typical snp cluster layouts.
  KnownGenoTypes_t m_KnownGenoTypes; ///< Map from probeset names to the appopriate m_InitialCalls;
  BrlmmParam m_Param; ///< Various parameters used by the brlmm algorithm.
  std::vector<affx::Gender> m_Genders;     ///< What gender is each sample?
  std::map<std::string,bool> m_HaploidSnps; ///< Index of probesets that are in a haploid regions (usually X SNPs).
  std::map<std::string,ClusterPrior> m_SnpPriorMap; ///< Customized priors for different snps indexed by probeset name.

  std::map<std::string,ClusterModel> m_SnpModelMap; ///< Customized models for snps indexed by probeset name.
  std::string m_ModelFileOutName; ///< File to output snp specific models to if that is requested.
  
  double m_BRLMMMaxScore;
  //
  // std::string m_DmCallsFileName; ///< File to output DM calls to if requested.
  //
public:
  affx::TsvReport m_SnpModelsTsv; ///< Output stream for printing snp specific models
  affx::TsvReport m_DmCallsTsv; ///< Output stream for printing snp specific models
};

#endif /* QUANTBRLMM_H */
