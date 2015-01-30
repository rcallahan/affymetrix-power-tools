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
 * @file   QuantLabelZ.h
 * @author Earl Hubbell 
 * @date   Oct 30 2006
 */

#ifndef _QUANTLABELZ_H_
#define _QUANTLABELZ_H_

//
#include "chipstream/BioTypes.h"
#include "chipstream/ClusterZ.h"
#include "chipstream/GenoUtility.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodExprReport.h"
//
#include "algorithm/covarnorm/covarnorm.h"
#include "algorithm/em/PrimeEM.h"
#include "algorithm/selector/ProbeSelector.h"
#include "file5/File5.h"
#include "label/snp.label.h"
#include "stats/stats.h"
#include "util/AffxArray.h"
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
#define QUANTBRLMMP "brlmm-p"

//////////

class snp_labeled_distribution
{
public:
    std::string probeset_id;

    snp_distribution Dist;
    //cluster_data bb;
    //cluster_data ab;
    //cluster_data aa;
    //double cv1; = xah
    //double cv2; = xab
    //double cv3; = xhb

    int compareTo(snp_labeled_distribution& that, int iCompareCode)
    {
        int iCompareResult = 0;
        switch (iCompareCode)
        {
        case 0:
            iCompareResult = probeset_id.compare(that.probeset_id);
            break;
        }
        return iCompareResult;
    }

    template<int k> struct ComparePred {
        bool operator()(const snp_labeled_distribution* lhs, const snp_labeled_distribution* rhs) const {
            Err::errAbort("snp_labeled_distribution: ComparePred instantiated with an invalid compare code = " + ToStr(k));
            return false;
        }
    };
};

template<> struct snp_labeled_distribution::ComparePred<0> {
    bool operator()(const snp_labeled_distribution* lhs, const snp_labeled_distribution* rhs) const {
        return lhs->probeset_id.compare(rhs->probeset_id) < 0;
    }
};

//////////

/**
 * This is the primary class wrapping the bayes_label_routine
 * it handles bookkeeping, gender determination, etc
 */
class QuantLabelZ : public QuantGTypeMethod {

public:

  /// What transformation of intensity estimates to feed into classifier.
  
  class BrlmmParam {
    
  public:

    BrlmmParam() {
      m_HetMultiplier = 1;
      m_Iterations = 0;
      m_IterationThreshold = 0.6;
      m_K = 1;
      m_PriorPseudoCount = 20;
      m_Transform = MvA;
      m_MinCluster = 2;
      m_LowPrecision = false;
    }
    
    double m_HetMultiplier;          ///< Number to balance het calls with to balance performance on het/hom calls.
    int m_Iterations;                ///< How many iterations to do of LabelZ, i.e. boostrapping from initial calls.
    double m_IterationThreshold;      ///< What is our threshold (maximum value) for confidences befor 
                                     ///  calling 'NN' in iterative mode?
    double m_K;                      ///< Scaling Parameter to use in transformations.
    enum Transformation m_Transform; ///< What transformation of initial data are we feeding into the classifier?
    int m_PriorPseudoCount;          ///< How much weight should the prior get?
    int m_MinCluster;                ///< Minimum examples for each genotype to be included in prior estimation.
    bool m_LowPrecision;             ///< Summary values fed into R are usually truncated, should we do 
                                     /// that to be compatible? (used for regression tests)
  };

  bool m_OutputProbabilities;
  std::vector<double> m_Txa;         /// probability of AA call
  std::vector<double> m_Txh;         /// probability of het call
  std::vector<double> m_Txb;         /// probability of BB call


  /** Constructor, currently creates prior estimates from pre-calculated data from R. */
  QuantLabelZ( enum Transformation transform,
               double K,
               bool lowPrecision);

  /** Fill in our self documentation from current state. */
  static void setupSelfDoc(SelfDoc &doc);
  /** Destructor. */
  virtual ~QuantLabelZ();

//  virtual void getPrior(string &TmpName, int copyIx, snp_labeled_distribution &objSearch, snp_param &tsp);

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

  void getSequentialPrior(std::string &name, snp_param &tsp, snp_labeled_distribution &sDist);
  
  void getPrior(std::string &TmpName, int presentCopyNumber, snp_param &tsp, snp_labeled_distribution &sDist);

  /**
   * @brief return comma-joined string of call probabilites: BB,AB,AA,Ocean
   */
  std::string getCallProbabilities(int idx);


  /** 
   * @brief Do the heavy lifting of estimation.
   */
  void computeEstimate();

  /**
   * @brief preprocess and normalize values
   */
  void PreprocessValues();
 
  /** 
   * @brief What is the name of the quantification method?
   * @return name of adjuster.
   */  
  std::string getType() {
    return std::string(QUANTBRLMMP);
  }

  /** 
   * What version of the BRLMM-p algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "1.0";
  }

  /** Set the quantification method to use for summarizing individual alleles */
  void setQuantExprMethod(QuantExprMethod *qMethod) {
    m_QuantMethod = qMethod;
  }

  /** Get the quantification method to use for summarizing individual alleles */
  QuantExprMethod *getQuantExprMethod() {
    return m_QuantMethod;
  }

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
  virtual inline size_t getNumCalls() { return m_Calls.size(); }

  /** 
   * Get our confidence value for a particular call in a particular sample.
   * @param index - sample of interest.
   * @return - Confidence in the predicted call.
   */
  double getConfidence(unsigned int index) { 
    if(index >= m_Confidences.size()) {
      Err::errAbort("Asking for call at index " + ToStr(index) + " when Probeset " + 
                    m_ProbesetName + " has only " + ToStr(m_Confidences.size()) + " calls.");
    }
    return m_Confidences[index];
  }

  /** 
   * Returns standardized distances to cluster centers for given genotype.
   * @param genoType -  Genotype of interest.
   * @param index - sample of interest.
   * @return - Standardized distance to cluster center.
   */
  virtual double getConfGenotype(affx::GType genoType, unsigned int index) {
      ///@todo we should either implement this, or not output 
      ///      it in QuantMethodReportGTypeChipsummary
    assert(index < m_Distances[0].size());
    if(genoType == affx::NN)
      return FLT_MAX;
   // this matrix is reversed from BRLMM for convenience
   // in computing the distances vs subsets of the data
   // which occurs in chrX/Y cases in males/females
    return m_Distances[genoType][index];

    //return 0;
  }
 
  /** 
   * Get the summary value after transformation for the A and B alleles respectively.
   * 
   * @param index - Which chip to get values for.
   * @param aValue - Filled in with transformed a value.
   * @param bValue - Filled in with transformed b value.
   */
  void getAlleleValues(unsigned int index, double &aValue, double &bValue) {
    assert(index < m_AValues.size());
    assert(m_AValues.size() == m_BValues.size());
    aValue = m_AValues[index];
    bValue = m_BValues[index];
  }

  /**
   * Transform allele estimates
   */
  void transform() {
      GenoUtility_transformData(m_AValues, m_BValues, m_Param.m_Transform,m_Param.m_K);
  }
 /**
   * Set probeset selection data as 'probes'
   *
   * @param aVals - vector of allele A estimates
   * @param bVals - vector of allele B estimates
   */
  void setProbeSetSelect(std::vector<double> &aVals, std::vector<double> &bVals){
    const size_t NTarg = aVals.size();
    int NFeat = 1; // just the probeset
        m_PM_A.clear();
        m_PM_A.push_back(vector<double>(NTarg));
     
        m_PM_B.clear();
        m_PM_B.push_back(vector<double>(NTarg));

        for (size_t i=0; i< NTarg; i++){
            m_PM_A[0][i] = aVals[i];
            m_PM_B[0][i] = bVals[i];
        }
        m_A_probeId.assign(NFeat,0);
        m_B_probeId.assign(NFeat,0);
  }
  /**
   * Set allele values
   * 
   * @param aVals - vector of allele A estimates
   * @param bVals - vector of allele B estimates
   */
  void setAlleleValues(std::vector<double> &aVals, std::vector<double> &bVals) {
      m_AValues = aVals;
      m_BValues = bVals;
      if (m_SnpProbeTsv.is_open())
        setProbeSetSelect(aVals,bVals);
  }

  /**
   * Set override flag
   * 
   * @param oflag - select an override option for hints
   *
   */
  void setOverride(int oflag) {
      m_SelectionOverride = oflag;
      setOptValue("override", ToStr(oflag));
  }


  inline void setOutputProbabilities(bool outFlag) {
      m_OutputProbabilities = outFlag;
  }
 
  /** 
   * Get the summary value names for the A and B alleles respectively.
   * 
   * @param aName - Filled in with a value name.
   * @param bName - Filled in with b value name.
   */
  void getAlleleValueNames(std::string &aName, std::string &bName)
  {
      pair<string,string> labels = columnNamesForTransformation(m_Param.m_Transform);
      aName = labels.first;
      bName = labels.second;
  }

  /** 
   * Get the name of the probeset that these calls are being made for. 
   * @return name of probeset.
   */
  const std::string &getProbeSetName() { return m_ProbesetName; }

  /** 
   * Set the probeset name
   * @param probeset name
   */
  void setProbeSetName(const std::string &probesetName) {
      m_ProbesetName = probesetName;
  }

  /** 
   * Add an expression reporter to output the summary values for each
   * allele, residuals, etc.
   * 
   * @param reporter - the object that will output all summary values.
   */
  void addExprReporter(QuantMethodReport *reporter) {
    m_Reporters.push_back(reporter);
  }

  /** 
   * Get the BRLMM-like parameters used for this run.
   */ 
  BrlmmParam getParam() {
    return m_Param;
  }
  /**
   * Set the snp-parameters for this run.
   */
  void setLabelParam(snp_param s){
    sp.copy(s);

    // keep our options in sync so that the state is
    // correctly reported
    s.setDocValues(*this);
  }

  /**
   * Get the parameters for this run of brlmm-p.
   */
  snp_param getLabelParam(){
    return sp;
  }

  /**
   * makes the normalization functions for all samples
   *
   * @param probeSets - my sample of snps to determine normalization function by intensity
   * @param layout - finding those probes on the chip
   * @param iMart - raw intensities
   * @param iTrans - transformations
   * @param pmAdjust - adjustments, such as background subtraction
   */
  vector<NormalizationPredictor> makeNorm(std::vector<std::string> probeSets, ChipLayout &layout,
                                          const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                                          PmAdjuster &pmAdjust); 

  /**
   * sets the normalization functions for all samples
   * 
   * @param Qmap - normalization functions
   */
  void setNorm(vector<NormalizationPredictor> &Qmap)
  {
    normMap.resize(Qmap.size());
    for (unsigned int i=0; i<normMap.size(); i++)
    {
      normMap[i].Copy(Qmap[i]);
    }
  }

  /** returns the max-score used for this analysis method */
  double getMaxScore(){ return(LabelzMaxScore);};
  /** sets the max-score used for this analysis method */
  void setMaxScore(double v){
    LabelzMaxScore=v;

    // keep our options in sync so that the state is
    // correctly reported
    setOptValue("MS", ToStr(LabelzMaxScore));
  };

  /**
   * which snps are haploid snps
   *
   * @param haploidSnps - map of string to chrX status 
   */
  void setHaploidSnps(std::map<std::string,bool> &haploidSnps) {
    m_HaploidSnps = haploidSnps;
  }

  /**
   * which snps are special
   *
   * @ param SpecialSnps - map of string to copy status
   */
  void setSpecialSnps(std::map<std::string,std::pair<int, int> > &SpecialSnps){
    m_SpecialSnps = SpecialSnps;
  }

  /**
   * which snps are special and sample-specific
   *
   * @ param SpecialSampleSnps - map of string to vector of copy status by sample
   */
  void setSpecialSampleSnps(std::map<std::string,std::vector<int> > &SpecialSampleSnps,
                            std::vector<std::string> &SpecialSampleNames){
    m_SpecialSampleSnps = SpecialSampleSnps;
    m_SpecialSampleNames = SpecialSampleNames;
  }

  /**
   * sets the genders for all samples
   *
   * @param genders - gender vector describing samples
   */
  void setGenders(const std::vector<affx::Gender> &genders) {
    m_Genders = genders;
  }

  void setInbredHetPenalty(const std::vector<double> &InbredPenalty){
	  m_InbredHetPenalty = InbredPenalty;
   }

  /**
   *  Extracts intensities, etc for a probe set, so further analysis can proceed
   *
   *  @param gtPs - pointer to probe set
   * @param layout - layout of probesets on array
   * @param iMart - method to grab raw intensities and serve them
   * @param iTrans - transformations of intensities
   * @param pmAdjust - background subtraction, etc
   * @param doReport - do we invoke the output reporters for summary values
   */
  bool SetUpProbeSet(const ProbeSet *gtPs,
                     const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                     PmAdjuster &pmAdjust,bool doReport);

  void SubsetSetup(vector<int> &samplesubset, int &subsetlength, int copyIx, vector<int>& SnpSampleCopyNumber);

  /**
   * Determines the gender of samples from chrX snps by EM.
   * 
   * @param genders - vector of genders for samples that this routine fills in
   * @param probeSets - sample of snps that happen to be chrX
   * @param layout - layout of probe sets on array
   * @param iMart - serves up raw intensities
   * @param iTrans - modifies the raw intensities
   * @param pmAdjuster - subtracts background
   */
  void FindGender(std::vector<affx::Gender> &genders, std::vector<ProbeSet *> probeSets, 
                    ChipLayout &layout, const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                    PmAdjuster &pmAdjust);

  /** 
   * Set the initial guesses of genotypes, usually from DM predictions.
   * 
   * @param genotypes vector of genotypes AA=0, AB=1, AB=2, NN=3,
   */
  void setKnownGenoTypes(std::map<const char *,std::vector<affx::GType>, Util::ltstr> *genotypes) {
    m_KnownGenoTypes = genotypes;
  }

  /**
   * Sets the EM specific parameters for calling genders from chrX snps
   *
   * @param em_t - threshold for em
   * @param em_c - cutoff for em
   * @param g_c - cutoff for choosing genders
   */
  void setGenderParam(double em_t,double em_c,double g_c)
  {
    em_thresh = em_t;
    em_cutoff = em_c;
    gender_cutoff = g_c;

    // keep our options in sync so that the state is
    // correctly reported
    setOptValue("em_thresh", ToStr(em_thresh));
    setOptValue("em_cutoff", ToStr(em_cutoff));
    setOptValue("gender_cutoff", ToStr(gender_cutoff));
  }

  /**
   * Selects good probes
   * 
   */
  void SelectGoodProbes();

  /**
   * Load up the initial call for the requested SNP
   * 
   * @param probeset - name of the probeset for which to load known genotypes
   */
  void loadInitialCallsForSnp(const string &probeset);

  typedef std::map<const char *,std::vector<affx::GType>, Util::ltstr >  KnownGenoTypes_t;
  typedef KnownGenoTypes_t::iterator KnownGenoTypes_iter_t;

  double getMinThresh() { return -DBL_MAX;}
  double getMaxThresh() { return LabelzMaxScore;}
  virtual bool prepare(const IntensityMart &iMart);
  virtual bool finish();

  void setCopyNumber(bool b) {
    m_bCopyNumber = b;
    setOptValue("copynumber",ToStr(b));
  }

  virtual void setParameters(PsBoard &board);

  // ArtifactReduction
  //do a dummy trick which will explode memory painfully
  //avert by loading in only probeset-id's that are relevant?
  //
  // given a sample:probeset:block trio, tell me the confidence
  std::map<string, int> m_ProbeSetTrust;  
  void EliminateUntrustedSamples(vector< int >&SnpSampleCopyNumber);
  void loadProbeSetTrustFromTmpFile();

  /**
   * Set TrustCheck flag
   * 
   * @param tc_flag - indicates whether to look for trusted probe file from
   * ArtifactReduction
   */
  void setTrustCheck(bool tc_flag) {
      m_bProbeSetTrustOpt = tc_flag;
      setOptValue("TrustCheck", ToStr(tc_flag));
  }

protected:
  std::string m_SnpPosteriorName;
  affx::TsvReport::TsvReportFmt_t m_PosteriorFormat;
  std::string m_SequentialModelFile;
  affx::TsvFile m_SequentialModelTsv;
  std::string m_Id, m_BB, m_AB, m_AA, m_CV;
  bool m_bCopyNumber; /// This meethod is begin called as part of a copynumber run.
  QuantExprMethod *m_QuantMethod; ///< Quantification method used for summarizing snps.
  BrlmmParam m_Param; ///< Various parameters used by the brlmm algorithm.
  std::vector<QuantMethodReport *> m_Reporters;   ///< Report object for outputting results
  //  std::map<std::string, snp_param> m_SnpPriorMap; ///< customized priors indexed by probeset name
  //AffxArray<snp_posterior> m_vectSnpPriors;
  snp_param sp; ///< parameters used by labeling
  std::vector<double> m_Confidences; ///< Our resulting confidences in those calls.
  AffxArray<snp_labeled_distribution> m_vectSnpPriors;
  std::vector<std::vector<double> > m_Distances; ///< standardized distances from AA, AB, and BB cluster centers
  std::map<std::string,std::pair<int, int> > m_SpecialSnps; ///< index of special snps and their copy numbers
  std::vector<affx::Gender> m_Genders;     ///< What gender is each sample?
  std::set<const char *, Util::ltstr> *m_ProbeSetsToReport;

private:

  double LabelzMaxScore; ///< parameter used for making no-calls when doing LabelZ

  // goofy I/O stuff
  // read >in< prior
  // process data 
  // >output< posterior which contains SNP specific tuning
  // when you read it in, it's a prior
  // when you output it, it's a posterior

  // normalization map
  // vector of normalization predictors
  // one per sample
  vector< NormalizationPredictor> normMap; ///< vector of normalization functions per sample

  std::map<std::string, vector<double> > m_SnpCovarMap; ///< customized covariates indexed by probeset name

  // handle haploid Snps and gender issues
  std::map<std::string,bool> m_HaploidSnps; ///< Index of probesets that are in a haploid regions (usually X SNPs).
  std::map<std::string,std::vector<int> > m_SpecialSampleSnps; ///< index of special sample-specific snps and their copy numbers
  std::vector<std::string> m_SpecialSampleNames; ///< index of sample names for special sample-specific snps
  double em_thresh, em_cutoff,gender_cutoff; ///< three parameters for calling genders

  // this is used as "hints" for generating the likelihood
  // these are not >necessary< as for BRLMM
  // but can be >helpful< in evaluating likelihood or generating specific priors
  // i.e. work from partial data to classify full data
  KnownGenoTypes_t *m_KnownGenoTypes; ///< Map from probeset names to the appopriate m_InitialCalls;
  std::vector<affx::GType> m_InitialCalls; ///< Genotype calls we are supplied with for learning.

  // Fix up for genotyping variably inbred samples
  // likelihood penalty for deviating from normal sample behavior
  // positive for penalizing hets, negative for increasing hets
  std::vector<double> m_InbredHetPenalty; ///< Sample specific covariate like gender, but affecting all SNPs

  // store the individual probe data for selection
  std::vector<std::vector<double> > m_PM_A; ///< individual intensities A
  std::vector<std::vector<double> > m_PM_B; ///< individual intensities B
  std::vector< unsigned int > m_A_probeId; ///< individual probe indexes A
  std::vector< unsigned int > m_B_probeId; ///< individual probe indexes B
  int m_SelectionOverride; ///< override calls with references for selection

public:

  // Artifact Reduction
  bool isTrusted(int);
  void EliminateUntrustedSamples();
  void readProbesetTrust(const std::string &fileName);
  std::string m_ProbeSetTrustTmpFileName;
  bool m_bProbeSetTrustRead;
  bool m_bProbeSetTrustOpt;
  bool m_ZWGenderCalling;
  ///< 
  affx::TsvReport m_SnpPosteriorTsv;

  int m_SnpPosteriorTsv_ver;
  ///<
  affx::TsvReport m_NormSummaryTsv;
  ///< 
  affx::TsvReport m_SnpProbeTsv;

  /**
   * writes out the contrast normalization functions
   *
   * @param nMap - function for normalization
   * @param fileName - which file to write to
   */
  static void writeNormMap(vector<NormalizationPredictor> &nMap, 
                    const std::string& fileName);

  ////
  void writeSnpPosteriorTsv(const std::string& fileName,
                            affx::TsvReport::TsvReportFmt_t format);
  ///
  void writeSnpPosteriorTsv(const std::string& fileName,
                            affx::TsvReport::TsvReportFmt_t format,
                            int ver);
  ///
  void writeSnpPosteriorValue(const std::string& probeset_id, snp_param &tsp);
 
  std::string getModelFile() {
    return m_SnpPosteriorTsv.getFilePath();
  }

  /**
   * open the stream for snp-specific probe selection
   */
  void readSnpProbeTsv(const std::string& fileName,
                       affx::TsvReport::TsvReportFmt_t format);
  void writeSnpProbeTsv(const std::string& fileName,
                        affx::TsvReport::TsvReportFmt_t format);

  void writeSnpProbeValue(affx::TsvReport& tsv,
                          int probenum,
                          unsigned int Aid,
                          unsigned int Bid,
                          const vector<double> &Summaries,
                          const vector<double> &m_Ax,
                          const vector<double> &m_Bx,
                          string& m_ProbeSetName);

  /**
   * open file for normalized values
   */
  void readNormSummaryTsv(const std::string& fileName, 
                          std::vector<std::string> &headers,
                          affx::TsvReport::TsvReportFmt_t format);
  void writeNormSummaryTsv(const std::string& fileName, 
                           std::vector<std::string> &headers,
                           affx::TsvReport::TsvReportFmt_t format);
  void writeNormSummaryValue(const std::string& m_ProbeSetName,
                             const std::vector<double> &m_AValues,
                             const std::vector<double> &m_BValues);
  /**
   * read in the snp specific priors from a file
   */
  void readSnpPriorMap(const std::string& fileName);
  void readSnpPriorMap_tsv5(affx::File5_Tsv* tsv5);

  virtual void registerProbeSetsToReport(std::set<const char *, Util::ltstr> *probeSetsToReport) {
    m_ProbeSetsToReport = probeSetsToReport;
  }


};

#endif /* _QUANTLABELZ_H_ */
