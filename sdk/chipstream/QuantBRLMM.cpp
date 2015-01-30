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
 * @file   QuantBRLMM.cpp
 * @author Chuck Sugnet
 * @date   Fri Feb 24 11:02:09 2006
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
 */

//
#include "chipstream/QuantBRLMM.h"
//
#include "chipstream/BioTypes.h"
#include "chipstream/MatrixUtil.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/QuantMethodFactory.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <sstream>

using namespace affx;

QuantBRLMM::BrlmmParam::BrlmmParam() {
  m_HetMultiplier = 1;
  m_Iterations = 0;
  m_IterationThreshold = 0.6;
  m_K = 1;
  m_PriorPseudoCount = 20;
  m_Transform = MvA;
  m_MinCluster = 2;
  m_LowPrecision = false;
}

/** Constructor, currently creates prior estimates from pre-calculated data from R. */
QuantBRLMM::QuantBRLMM(double hetMultiplier, int iterations, double iterationThreshold,
                       double K, enum QuantBRLMM::Transformation transform,
                       int priorWeight, int priorMinCall, bool lowPrecision) {
  m_Param.m_HetMultiplier = hetMultiplier;
  m_Param.m_Iterations = iterations;
  m_Param.m_K = K;
  if (iterationThreshold < 0 || iterationThreshold > 1) {
    Err::errAbort("QuantBRLMM() - iter-thresh must be >= 0 and <= 1");
  }
  m_Param.m_IterationThreshold = iterationThreshold;
  if (transform > CCS || transform < MvA) {
    Err::errAbort("QuantBRLMM() - Expecting transform to be between: " + ToStr(MvA) + " (MvA) and: " +
                  ToStr(CCS) + " (CCS), but got: " + ToStr(transform));
  }
  if (priorWeight < 0) {
    Err::errAbort("QuantBRLMM() - Expecting a non-negative prior weight.");
  }
  m_Param.m_PriorPseudoCount = priorWeight;
  if (priorMinCall < 2) {
    Err::errAbort("QuantBRLMM() - Must have at least 2 examples per genotype for "
                  "use in prior generation. Current priorMinCall is: " + ToStr(priorMinCall));
  }
  m_Param.m_MinCluster = priorMinCall;
  m_Param.m_LowPrecision = lowPrecision;
  m_Param.m_Transform = transform;
  m_QuantMethod = NULL;
  m_Prior = NULL;
  m_CurrentPrior = NULL;
  m_BRLMMMaxScore = -1;

  setupSelfDoc(*this);

  setOptValue("het-mult", ToStr(m_Param.m_HetMultiplier));
  setOptValue("iterations", ToStr(m_Param.m_Iterations));
  setOptValue("iter-thresh", ToStr(m_Param.m_IterationThreshold));
  setOptValue("K", ToStr(m_Param.m_K));
  setOptValue("transform", stringForTransformation(m_Param.m_Transform));
  setOptValue("prior-weight", ToStr(m_Param.m_PriorPseudoCount));
  setOptValue("prior-mincall", ToStr(m_Param.m_MinCluster));
  setOptValue("lowprecision", m_Param.m_LowPrecision);
}

/** Fill in our self documentation from current state. */
void QuantBRLMM::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(QUANTBRLMM);
  doc.setDocDescription("Do genotyping calls using the BRLMM (Bayesian RLMM) algorithm.");
  doc.setDocOptions(getDefaultDocOptions());
}

/** Destructor. */
QuantBRLMM::~QuantBRLMM() {
  if (m_SnpModelsTsv.is_open()) {
    m_SnpModelsTsv.close();
  }
  if (m_DmCallsTsv.is_open()) {
    m_DmCallsTsv.close();
  }
  for (std::vector<QuantMethodReport *>::iterator it = m_Reporters.begin(); it != m_Reporters.end(); it++)
    {
      delete *it;
    }
  m_Reporters.clear();
  Freez(m_QuantMethod);
}

/**
 * Load up the snp models from a text file.
 * @param fileName - text file to read models from.
 */
void QuantBRLMM::readSnpModels(const std::string& fileName) {
  loadSnpModels(m_SnpModelMap, fileName);
}

void QuantBRLMM::transformData(std::vector<double> &aVals,
                               std::vector<double> &bVals,
                               QuantBRLMM::BrlmmParam &param) {
  assert(aVals.size() == bVals.size());
  double x = 0, y = 0;
  for (unsigned int i = 0; i < aVals.size(); i++) {
    valueForTransformation(aVals[i], bVals[i], param.m_Transform, param.m_K, x, y);
    aVals[i] = x;
    bVals[i] = y;
  }
}

/**
 * Specify the prior that will be used for predictions.
 *
 * @param prior - Prior data to use.
 */
void QuantBRLMM::setPrior(const ClusterPrior &prior) {
  delete m_Prior;
  m_Prior = new ClusterPrior();
  m_Prior->id = prior.id;
  m_Prior->covars = prior.covars;
  m_Prior->centers = prior.centers;
  m_Prior->centerVars = prior.centerVars;
}

/**
 * @brief What is the name of the quantification method?
 * @return name of adjuster.
 */
std::string QuantBRLMM::getType() {
  return std::string(QUANTBRLMM);
}

/* Set the quantification method to use for summarizing individual alleles */
void QuantBRLMM::setQuantExprMethod(QuantExprMethod *qMethod) {
  Freez(m_QuantMethod);
  m_QuantMethod = qMethod;
}
/* Get the quantification method to use for summarizing individual alleles */
QuantExprMethod *QuantBRLMM::getQuantExprMethod() {
  return m_QuantMethod;
}

/**
 * What version of the BRLMM algorithm are we implementing?
 * @return version string.
 */
std::string QuantBRLMM::getVersion() {
  return QUANTBRLMM_VERSION;
}

/** 
 * Set the initial guesses of genotypes, usually from DM predictions.
 * 
 * @param genotypes vector of genotypes AA=0, AB=1, AB=2, NN=3,
 */
void QuantBRLMM::setKnownGenoTypes(std::map<const char *,std::vector<affx::GType>, Util::ltstr> &genotypes) {
  m_KnownGenoTypes = genotypes;
}

/**
 * Clear out all of our results and values.
 */
void QuantBRLMM::blankSelf() {
  clearProbeSet(m_Aallele);
  clearProbeSet(m_Ballele);
  m_AValues.clear();
  m_BValues.clear();
  m_InitialCalls.clear();
  m_Calls.clear();
  m_Confidences.clear();
  m_Distances.clear();
};

/**
 * How many genotyping calls do we have? Also indicates how many
 * samples there are.
 *
 * @return Number of genotyping calls made.
 */
size_t QuantBRLMM::getNumCalls() { return m_Calls.size(); }

/**
 * Get our confidence value for a particular call in a particular sample.
 * @param index - sample of interest.
 * @return - Confidence in the predicted call.
 */
double QuantBRLMM::getConfidence(unsigned int index) {
  if (index >= m_Confidences.size()) {
    Err::errAbort("Asking for call at index " + ToStr(index) + " when Probeset " +
                  m_ProbesetName + " has only " + ToStr(m_Confidences.size()) + " calls.");
  }
  return m_Confidences[index];
}

/**
 * Get our confidence value for a particular call in a particular sample.
 * @param genoType -  Genotype of interest.
 * @param index - sample of interest.
 * @return - Confidence in the predicted call.
 */
double QuantBRLMM::getConfGenotype(affx::GType  genoType, unsigned int index) {
  assert(index < m_Distances.size());
  if (genoType == affx::NN)
    return FLT_MAX;
  return m_Distances[index][genoType];
}

/**
 * Get the summary value after transformation for the A and B alleles respectively.
 *
 * @param index - Which chip to get values for.
 * @param aValue - Filled in with transformed a value.
 * @param bValue - Filled in with transformed b value.
 */
void QuantBRLMM::getAlleleValues(unsigned int index, double &aValue, double &bValue) {
  assert(index < m_AValues.size());
  assert(m_AValues.size() == m_BValues.size());
  aValue = m_AValues[index];
  bValue = m_BValues[index];
}

/**
 * Get the summary value names for the A and B alleles respectively.
 *
 * @param aName - Filled in with a value name.
 * @param bName - Filled in with b value name.
 */
void QuantBRLMM::getAlleleValueNames(std::string &aName, std::string &bName)
{
  aName = "Contrast";
  bName = "Strength";
}

/**
 * Get the name of the probeset that these calls are being made for.
 * @return name of probeset.
 */
const std::string &QuantBRLMM::getProbeSetName() { return m_ProbesetName; }

/**
 * Add an expression reporter to output the summary values for each
 * allele, residuals, etc.
 *
 * @param reporter - the object that will output all summary values.
 */
void QuantBRLMM::addExprReporter(QuantMethodReport *reporter) {
  m_Reporters.push_back(reporter);
}

/**
 * Get the parameters used for this brlmm run.
 */
QuantBRLMM::BrlmmParam QuantBRLMM::getParam() {
  return m_Param;
}

void QuantBRLMM::setHaploidSnps(std::map<std::string,bool> &haploidSnps) {
  m_HaploidSnps = haploidSnps;
}

void QuantBRLMM::setSnpSpecificPriors(std::map<std::string,ClusterPrior> &snpPriorMap) {
  m_SnpPriorMap = snpPriorMap;
}

void QuantBRLMM::setGenders(std::vector<affx::Gender> &genders) {
  m_Genders = genders;
}

double QuantBRLMM::getMaxScore() { return(m_BRLMMMaxScore); };
void QuantBRLMM::setMaxScore(float v){
  m_BRLMMMaxScore=v;

  // keep our options in sync so that the state is
  // correctly reported
  setOptValue("MS", ToStr(m_BRLMMMaxScore));
};


double QuantBRLMM::getMinThresh() { return -DBL_MAX;}
double QuantBRLMM::getMaxThresh() { return m_BRLMMMaxScore;}

void QuantBRLMM::setDmCallsOut(const std::string& fileName, std::vector<std::string> &col_names) {
  //    m_DmCallsFileName = fileName;
  //    if (m_DmCallsOut.is_open())
  //      m_DmCallsOut.close();
  //    Util::mustOpenToWrite(m_DmCallsOut, m_DmCallsFileName.c_str());
  //    m_DmCallsOut << "#Calls: -1=NN, 0=AA, 1=AB, 2=BB" << endl;
  //    m_DmCallsOut << "probeset_id";
  //    for (unsigned int i = 0; i < headers.size(); i++) {
  //      std::string header = Util::fileRoot(headers[i]);
  //      m_DmCallsOut << "\t" << header;
  //    }
  //    m_DmCallsOut << endl;
  m_DmCallsTsv.setFilename(fileName);
  m_DmCallsTsv.setFormat(affx::TsvReport::FMT_TSV);
  m_DmCallsTsv.addHeaderComment("Calls: -1=NN, 0=AA, 1=AB, 2=BB");
  m_DmCallsTsv.defineStringColumn(0,0,"probeset_id",TSVREPORT_PROBESET_STRLEN);
  m_DmCallsTsv.defineColumns(Fs::basename(col_names),affx::FILE5_DTYPE_INT,0);
  //
  m_DmCallsTsv.writeTsv_v1();
}

bool QuantBRLMM::prepare(const IntensityMart &iMart) {
  bool result = true;
  for (int i=0; i<m_Reporters.size(); i++)
    result = result && m_Reporters[i]->prepare(*this, iMart);
  return result;
}

bool QuantBRLMM::finish() {
  bool result = true;
  for (int i=0; i<m_Reporters.size(); i++)
    result = result && m_Reporters[i]->finish(*this);
  return result;
}

/**
 * @brief Do the heavy lifting of estimation.
 */
void QuantBRLMM::computeEstimate() {
  transformData(m_AValues, m_BValues, m_Param);
  ClusterPrior *prior = priorForSnp();
  if (m_HaploidSnps.find(m_ProbesetName) == m_HaploidSnps.end()) {
    ClusterModel model;
    /* Check to see if there is a precomputed model. */
    if (!m_SnpModelMap.empty()) {
      std::map<std::string,ClusterModel>::iterator iter;
      iter = m_SnpModelMap.find(m_ProbesetName + "-2");
      if (iter == m_SnpModelMap.end()) {
        Err::errAbort("Can't find model for SNP: " + m_ProbesetName);
      }
      model = iter->second;
    }
    else {
      model = fitModel(m_AValues, m_BValues, m_Distances,
                       m_InitialCalls, *prior, m_ProbesetName, m_Param);
    }
    /* If we are saving models write this one out. */
    if (m_SnpModelsTsv.is_open()) {
      writeModelOut(model,m_SnpModelsTsv, m_ProbesetName, 2);
    }
    brlmmCall(m_AValues, m_BValues, model, m_Calls, m_Confidences, m_Distances, m_Param);
  }
  else {
    callHaploidSnp(*prior);
  }
}

/**
 * Remove the het cluster from a prior by setting center to be very
 * large negative number and removing covariance with other
 * clusters.
 * @param prior - Prior to remove het cluster from.
 * @param trans - transformation used for clustering.
 */
void QuantBRLMM::removeHetFromPrior(ClusterPrior &prior, enum QuantBRLMM::Transformation trans) {

  /* Adjust het center to very large negative number. */
  if (trans == MvA) {
    prior.centers.element(2) = -1 * FLT_MAX; // -5000 in R code.
    prior.centers.element(3) = 0;
  }
  else if (trans == RvT || trans == CCS || trans == CES) {
    prior.centers.element(2) = 0;
    prior.centers.element(3) = -1 * FLT_MAX; // -5000 in R code
  }
  else {
    Err::errAbort("QuantBRLMM::removeHetFromPrior() - Don't recognize transform type: " + ToStr(trans));
  }

  /* cancel het covariance with other centers. */
  for (unsigned int i = 0; i < prior.covars.Nrows(); i++) {
    prior.covars.element(i,2) = 0;
    prior.covars.element(2,i) = 0;
    prior.covars.element(i,3) = 0;
    prior.covars.element(3,i) = 0;
  }

  /* make hets not vary particularly much, .01 comes directly from R code. */
  prior.covars.element(2,2) = .01;
  prior.covars.element(3,3) = .01;
}

/**
 * This is the fix for the snps on chromosome X that don't appear
 * diploid in men (i.e. aren't in the pseudoautosomal region). The
 * concept is to fit and call the model separately for men and
 * women. Implementation is to set all the male initial calls to NN
 * when fitting for women. For fitting and calling men the female
 * samples are set to NN and further the prior is ajusted such that
 * het calls are impossible to get.
 */
void QuantBRLMM::callHaploidSnp(ClusterPrior &prior) {
  vector<affx::GType> xxCalls, xyCalls;
  vector<double> xxConfidences, xyConfidences;
  vector<vector<double> > xxDistances, xyDistances;
  vector<affx::GType> xxInitCalls(m_InitialCalls), xyInitCalls(m_InitialCalls);
  ClusterModel model;
  std::map<std::string,ClusterModel>::iterator iter;
  /* Sanity check. */
  Err::check(m_Genders.size() == m_InitialCalls.size(),
             "Error: QuantBRLMM::callHaploidSnp() - m_Gender and m_InitialCalls are different sizes.");

  /* Do the female fitting and calls. */
  for (unsigned int i = 0; i < m_InitialCalls.size(); i++) {
    if (m_Genders[i] == affx::Male) {
      xxInitCalls[i] = NN;
    }
  }
  ClusterModel xxModel;
  /* Check to see if there is a precomputed model. */
  if (!m_SnpModelMap.empty()) {
    string name = m_ProbesetName + "-2";
    iter = m_SnpModelMap.find(name);
    if (iter == m_SnpModelMap.end()) {
      Err::errAbort("Can't find model for SNP: " + m_ProbesetName + "-2");
    }
    xxModel = iter->second;
  }
  else {
    xxModel = fitModel(m_AValues, m_BValues, xxDistances,
                 xxInitCalls, prior, m_ProbesetName, m_Param);

  }
  /* If we are saving models write this one out. */
  if (m_SnpModelsTsv.is_open()) {
    writeModelOut(xxModel, m_SnpModelsTsv, m_ProbesetName, 2);
  }
  brlmmCall(m_AValues, m_BValues, xxModel, xxCalls, xxConfidences, xxDistances, m_Param);

  /* Do the male fitting and calls. */
  for (unsigned int i = 0; i < m_InitialCalls.size(); i++) {
      if (m_Genders[i] == affx::Female || m_Genders[i] == affx::UnknownGender) {
      xyInitCalls[i] = NN;
    }
  }
  ClusterPrior xyPrior = prior;
  /* This is where we adjust the prior to remove het cluster. */
  removeHetFromPrior(xyPrior, m_Param.m_Transform);
  /* Then fit as we usually would. */
  ClusterModel xyModel;
  /* Check to see if there is a precomputed model. */
  if (!m_SnpModelMap.empty()) {
    string name = m_ProbesetName + "-1";
    iter = m_SnpModelMap.find(name);
    if (iter == m_SnpModelMap.end()) {
      Err::errAbort("Can't find model for SNP: " + m_ProbesetName + "-1");
    }
    xyModel = iter->second;
  }
  else {
    xyModel = fitModel(m_AValues, m_BValues, xyDistances,
                       xyInitCalls, xyPrior, m_ProbesetName, m_Param);
  }
  /* If we are saving models write this one out. */
  if (m_SnpModelsTsv.is_open()) {
    writeModelOut(xyModel, m_SnpModelsTsv, m_ProbesetName, 1);
  }
  brlmmCall(m_AValues, m_BValues, xyModel, xyCalls, xyConfidences, xyDistances, m_Param);

  m_Calls.resize(xyCalls.size());
  m_Confidences.resize(xyConfidences.size());
  m_Distances.resize(xyDistances.size());
  /* fill in the results appropriate for each gender. */
  for (unsigned int i = 0; i < m_Genders.size(); i++) {
      if (m_Genders[i] == affx::Female || m_Genders[i] == affx::UnknownGender) {
      m_Calls[i] = xxCalls[i];
      m_Confidences[i] = xxConfidences[i];
      m_Distances[i] = xxDistances[i];
    }
    else if (m_Genders[i] == affx::Male) {
      m_Calls[i] = xyCalls[i];
      m_Confidences[i] = xyConfidences[i];
      m_Distances[i] = xyDistances[i];
    }
    else {
      Err::errAbort("QuantBRLMM::callHaploidSnp() - Don't recognize gender: " + ToStr(m_Genders[i]));
    }
  }
}

/**
 * Write out a text representation of the prior to the file specified.
 * @param prior - Prior to be written.
 * @param out - ostream to write to.
 */
void QuantBRLMM::writePriorOut(ClusterPrior &prior, std::ostream &out,
                               const std::string& fieldDelim, const std::string& sep) {
  printColumnVector(prior.centers, &out, sep);
  out << fieldDelim;
  printColumnVector(prior.centerVars, &out, sep);
  out << fieldDelim;
  printMatrix(prior.covars, &out, sep);
  out << endl;
}

void QuantBRLMM::setSnpModelOut(const std::string& fileName,affx::TsvReport::TsvReportFmt_t format)
{
  //
  // printf("### QuantBRLMM::setSnpModelOut('%s',%d)\n",fileName.c_str(),format);
  //
  m_SnpModelsTsv.setFilename(fileName);
  m_SnpModelsTsv.setFormat(format);
  //
  m_SnpModelsTsv.defineStringColumn(0,0,"probeset_id",TSVREPORT_PROBESET_STRLEN);
  m_SnpModelsTsv.defineColumn(0,1,"copy_number",affx::FILE5_DTYPE_INT);
  //
  m_SnpModelsTsv.defineStringColumn(0,2,"center",TSVREPORT_PROBESET_STRLEN);
  m_SnpModelsTsv.defineStringColumn(0,3,"var",TSVREPORT_PROBESET_STRLEN);
  // Add headers here.
  m_SnpModelsTsv.writeTsv_v1();
}

/**
 * Write out a text representation of the model to the file specified.
 * @param model - Model to be written.
 * @param out - ostream to write to.
 */
void QuantBRLMM::writeModelOut(ClusterModel& model,affx::TsvReport& tsv, const std::string& name,int copies)
{
  tsv.set_string(0,0,name); // "probeset_id"
  tsv.set_i(0,1,copies);    // "copy_number"
  //
  tsv.set_string(0,2,stringMatrix(model.clusterMeans,",")); // "center"
  tsv.set_string(0,3,stringMatrix(model.clusterVars, ",")); // "var"
  //
  tsv.writeLevel(0);
}

/**
 * Write out a text representation of the prior to the file specified.
 * @param prior
 * @param fileName
 */

// Why is this singular
void QuantBRLMM::writePrior(ClusterPrior &prior, const std::string& fileName) {
  std::string fieldDelim = "\t";
  std::string sep=",";
  ofstream out;
  Fs::mustOpenToWrite(out, fileName);
  out.setf(ios::fixed, ios::floatfield);
  out.precision(10);
  out << "center" << fieldDelim << "var" << fieldDelim << "center.var" << endl;
  QuantBRLMM::writePriorOut(prior, out, fieldDelim, sep);
}

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
ClusterPrior QuantBRLMM::clusterPriorFromStrings(const std::string& fileName, std::string &id, std::string &center,
                                                 std::string &var, std::string &covar) {
  ClusterPrior prior;
  vector<double> vec;
  vector<string> words;

  prior.id = id;

  /* Read in the centers. */
  Util::chopString(center, ',', words);
  vec.clear();
  if (words.size() != 6)
    Err::errAbort("Expecting 6 entires in 'center' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for (unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.centers = columnVectorFromArray(&vec[0], (int)vec.size());

  /* Read in the variances. */
  Util::chopString(var, ',', words);
  vec.clear();
  if (words.size() != 9)
    Err::errAbort("Expecting 9 entires in 'var' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for (unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.centerVars = columnVectorFromArray(&vec[0], (int)vec.size());

  /* Do the covariacne matrix. */
  Util::chopString(covar, ',', words);
  vec.clear();
  if (words.size() != 36)
    Err::errAbort("Expecting 36 entires in 'centers.var' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for (unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.covars = matrixFromArray(&vec[0], (int)vec.size(), 6, 6);
  return prior;

}

/**
 * Covert string representations into ClusterModel
 *
 * @param id - identifying string.
 * @param center - centers.
 * @param var - variance of center variances.
 *
 * @return - ClusterModel
 */
ClusterModel QuantBRLMM::clusterModelFromStrings(const std::string& fileName, std::string &id,
                                                 std::string &center, std::string &var) {
  ClusterModel model;
  vector<double> vec(9);
  vector<string> words;
  /* Read in the centers. */
  Util::chopString(center, ',', words);
  vec.clear();
  if (words.size() != 6)
    Err::errAbort("Expecting 6 entires in 'center' entry in file: " + ToStr(fileName) +
                  " got: " + ToStr(words.size()) + " for probest: " + id);
  for (unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  model.clusterMeans = matrixFromArray(&vec[0], (int)vec.size(), 6, 1);

  /* Read in the variances. */
  Util::chopString(var, ',', words);
  vec.clear();
  if (words.size() != 9)
    Err::errAbort("Expecting 9 entires in 'var' entry in file: " + ToStr(fileName) +
                  " got: " + ToStr(words.size()) + " for probest: " + id);
  for (unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  model.clusterVars = matrixFromArray(&vec[0], (int)vec.size(), 9, 1);
  return model;
}


/**
 * Read the prior from the test file specified.
 * @param fileName - Name of file to read from.
 * @return - The prior read in from file.
 */
ClusterPrior QuantBRLMM::loadPrior(const std::string& fileName) {
  affx::TsvFile tsv;
  string center, centerVar, var, id="unknown";
  ClusterPrior prior;
  vector<double> vec;
  tsv.open(fileName);
  //tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center.var", &centerVar, affx::TSV_BIND_REQUIRED);
  if (tsv.nextLevel(0) != affx::TSV_OK) {
    Err::errAbort("Didnt' get an entry in file: " + ToStr(fileName));
  }
  prior = QuantBRLMM::clusterPriorFromStrings(fileName, id, center, var, centerVar);
  tsv.close();
  return prior;
}

/**
 * Generate a string representation of a snp model.
 * @return string - text description of snp clusters.
 */
std::string QuantBRLMM::getModelString() {
  int i = 0, j = 0;
  std::string s = m_ProbesetName;
  s += "\t";
  for (i = 0; i < m_CurrentPrior->centers.Nrows() - 1; i++) {
    s += ToStr(m_CurrentPrior->centers.element(i));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->centers.element(i));
  s += "\t";
  for (i = 0; i < m_CurrentPrior->centerVars.Nrows() - 1; i++) {
    s += ToStr(m_CurrentPrior->centerVars.element(i));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->centerVars.element(i));
  s += "\t";
  for (i = 0; i < m_CurrentPrior->covars.Nrows()-1; i++) {
    for (j = 0; j < m_CurrentPrior->covars.Ncols(); j++) {
      s += ToStr(m_CurrentPrior->covars.element(i,j));
      s += ",";
    }
  }
  for (j = 0; j < m_CurrentPrior->covars.Ncols() - 1; j++) {
    s += ToStr(m_CurrentPrior->covars.element(i,j));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->covars.element(i,j));
  return s;
}

/**
 * Load in a number of snp specific models.
 * @param snpModels - map of snp probeset names to models.
 * @param fileName - filename to read models from.
 */
void QuantBRLMM::loadSnpModels(std::map<std::string, ClusterModel> &snpModels,const std::string& fileName) {
  affx::TsvFile tsv;
  string center, var, id, copies;
  ClusterModel model;
  snpModels.clear();
  tsv.open(fileName);

  // we want to allow either "probeset_id" or "id"...
  int probeset_id_cidx;
  probeset_id_cidx=tsv.cname2cidx(0,"probeset_id","id");
  if (probeset_id_cidx<0) {
    Err::errAbort("QuantBRLMM::loadSnpModels: no 'probeset_id' or 'id' column");
  }
  int copy_number_cidx;
  copy_number_cidx=tsv.cname2cidx(0,"copy_number","copies");
  if (copy_number_cidx<0) {
    Err::errAbort("QuantBRLMM::loadSnpModels: no 'copy_number' or 'copies' column");
  }
  //
  tsv.bind(0,probeset_id_cidx, &id    , affx::TSV_BIND_REQUIRED);
  tsv.bind(0,copy_number_cidx, &copies, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  //
  while (tsv.nextLevel(0) == affx::TSV_OK) {
    model = clusterModelFromStrings(fileName, id, center, var);
    if (snpModels.find(id) != snpModels.end()) {
      Err::errAbort("QuantBRLMM::loadSnpModels() - Id: " + ToStr(id) + " seen multiple times.");
    }
    snpModels[id + "-" + copies] = model;
  }
  tsv.close();
}

/**
 * Load in a number of snp specific priors.
 * @param snpPriors - map of snp probeset names to priors.
 * @param fileName - filename to read priors from.
 */
void QuantBRLMM::loadSnpPriors(std::map<std::string,ClusterPrior> &snpPriors,
                               const std::string& fileName) {
  affx::TsvFile tsv;
  string center, centerVar, var, id="unknown";
  ClusterPrior prior;
  snpPriors.clear();
  tsv.open(fileName);
  tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center.var", &centerVar, affx::TSV_BIND_REQUIRED);
  while (tsv.nextLevel(0) == affx::TSV_OK) {
    prior = clusterPriorFromStrings(fileName, id, center, var, centerVar);
    if (snpPriors.find(prior.id) != snpPriors.end()) {
      Err::errAbort("QuantBRLMM::loadSnpPriors() - Id: " + ToStr(prior.id) + " seen multiple times.");
    }
    snpPriors[prior.id] = prior;
  }
  tsv.close();
}


/**
 * Compute the prior from a collection of snps and fill in the prior.
 *
 * @param snpStats - Stats for a relatively large number of snps
 * @param prior - Fill in the prior.
 */
void QuantBRLMM::computeClusterPrior(const std::vector<ClusterStats> &snpStats, ClusterPrior &prior) {
  vector<vector<double> > means(6);
  vector<vector<double> > vars(9);
  prior.snpCount = snpStats.size();
  /* Save up vectors of the values for calculating stats. */
  for (unsigned int snpIx = 0; snpIx < snpStats.size(); snpIx++) {
    const ClusterStats &stats = snpStats[snpIx];
    for (unsigned int meanIx = 0; meanIx < stats.means.size(); meanIx++) {
      means[meanIx].push_back(stats.means[meanIx]);
    }
    for (unsigned int varIx = 0; varIx < stats.vars.size(); varIx++) {
      vars[varIx].push_back(stats.vars[varIx]);
    }
  }
  prior.centers.ReSize(6);
  prior.centerVars.ReSize(9);
  prior.covars.ReSize(6,6);
  for (unsigned int meanIx = 0; meanIx < means.size(); meanIx++) {
    prior.centers.element(meanIx) = Average<vector<double>::iterator, double>()(means[meanIx].begin(),
                                                                               means[meanIx].end());
    for (unsigned int nIx = 0; nIx < means.size(); nIx++) {
      prior.covars.element(meanIx, nIx) = UnbiasedCovariance<vector<double>::iterator, double>()(means[meanIx].begin(),
                                                                                                means[meanIx].end(),
                                                                                                means[nIx].begin());
    }
  }
  for (unsigned int varIx = 0; varIx < vars.size(); varIx++)
    prior.centerVars.element(varIx) = Average<vector<double>::iterator, double>()(vars[varIx].begin(), vars[varIx].end());
}

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
ClusterPrior QuantBRLMM::makePrior(std::vector<string>& probeSetNames,
                                   ChipLayout& layout,
                                   IntensityMart& iMart,
                                   std::vector<ChipStream *>& iTrans,
                                   PmAdjuster &pmAdjust)
{
  ClusterPrior prior;
  vector<ClusterStats> snpStats;
  int minCalls = m_Param.m_MinCluster;
  Err::check(minCalls >= 2, "QuantBRLMM::makePrior() - Threshold must have at least two calls per genotype.");
  try {
    for (unsigned int psIx = 0; psIx < probeSetNames.size(); psIx++) {
      ProbeListPacked pList = layout.getProbeListByName(probeSetNames[psIx]);
      ///@todo deal with missing probes from kill list
      if (pList.isNull()) {
         Err::errAbort("QuantBRLMM::makePrior(): Cannot deal with missing probes from kill list");
      }
      ProbeSet *gtPs = ProbeListFactory::asProbeSet(pList);
      // Have to skip probesets that don't have 2 or 4 groups as it
      // appears that some are mislabelled in the cdf file.
      if (gtPs->psType != ProbeSet::GenoType || (gtPs->numGroups != 2 && gtPs->numGroups != 4))
        continue;
      m_AValues.clear();
      m_BValues.clear();
      /* Load up the postulated genotypes for this snp. */
      //map<string, vector<char> >::iterator mapIter = m_KnownGenoTypes.find(gtPs->name);
      KnownGenoTypes_iter_t mapIter = m_KnownGenoTypes.find(gtPs->name);
      if (mapIter == m_KnownGenoTypes.end()) {
        Err::errAbort("QuantBRLMM::makePrior() - Can't find genotype data for SNP: " + ToStr(gtPs->name));
      }
      vector<GType> &gType = mapIter->second;
      m_InitialCalls.clear();
      for (unsigned int i = 0; i < gType.size(); i++) {
        m_InitialCalls.push_back(gType[i]);
      }
      fillInAlleleProbeSets(*gtPs, m_Aallele, m_Ballele);
      if (!summarizeAllele(&m_Aallele, m_AValues, iMart, iTrans, pmAdjust, m_QuantMethod, m_Param.m_LowPrecision))
        Err::errAbort("Couldn't calculate A summary values for SNP: " + ToStr(gtPs->name));
      if (!summarizeAllele(&m_Ballele, m_BValues, iMart, iTrans, pmAdjust, m_QuantMethod, m_Param.m_LowPrecision))
        Err::errAbort("Couldn't calculate B summary values for SNP: " + ToStr(gtPs->name));
      transformData(m_AValues, m_BValues, m_Param);
      ClusterStats stats = computeClusterStats(m_AValues, m_BValues, m_InitialCalls,
                                               m_Param.m_Transform, m_Param.m_K);
      /* Only use these stats if there at least two data points per snp. */
      if (stats.counts.size() == 3 && stats.counts[0] >= minCalls &&
         stats.counts[1] >= minCalls && stats.counts[2] >= minCalls) {
        snpStats.push_back(stats);
      }
      delete gtPs;
    }
    computeClusterPrior(snpStats, prior);
  }
  catch(...) {
    Err::errAbort("exception while calculating prior.");
  }
  Verbose::out(1, ToStr(snpStats.size()) + " of " + ToStr(probeSetNames.size()) + " (" +
               ToStr((float)snpStats.size()/probeSetNames.size()*100.0) +
               "%) of SNPs had at least 2 observations per genotype.");
  return prior;
}

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
bool QuantBRLMM::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                       std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {
  const ProbeSet *gtPs = NULL; // Genotyping probeset
  blankSelf(); // Make sure we clear out the past analysis before starting this one.
  bool success = true;
  /* Sanity checks about probesets. */
  if (psGroup.probeSets.empty())
    Err::errAbort("Zero probesets in ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");
  /* Remember this probeset. */
  gtPs = psGroup.probeSets[0];
  m_GtProbeSet = gtPs;
  m_ProbesetName = gtPs->name;
  if (!canSetUpProbeSet(gtPs)) {return false;}
  if (psGroup.probeSets.size() > 1)
    Err::errAbort("Can't have multiple probesets in a genotyping ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");

  /* Load up the postulated genotypes for this snp. */
  KnownGenoTypes_t::iterator mapIter = m_KnownGenoTypes.find(gtPs->name);
  if (mapIter == m_KnownGenoTypes.end()) {
    blankSelf();
    Verbose::out(2, "QuantBRLMM::setUp() - Can't find genotype data for SNP: (probeset: name=" + ToStr(gtPs->name) + ", type=" + ProbeSet::typeEnumToString(gtPs->psType) + ")");
    return false;
  }
  if (m_QuantMethod == NULL) {
    QuantMethodFactory factory(QuantMethodFactory::Expression);
    ChipLayout layout;
    string spec = "plier.optmethod=1";
    m_QuantMethod = factory.quantExprMethodForString(spec, layout, QuantMethodFactory::Expression);
  }

  vector<affx::GType> &gType = mapIter->second;
  m_InitialCalls.clear();
  for (unsigned int i = 0; i < gType.size(); i++) {
     m_InitialCalls.push_back(gType[i]);
   }
  fillInAlleleProbeSets(*gtPs, m_Aallele, m_Ballele);

  /* Get summaries for each allele. */
  success &= summarizeAllele(&m_Aallele, m_AValues, iMart, iTrans, pmAdjust,
                             m_QuantMethod, m_Param.m_LowPrecision, true, m_Reporters);
  if (gtPs->psType == ProbeSet::Copynumber) {
    blankSelf();
    return false;
  }
  success &= summarizeAllele(&m_Ballele, m_BValues, iMart, iTrans, pmAdjust,
                             m_QuantMethod, m_Param.m_LowPrecision, true, m_Reporters);
  //
  if (m_DmCallsTsv.is_open()) {
    m_DmCallsTsv.set_string(0,0,m_ProbesetName);
    for (unsigned int i = 0; i < m_InitialCalls.size(); i++) {
      m_DmCallsTsv.set_i(0,1+i,affx::GType_to_int(m_InitialCalls[i]));
    }
    m_DmCallsTsv.writeLevel(0);
  }

  if (!success)
    blankSelf();
  return success;
}

ClusterPrior *QuantBRLMM::priorForSnp() {
  if (m_SnpPriorMap.empty()) {
    m_CurrentPrior = m_Prior;
    return m_Prior;
  }
  else {
    std::map<std::string,ClusterPrior>::iterator mapIter = m_SnpPriorMap.find(m_ProbesetName);
    if (mapIter == m_SnpPriorMap.end()) {
      m_CurrentPrior = m_Prior;
      return m_Prior;
      //      Err::errAbort("QuantBRLMM::priorForSnp() - Can't find prior for SNP: " + m_ProbesetName);
    }
    m_CurrentPrior = &mapIter->second;
    return &mapIter->second;
  }
}

/** Conversion utility for transforming vectors to column vectors. */
ColumnVector QuantBRLMM::columnVectorFromVector(std::vector<double> &v) {
  ColumnVector x;
  x.ReSize(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    x.element(i) = v[i];
  return x;
}

/** Conversion utility for transforming arrays to column vectors. */
ColumnVector QuantBRLMM::columnVectorFromArray(double *array, int count) {
  ColumnVector x;
  x.ReSize(count);
  for (int i = 0; i < count; i++)
    x.element(i) = array[i];
  return x;
}

/** Conversion utility for transforming arrays to matrices. */
Matrix QuantBRLMM::matrixFromArray(double *array, int count, int nRow, int nCol) {
  Matrix m(nRow, nCol);
  assert(nRow * nCol == count);
  for (int i = 0; i < count; i++) {
    int col = i / nRow;
    int row = i - col * nRow;
    m.element(row, col) = array[i];
  }
  return m;
}


/** Make banded count matrix. */
Matrix QuantBRLMM::makeCountMat(std::vector<int> &counts) {
  Matrix m(6,6);
  m = 0;
  assert(counts.size() == 3);
  m.element(0,0) = counts[0];
  m.element(0,1) = counts[0];
  m.element(1,0) = counts[0];
  m.element(1,1) = counts[0];

  m.element(2,2) = counts[1];
  m.element(2,3) = counts[1];
  m.element(3,2) = counts[1];
  m.element(3,3) = counts[1];

  m.element(4,4) = counts[2];
  m.element(4,5) = counts[2];
  m.element(5,4) = counts[2];
  m.element(5,5) = counts[2];
  return m;
}

/** Reformat variance/covariance vector as a banded matrix. */
Matrix QuantBRLMM::formatAsBigVarMat(ColumnVector v) {
  Matrix m(6, 6);
  m = 0;
  int nRows = v.Nrows();
  assert(nRows == 9);
  // Diagonal are m and a variances.
  m.element(0,0) = v.element(0);
  m.element(1,1) = v.element(1);
  m.element(2,2) = v.element(3);
  m.element(3,3) = v.element(4);
  m.element(4,4) = v.element(6);
  m.element(5,5) = v.element(7);
  // Fill in the covariances appropriately
  m.element(0,1) = m.element(1,0) = v.element(2);
  m.element(2,3) = m.element(3,2) = v.element(5);
  m.element(4,5) = m.element(5,4) = v.element(8);
  /// @todo what about other points like 1,3 or 1,5?
  return m;
}

/** Bayesian update of the mean vector, this is where the voodoo
    happens. Have to query Simon more about the motivation. */
Matrix QuantBRLMM::bayesianMeanUpdate(ColumnVector &newVar, std::vector<int> &counts, ColumnVector &means,
                                             Matrix &priorCenterVar, ColumnVector &priorCenter) {
  Matrix sigmaInv = QuantBRLMM::formatAsBigVarMat(newVar);
  sigmaInv = sigmaInv.i();
  Matrix countMatrix = QuantBRLMM::makeCountMat(counts);
  Matrix centerVarInv = priorCenterVar.i();
  Matrix elementMult = elementMultiplication(countMatrix, sigmaInv);
  Matrix t1 = centerVarInv + elementMult;
  t1 = t1.i();
  Matrix t2 = centerVarInv * priorCenter;
  Matrix temp = elementMult * means;
  t2 = t2 + temp;
  Matrix newMean = t1 * t2;
  return newMean;
}

/// @todo should degFree be weight or counts rather than degrees of freedom?
ColumnVector QuantBRLMM::shrinkVariance(ClusterStats &empirical, ClusterPrior &prior, int degFree) {
  ColumnVector newVar;
  newVar.ReSize(empirical.vars.size());
  /// @todo improve handling of NA rather than setting to zero
  for (int i = 0; i < empirical.vars.size(); i++) {
    int counts = empirical.counts[i/3];
    /* empirical variance gets set to zero if there weren't enough
       observations for variance estimations (i.e. less than 2). If we
       don't have enough observations to compute a variance set the counts
       to zero so that we just use the prior estimate of variance. */
    if (counts < 2) {
      counts = 0;
    }
    int totalCount = counts + degFree;
    double result = ((counts * empirical.vars[i]) +
                    (degFree * prior.centerVars.element(i)))/totalCount;
    newVar.element(i) = result;
  }
  return newVar;
}

/** Wrapper function for getting bayesian adjusted parameters given data. */
ClusterModel QuantBRLMM::fitModel(std::vector<double> &aValues,
                                  std::vector<double> &bValues,
                                  std::vector<std::vector<double> > &genoDist,
                                  std::vector<affx::GType> &genoTypes,
                                  ClusterPrior &prior, const std::string &name,
                                  BrlmmParam &param) {

  enum Transformation transform = param.m_Transform;
  int iterations = param.m_Iterations;
  double iterationThreshold = param.m_IterationThreshold;
  double K = param.m_K;

  ClusterStats stats = QuantBRLMM::computeClusterStats(aValues, bValues, genoTypes, transform, K );

  ClusterModel model;
  ColumnVector newVar = shrinkVariance(stats, prior, param.m_PriorPseudoCount);
  model.clusterVars = newVar;
  ColumnVector means = columnVectorFromVector(stats.means);
  model.clusterMeans = QuantBRLMM::bayesianMeanUpdate(newVar, stats.counts, means, prior.covars, prior.centers);
  /* If iterations are greater than zero then we are going to recall
     using current model and relearn a new model using those calls
     instead of original calls. */
  while (iterations > 0) {
    iterations--;
    vector<affx::GType> newCalls;  // brlmm calls for current iteration.
    vector<double> newConfidences; // brlmm confidences for current iteration.
    /* Make initial call using current model. */
    brlmmCall(aValues, bValues, model, newCalls, newConfidences, genoDist, param);

    /* Set calls that are less than our threshold to NN */
    assert(newCalls.size() == newConfidences.size());
    for (unsigned int i = 0; i < newCalls.size(); i++) {
      if (newConfidences[i] > iterationThreshold) {
        newCalls[i] = NN;
      }
    }
    /* Compute a new model with updated calls. */
    ClusterStats stats = computeClusterStats(aValues, bValues, newCalls, transform, K);
    ColumnVector newVar = shrinkVariance(stats, prior, param.m_PriorPseudoCount);
    model.clusterVars = newVar;
    ColumnVector means = columnVectorFromVector(stats.means);
    model.clusterMeans = QuantBRLMM::bayesianMeanUpdate(newVar, stats.counts, means, prior.covars, prior.centers);
  }
  return model;
}

void QuantBRLMM::valueForTransformation(double aVal, double bVal, enum QuantBRLMM::Transformation trans,
                                         double K, double &xVal, double &yVal) {
  double c = 0.001;
  double denom = 0;
  double aV = 0.0;
  double bV = 0.0;
  double r = 0.0;
  switch (trans) {
  case MvA:
    if (aVal == 0) aVal = c;
    if (bVal == 0) bVal = c;
    xVal = (log2(bVal) + log2(aVal)) / 2;
    yVal = log2(bVal) - log2(aVal);
    break;
  case RvT:
    // Check to make sure we're not got to take log of 0
    denom = bVal;
    if (bVal == 0) bVal = c;
    xVal = atan(aVal/ denom);
    yVal = log( sqrt((aVal * aVal) + (bVal * bVal) + c) );
    break;
  case  CES:
    // Check to make sure we're not got to take log of 0
    denom = bVal + aVal;
    if (denom == 0) denom = c;
    xVal = sinh( K * (aVal - bVal)/denom) / sinh(K);
    yVal = log2(denom);
    break;
  case CCS:
    denom = bVal + aVal;
    if (denom == 0) denom = c;
    aV = aVal;
    if (aV == 0) aV = c;
    bV = bVal;
    if (bV == 0) bV = c;
    r = K * (aV - bV)/denom;
    r = log(r + sqrt(r*r+1));
    xVal = r / log(K + sqrt(K*K+1));
    //r = K * (aV - bV)/denom;
    // r = log(r + sqrt(r*r+1));
    //    xVal = (log(r + sqrt(r*r+1))) / log(K + sqrt(K*K+1));
    // Note that we are using natural logarithm here (log_e), not log2
    // xVal = asinh(K * tanh( .5 * (log(aV) - log(bV)))) / asinh(K);
    yVal = log2(denom);
    break;
  default:
    Err::errAbort("QuantBRLMM::fillInTransformedData() - Don't recognize transformation type: " + ToStr(trans));
  }
}

/**
 * Transform the data in aValues and bValues using the method specified by trans.
 *
 * @param aValues - Summary values for A allele.
 * @param bValues - Summary values for B allele.
 * @param genoTypes - Genotypes for each sample.
 * @param Transformation - Type of transformation to use.
 * @param xDim - two dimensional matrix for x dimension with rows
 * (index 0) being genotype (0,1,2) and columns being values observed.
 * @param yDim - two dimensional matrix for y dimension with rows
 * (index 0) being genotype (0,1,2) and columns being values observed.
 */
void QuantBRLMM::fillInTransformedData(std::vector<double> &aValues,
                                       std::vector<double> &bValues,
                                       std::vector<affx::GType> &genoTypes,
                                       enum QuantBRLMM::Transformation trans,
                                       double K,
                                       std::vector<std::vector<double> > &xDim,
                                       std::vector<std::vector<double> > &yDim) {
  /* Sanity check. */
  assert(aValues.size() == bValues.size());
  assert(aValues.size() == genoTypes.size());

  /* Clean out any previous results. */
  xDim.clear();
  yDim.clear();

  xDim.resize(3); // Number of genotypes 0, 1, 2
  yDim.resize(3); // Number of genotypes 0, 1, 2
  //  double xVal = 0, yVal = 0;
  for (unsigned int i = 0; i < aValues.size(); i++) {
    if (genoTypes[i] == NN)
      continue;
    //    valueForTransformation(aValues[i], bValues[i], trans, K, xVal, yVal);
    xDim[genoTypes[i]].push_back( aValues[i] );
    yDim[genoTypes[i]].push_back( bValues[i] );
  }
}

/** Count up some elementary statistics about a particular snp over
    multiple experiments. */
ClusterStats QuantBRLMM::computeClusterStats(std::vector<double> &aValues,
                                             std::vector<double> &bValues,
                                             std::vector<affx::GType> &genoTypes,
                                             enum QuantBRLMM::Transformation trans,
                                             double K) {
  // Sanity checks
  assert(aValues.size() > 0);
  Err::check(aValues.size() == bValues.size(), "computClusterStats() - Expression should be same size.");
  Err::check(bValues.size() == genoTypes.size(), "computClusterStats() - Expression should be same size as genotyping.");
  // Indexes are AA = 0, AB=1, BB=2 dType = "difference" bType =
  // "intensity" Normally these would be M & A (like MA plots), but
  // A gets confusing with genotype 'A'
  vector<vector <double> > dTypeValues(3);
  vector<vector <double> > iTypeValues(3);
  //  bool allGood = true;
  ClusterStats stats;
  // Loop through and group into genotype specific data.
  fillInTransformedData(aValues, bValues, genoTypes, trans, K, iTypeValues, dTypeValues);

  // Compute the means of the clusters for each genotype. Note that if
  // we don't have enough samples to compute the variance it is set to
  // zero as a placeholder, but it is necessary to check the counts
  // later to make sure we're not using an invalid number.
  for (unsigned int i = 0; i < 3; i++) {
    stats.counts[i] = iTypeValues[i].size();
    if (stats.counts[i] > 0) {
      stats.mMeans[i] = Average<vector<double>::iterator, double>()(iTypeValues[i].begin(), iTypeValues[i].end());
      stats.aMeans[i] = Average<vector<double>::iterator, double>()(dTypeValues[i].begin(), dTypeValues[i].end());
      if (stats.counts[i] > 1) {
        stats.mVars[i] = UnbiasedVariance<vector<double>::iterator, double>()(iTypeValues[i].begin(), iTypeValues[i].end());
        stats.aVars[i] = UnbiasedVariance<vector<double>::iterator, double>()(dTypeValues[i].begin(), dTypeValues[i].end());
        stats.covars[i] = UnbiasedCovariance<vector<double>::iterator, double>()(dTypeValues[i].begin(), dTypeValues[i].end(), iTypeValues[i].begin());
      }
      else {
        stats.mVars[i] = 0;
        stats.aVars[i] = 0;
        stats.covars[i] = 0;
      }
    }
    else {
      stats.mMeans[i] = 0;
      stats.mVars[i] = 0;
      stats.aMeans[i] = 0;
      stats.aVars[i] = 0;
      stats.covars[i] = 0;
    }
    stats.means.push_back(stats.mMeans[i]);
    stats.means.push_back(stats.aMeans[i]);
    stats.vars.push_back(stats.mVars[i]);
    stats.vars.push_back(stats.aVars[i]);
    stats.vars.push_back(stats.covars[i]);
  }
  return stats;
}

/* Make the inverse matrix from the variance and covariances. */
Matrix QuantBRLMM::makeVarInv(Matrix &clusterVars, int index) {
  Matrix var(2,2);
  index = 3 * index;
  var.element(0,0) = clusterVars.element(index,0);
  var.element(1,1) = clusterVars.element(index+1,0);
  var.element(0,1) = clusterVars.element(index+2,0);
  var.element(1,0) = clusterVars.element(index+2,0);
  return var.i();
}

/* Take the model derived in fitModel and use it to do predictions. */
void QuantBRLMM::brlmmCall(const std::vector<double> &aValues,
                           const std::vector<double> &bValues,
                           ClusterModel &model,
                           std::vector<affx::GType> &calls,
                           std::vector<double> &confidences,
                           std::vector<std::vector<double> > &genoDist,
                           BrlmmParam &param)
{
  //
  assert(aValues.size() > 0);
  assert(aValues.size() == bValues.size());
  // JHG - debugging
  // printVec("aValues",aValues);
  // printVec("aValues",bValues);
  //
  double hetMult = param.m_HetMultiplier;
  Matrix obs(aValues.size(),2); // Our observed means in both M and A dimensions.
  Matrix distance(aValues.size(),3); // Our distances to clusters.
  calls.clear();
  confidences.clear();

  /* Load up the observed means into matrix. */
  for (unsigned int i = 0; i < aValues.size(); i++) {
    obs.element(i, 0) = aValues[i];
    obs.element(i, 1) = bValues[i];
  }

  /* Calculate the distances for each genotype. */
  for (int colIx = 0; colIx < distance.Ncols(); colIx++) {
    Matrix varInv = QuantBRLMM::makeVarInv(model.clusterVars, colIx);
    // I'm thinking these are the covariances for the mahalanobis distance
    double t1 = varInv.element(0,0);
    double t2 = varInv.element(0,1) + varInv.element(1,0); /// @todo why is this additive?
    double t3 = varInv.element(1,1);

    /* Calculate the distance for different samples for this genotype. */
    for (int rowIx = 0; rowIx < distance.Nrows(); rowIx++) {
      double iDist =  obs.element(rowIx,0) - model.clusterMeans.element(colIx * 2, 0);
      double dDist =  obs.element(rowIx,1) - model.clusterMeans.element(colIx * 2 + 1, 0);
      double dist = iDist*iDist*t1 + iDist*dDist*t2 + dDist*dDist*t3;
      /* If we are calculating distance for AB genotype (column 1) then multiply by our
         het multiplier as a way to balance het/hom accurancy. */
      if (colIx == 1) {
        dist = dist * hetMult;
      }
      distance.element(rowIx, colIx) = dist;
    }
  }
  genoDist.resize(distance.Nrows());

  /* for each sample predict based on which cluster is closest. */
  for (int rowIx = 0; rowIx < distance.Nrows(); rowIx++) {
    double minDistance = FLT_MAX;
    int minDistIx = -1;
    double secondMin = FLT_MAX;
    for (int colIx = 0; colIx < distance.Ncols(); colIx++) {
      double dist = distance.element(rowIx, colIx);
      genoDist[rowIx].push_back(dist);
      // JHG - debug
      // printf("dist=%.10f\n",dist);
      if (dist < minDistance) {
        secondMin = minDistance;
        minDistance = dist;
        minDistIx = colIx;
      }
      else if (dist < secondMin) {
        secondMin = dist;
      }
    }
    calls.push_back((affx::GType)minDistIx);
    confidences.push_back(minDistance/secondMin);
  }
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> QuantBRLMM::getDefaultDocOptions() {
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt hetMult = {"het-mult", SelfDoc::Opt::Double, "1.0", "1.0", "0", "1.0",
                          "Number to balance het calls with to balance performance on het/hom calls."};
  opts.push_back(hetMult);
  SelfDoc::Opt iterations = {"iterations", SelfDoc::Opt::Integer, "0", "0", "0", "NA",
                             "Number of times to iterate BRLMM classifier, feeding in new calls from previous iteration."};
  opts.push_back(iterations);
  SelfDoc::Opt iterThresh = {"iter-thresh", SelfDoc::Opt::Double, "0.3", "0.3", "0", "1",
                             "Maximum confidence score to use when doing iterations internally [0,1]."};
  opts.push_back(iterThresh);
  SelfDoc::Opt k = {"K", SelfDoc::Opt::Double, "4.0", "4.0", "0", "NA",
                    "Scale parameter used used in CCS and CES transformations."};
  opts.push_back(k);
  SelfDoc::Opt transform = {"transform", SelfDoc::Opt::String, "CCS", "CCS", "NA", "NA",
                            "Transformation of initial data are we feeding into the classifier? {'CCS', 'CES', 'MvA','RvT'}"};
  opts.push_back(transform);
  SelfDoc::Opt priorWeight = {"prior-weight", SelfDoc::Opt::Integer, "40", "40", "0", "NA",
                              "Psuedocount weight should the prior have? Also known as 'degrees of freedom' in R code."};
  opts.push_back(priorWeight);
  SelfDoc::Opt priorMinCall = {"prior-mincall", SelfDoc::Opt::Integer, "2", "2", "2", "NA",
                               "Minimum number of genotypes per cluster for inclusion in prior estimation, must be >= 2."};
  opts.push_back(priorMinCall);
  SelfDoc::Opt lowprecision = {"lowprecision", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                               "R prototype uses summary values rounded to first decimal place. Use this flag to be simulate behavior."};
  opts.push_back(lowprecision);
 SelfDoc::Opt maxscore = {"MS",SelfDoc::Opt::Float, "0.5","0.5", "0", "2",
				"Threshold for making no-calls"};
	opts.push_back(maxscore);
 return opts;
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc QuantBRLMM::explainSelf() {
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}

/**
 * Convert a (text insensitive) text repesentation of a transformation into
 * the enumration.
 * @param s - string specifying transformation.
 * @return - Enumeration version of transformation.
 */
enum QuantBRLMM::Transformation QuantBRLMM::transformationForString(const std::string& str) {
  std::string lowerS = Util::downcaseString(str);

  if (lowerS == "mva")
    return MvA;
  else if (lowerS == "rvt")
    return RvT;
  else if (lowerS == "ssf" || lowerS == "ces")
    return CES;
  else if (lowerS == "assf" || lowerS == "ccs")
    return CCS;
  else
    Err::errAbort("QuantBRLMM::transformationForString() - Don't recognize transformation: " + str);
  return MvA;
}

/**
 * Convert a transformation enum to string representation
 * @param t - tranformation enum value.
 * @return - text description of enumeration.
 */
std::string QuantBRLMM::stringForTransformation(enum QuantBRLMM::Transformation t) {
  switch(t) {
  case MvA :
    return "MvA";
    break;
  case RvT :
    return "RvT";
    break;
  case CES :
    return "CES";
    break;
  case CCS:
    return "CCS";
    break;
  default:
    Err::errAbort("QuantBRLMM::stringForTransformation() - Don't recognize type: " + ToStr(t));
  }
  Err::errAbort("QuantBRLMM::stringForTransformation() - Should never reach this point.");
  return "";
}


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
SelfCreate *QuantBRLMM::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  double hetMult = 1.0, K = 4;
  double iterationThreshold = 0.3;
  int iterations = 0;
  bool lowprecision = false;
  enum QuantBRLMM::Transformation transform = CCS;
  string transformStr;
  map<string,string>::iterator iter;
  int priorWeight = 40;
  int priorMinCall = 2;
  fillInValue(hetMult, "het-mult", param, doc);
  fillInValue(iterations, "iterations", param, doc);
  fillInValue(iterationThreshold, "iter-thresh", param, doc);
  fillInValue(K, "K", param, doc);
  fillInValue(transformStr, "transform", param, doc);
  transform = transformationForString(transformStr);
  fillInValue(priorWeight, "prior-weight", param, doc);
  fillInValue(priorMinCall, "prior-mincall", param, doc);
  fillInValue(lowprecision, "lowprecision", param, doc);
  QuantBRLMM *brlmm = new QuantBRLMM(hetMult, iterations, iterationThreshold, K, transform,
                                     priorWeight, priorMinCall, lowprecision);
	double score;
	fillInValue(score,"MS",param,doc);
	brlmm->setMaxScore(score);
  return brlmm;
}
