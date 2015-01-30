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

#include "chipstream/BioTypes.h"
#include "chipstream/ClusterZ.h"
#include "chipstream/MatrixUtil.h"
#include "chipstream/ProbeSet.h"
#include "file/TsvFile/TsvFile.h"
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include "newmat.h"
//
#include <cfloat>
#include <cstring>
#include <string>
#include <vector>
//

using namespace affx;

/* Set all our counts to three for the three possible genotypes:
   AA,AB,BB */
ClusterStats::ClusterStats() {
  counts.resize(3);
  aMeans.resize(3);
  aVars.resize(3);
  mMeans.resize(3);
  mVars.resize(3);
  covars.resize(3);
}

/**
 * Write out a text representation of the prior to the file specified.
 * @param prior - Prior to be written.
 * @param out - ostream to write to.
 */
void x_writePriorOut(ClusterPrior &prior, std::ostream &out,
                               const char *fieldDelim, const char *sep) {
  printColumnVector(prior.centers, &out, sep);
  out << fieldDelim;
  printColumnVector(prior.centerVars, &out, sep);
  out << fieldDelim;
  printMatrix(prior.covars, &out, sep);
  out << endl;
}

/**
 * Generate a string representation of a snp model.
 * @return string - text description of snp clusters.
 */
std::string x_getModelString(string m_ProbesetName, ClusterPrior *m_CurrentPrior) {
  int i = 0, j = 0;
  std::string s = m_ProbesetName;
  s += "\t";
  for(i = 0; i < m_CurrentPrior->centers.Nrows() - 1; i++) {
    s += ToStr(m_CurrentPrior->centers.element(i));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->centers.element(i));
  s += "\t";
  for(i = 0; i < m_CurrentPrior->centerVars.Nrows() - 1; i++) {
    s += ToStr(m_CurrentPrior->centerVars.element(i));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->centerVars.element(i));
  s += "\t";
  for(i = 0; i < m_CurrentPrior->covars.Nrows()-1; i++) {
    for(j = 0; j < m_CurrentPrior->covars.Ncols(); j++) {
      s += ToStr(m_CurrentPrior->covars.element(i,j));
      s += ",";
    }
  }
  for(j = 0; j < m_CurrentPrior->covars.Ncols() - 1; j++) {
    s += ToStr(m_CurrentPrior->covars.element(i,j));
    s += ",";
  }
  s += ToStr(m_CurrentPrior->covars.element(i,j));
  return s;
}


/**
 * Write out a text representation of the prior to the file specified.
 * @param prior
 * @param fileName
 */
void x_writePrior(ClusterPrior &prior, const char *fileName) {
  static const char *fieldDelim = "\t", *sep=",";
  ofstream out;
  Fs::mustOpenToWrite(out, fileName);
  out.setf(ios::fixed, ios::floatfield);
  out.precision(10);
  out << "center" << fieldDelim << "var" << fieldDelim << "center.var" << endl;
  x_writePriorOut(prior, out, fieldDelim, sep);
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
ClusterPrior x_clusterPriorFromStrings(const char *fileName, std::string &id, std::string &center,
                                                 std::string &var, std::string &covar) {
  ClusterPrior prior;
  vector<float> vec;
  vector<string> words;

  prior.id = id;

  /* Read in the centers. */
  Util::chopString(center, ',', words);
  vec.clear();
  if(words.size() != 6)
    Err::errAbort("Expecting 6 entires in 'center' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for(unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.centers = columnVectorFromArray(&vec[0], (int)vec.size());

  /* Read in the variances. */
  Util::chopString(var, ',', words);
  vec.clear();
  if(words.size() != 9)
    Err::errAbort("Expecting 9 entires in 'var' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for(unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.centerVars = columnVectorFromArray(&vec[0], (int)vec.size());

  /* Do the covariacne matrix. */
  Util::chopString(covar, ',', words);
  vec.clear();
  if(words.size() != 36)
    Err::errAbort("Expecting 36 entires in 'centers.var' entry in file: " + ToStr(fileName) + " got: " + ToStr(words.size()));
  for(unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  prior.covars = matrixFromArray(&vec[0], (int)vec.size(), 6, 6);
  return prior;

}
/**
 * Read the prior from the test file specified.
 * @param fileName - Name of file to read from.
 * @return - The prior read in from file.
 */
ClusterPrior x_loadPrior(const char *fileName) {
  affx::TsvFile tsv;
  string center, centerVar, var, id="unknown";
  ClusterPrior prior;
  vector<float> vec;
  tsv.open(fileName);
  tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center.var", &centerVar, affx::TSV_BIND_REQUIRED);
  if(tsv.nextLevel(0) != affx::TSV_OK) {
    Err::errAbort("Didnt' get an entry in file: " + ToStr(fileName));
  }
  prior = x_clusterPriorFromStrings(fileName, id, center, var, centerVar);
  tsv.close();
  return prior;
}


/**
 * Load in a number of snp specific priors.
 * @param snpPriors - map of snp probeset names to priors.
 * @param fileName - filename to read priors from.
 */
void x_loadSnpPriors(std::map<std::string, ClusterPrior> &snpPriors,
                          const char *fileName) {
  affx::TsvFile tsv;
  string center, centerVar, var, id="unknown";
  ClusterPrior prior;
  snpPriors.clear();
  tsv.open(fileName);
  tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center.var", &centerVar, affx::TSV_BIND_REQUIRED);
  while(tsv.nextLevel(0) == affx::TSV_OK) {
    prior = x_clusterPriorFromStrings(fileName, id, center, var, centerVar);
    if(snpPriors.find(prior.id) != snpPriors.end()) {
      Err::errAbort("ClusterZ::x_loadSnpPriors() - Id: " + ToStr(prior.id) + " seen multiple times.");
    }
    snpPriors[prior.id] = prior;
  }
  tsv.close();
}


/**
 * Write out a text representation of the model to the file specified.
 * @param model - Model to be written.
 * @param out - ostream to write to.
 */
void x_writeModelOut(ClusterModel &model, ostream &out, const std::string &name,
                               int copies, const char *fieldDelim, const char *sep) {
  out << name << fieldDelim << copies << fieldDelim;
  printMatrix(model.clusterMeans, &out, sep);
  out << fieldDelim;
  printMatrix(model.clusterVars, &out, sep);
  out << endl;
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
ClusterModel x_clusterModelFromStrings(const char *fileName, std::string &id,
                                                 std::string &center, std::string &var) {
  ClusterModel model;
  vector<float> vec(9);
  vector<string> words;
  /* Read in the centers. */
  Util::chopString(center, ',', words);
  vec.clear();
  if(words.size() != 6)
    Err::errAbort("Expecting 6 entires in 'center' entry in file: " + ToStr(fileName) +
                  " got: " + ToStr(words.size()) + " for probest: " + id);
  for(unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  model.clusterMeans = matrixFromArray(&vec[0], (int)vec.size(), 6, 1);

  /* Read in the variances. */
  Util::chopString(var, ',', words);
  vec.clear();
  if(words.size() != 9)
    Err::errAbort("Expecting 9 entires in 'var' entry in file: " + ToStr(fileName) +
                  " got: " + ToStr(words.size()) + " for probest: " + id);
  for(unsigned int i = 0; i < words.size(); i++) {
    vec.push_back(Convert::toDouble(words[i].c_str()));
  }
  model.clusterVars = matrixFromArray(&vec[0], (int)vec.size(), 9, 1);
  return model;
}

/**
 * Load in a number of snp specific models.
 * @param snpModels - map of snp probeset names to models.
 * @param fileName - filename to read models from.
 */
void x_loadSnpModels(std::map<std::string, ClusterModel> &snpModels,const char *fileName) {
  affx::TsvFile tsv;
  string center, var, id, copies;
  ClusterModel model;
  snpModels.clear();
  tsv.open(fileName);
  tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"copies", &copies, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"center", &center, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"var", &var, affx::TSV_BIND_REQUIRED);
  while(tsv.nextLevel(0) == affx::TSV_OK) {
    model = x_clusterModelFromStrings(fileName, id, center, var);
    if(snpModels.find(id) != snpModels.end()) {
      Err::errAbort("ClusterZ::x_loadSnpModels() - Id: " + ToStr(id) + " seen multiple times.");
    }
    snpModels[id + "-" + copies] = model;
  }
  tsv.close();
}
