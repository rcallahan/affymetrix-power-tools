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

//
#include "chipstream/AnalysisStreamExpPcaSel.h"
//
#include "algorithm/spectclust/SpectClust.h"
#include "portability/affy-base-types.h"
#include "util/Fs.h"
//
#include "newmat.h"
#include "newmatio.h"
//
#include <cfloat>
#include <cmath>
#include <set>
#include <vector>

using namespace std;

/*
questions for Jim:
- Subtract median and divide by MAD for probes? madstdz()? Also for chips?
- What does =~ do?
*/

enum AnalysisStreamExpPcaSel::InfoFilter
AnalysisStreamExpPcaSel::stringToFilter(const std::string& filter) {
  if(filter=="none") {
    return NoFilter;
  }
  else if(filter=="bic") {
    return BIC;
  }
  else if(filter=="aic") {
    return AIC;
  }
  //
  Err::errAbort("Don't recognize InfoFilter of type: '" + ToStr(filter) + "'");
  return NoFilter; // for compiler...
}

  /** Constructor. */
AnalysisStreamExpPcaSel::AnalysisStreamExpPcaSel(bool doLog,
                                                 bool doCorr,
                                                 bool doDebug,
                                                 const std::string &infoFilter,
                                                 int hardMin, double minProportion, bool qnormOnly)
  : AnalysisStreamExpression()
{
  setupSelfDoc(*this);
  m_DoLog = doLog;
  m_DoCorrelation = doCorr;
  m_Debug = doDebug;
  m_HardMinimum = hardMin;
  m_HardProportion = minProportion;
  m_InfoFilter = AnalysisStreamExpPcaSel::stringToFilter(infoFilter.c_str());
  m_QuantNormOnly = qnormOnly;
  m_QNorm = NULL;
  setOptValue("log", m_DoLog);
  setOptValue("corr", m_DoCorrelation);
  setOptValue("debug", m_Debug);
  setOptValue("info-criterion", infoFilter);
  setOptValue("min-percent", ToStr(m_HardProportion));
  setOptValue("hard-min", ToStr(m_HardMinimum));
  setOptValue("qnorm-only", ToStr(m_QuantNormOnly));
}

AnalysisStreamExpPcaSel::~AnalysisStreamExpPcaSel() {
  Fs::carefulClose(m_Log);
  delete m_QNorm;
}

void AnalysisStreamExpPcaSel::setQuantNorm(SketchQuantNormTran *qnorm) {
  m_QNorm = qnorm;
}

void AnalysisStreamExpPcaSel::registerChipStreamObjs(IntensityReader &reader) {
  if(getChipStreamHead() != NULL)
    reader.registerStream(getChipStreamHead());
  if(m_QuantNormOnly) {
    if(m_QNorm == NULL) {
      Err::errAbort("AnalysisStreamExpPcaSel::registerChipStreamObjs() - Can't have NULL SketchQuantNormTran and set m_QuantNormOnly");
    }
    reader.registerStream(m_QNorm);
  }
}

void AnalysisStreamExpPcaSel::fillInPmData(Matrix &mat,
                                           std::vector<probeid_t> &probeIds,
                                           ProbeSetGroup &psGroup,
                                           std::vector<ChipStream *> &cStream,
                                           IntensityMart &iMart,
                                           bool doLog) {
  unsigned int psIx = 0;
  unsigned int atomPmCount = 0;
  unsigned int chipCount = iMart.getCelFileCount();
  /* Loop through and fill in the data as coming from the intensity mart. */
  for(psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    if(ps == NULL)
      continue;
    for(uint32_t atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom &atom = *(ps->atoms[atomIx]);
      unsigned int channelIx = atom.getChannelCode();
      for(uint32_t probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
        Probe *p = atom.probes[probeIx];
        if(p->type == Probe::PMST || p->type == Probe::PMAT) {
          probeid_t probeId = p->id;
          probeIds[atomPmCount] = p->id;
          for(uint32_t chipIx = 0; chipIx < chipCount; chipIx++) {
            /* PM transformation. */
            double intensity = QuantMethod::transformPrimaryData(probeId, chipIx, iMart, cStream, channelIx);
            ///@todo update mat to handle multi channel expression probesets
            mat.element(atomPmCount, chipIx) = intensity;
          }
          atomPmCount++;
        }
      }
    } /* atoms. */
  } /* probe sets. */
}

// void AnalysisStreamExpPcaSel::subColAvg(Matrix &M) {
//   ColumnVector Ones(M.Nrows());
//   Ones = 1.0;
//   // Calculate the average for each column matrix style.
//   Matrix ColAvg = Ones.t() * M;
//   ColAvg = ColAvg / M.Nrows();
//   // Make a matrix with row averages for each column and
//   // subtract from original matrix M
//   M = M - (Ones * ColAvg);
// }

// void AnalysisStreamExpPcaSel::MatrixScatter(Matrix &M, Matrix &C) {
//   Matrix X(M);
//   SpectClust::subColAvg(X); // subtract of mean of each column
//   C = (X.t() * X); // expected value of (X - xmu)(Y - ymu)
// }


// void AnalysisStreamExpPcaSel::MatrixCor(Matrix &M, Matrix &C) {
//   Matrix X(M);
//   SpectClust::MatrixCov(M,C);
//   DiagonalMatrix V(C.Nrows());
//   V << C;
//   for(int i = 0; i < V.Nrows(); i++) {
//     V.element(i,i) = 1 / sqrt(V.element(i,i));
//   }
//   C = V * C * V;
// }

// void AnalysisStreamExpPcaSel::MatrixCov(Matrix &M, Matrix &C) {
//   MatrixScatter(M,C);
//   C = C / (M.Nrows() - 1); // expected value of (X - xmu)(Y - ymu)
// }

// void AnalysisStreamExpPcaSel::MaxEigen(Matrix &M, double &maxValue, ColumnVector &MaxVec) {
//   double maxDelta = 1e-6;
//   int maxIterations = 50;
//   int i = 0;
//   if(M.Ncols() != M.Nrows())
//     Err::errAbort("MaxEigen() - Can't get eigen values of non square matrices.");
//   if(M.Ncols() <= 0)
//     Err::errAbort("MaxEigen() - Must have positive number of rows and columns.");
//   ColumnVector V(M.Ncols());
//   V = 1.0 / V.Nrows(); // any vector really...
//   V = V / Norm1(V);
//   for(i = 0; i < maxIterations; i++) {
//     //    MaxVec = M * V;
//     SpectClust::multByMatrix(MaxVec, V, M);
//     // cout << "MaxVec: " << endl << MaxVec << endl;
//     double norm = sqrt(SumSquare(MaxVec));
//     MaxVec = MaxVec / norm; // scale so we don't get too big.
//     // cout << "Norm: " << norm << endl;
//     // cout << "Normed MaxVec: " << endl << MaxVec << endl;
//     double delta = SumAbsoluteValue(MaxVec - V);
//     if(delta < maxDelta)
//       break; // we've already converged to eigen vector.
//     V = MaxVec;
//   }
//   // calculate approximate max eigen value using Rayleigh quotient (x'*M*x/x'*x).
//   Matrix num = (MaxVec.t() * M * MaxVec);
//   Matrix denom =  (MaxVec.t() * MaxVec);
//   maxValue = num.element(0,0) / denom.element(0,0);
// }


void AnalysisStreamExpPcaSel::doFeatureSelection(ColumnVector &W, std::set<probeid_t> &goodIds,
                                                 Matrix &PM, std::vector<probeid_t> &probeIds,
                                                 bool doCorr, std::ofstream *out, const char *psName) {
  double mEigVal = 0; // Maximum eigen value
  int atomPmCount = PM.Ncols();
  Matrix Cov(atomPmCount, atomPmCount);
  if(doCorr) {
    SpectClust::MatrixCor(PM, Cov);
  }
  else {
    SpectClust::MatrixScatter(PM, Cov);
  }
  ColumnVector mEig(atomPmCount); // Max eigen vector.
  try {
    SpectClust::MaxEigen(Cov, mEigVal, mEig, 50);
    // If eigenvector is negative overall make it positive.
    Matrix ToSum = Cov * mEig; // Why can't we just sum the vector
                               // itself? If it is really an
                               // eigenvector, should just be a constant
                               // times sum and not change sign
    if(ToSum.Sum() < 0) {
      mEig = mEig * -1;
    }
  }
  catch(BaseException &e) {
    Err::errAbort("Newmat Err: " + ToStr(e.what()));
  }
  catch(...) {
    Err::errAbort("Unknown newmat Error.");
  }

  // Set negative weights to zero.
  for(int i = 0; i < mEig.Nrows(); i++) {
    mEig.element(i) = Max((double)mEig.element(i), 0.0);
  }
  // Divide all weights by the norm of the vector.
  mEig = mEig / sqrt(mEig.SumSquare());

  // Use each component of the eigenvector as the weight for a particular probe.
  W = mEig;
  for(int i = 0; i < W.Nrows(); i++) {
    W.element(i) = W.element(i) * W.element(i);
  }
  vector<double> weights(atomPmCount);
  for(int i = 0; i < W.Nrows(); i++) {
    weights[i] = W.element(i);
  }
  if(out != NULL && out->is_open()) {
    for(int i = 0; i < weights.size(); i++) {
      (*out) << psName << "\t" << probeIds[i] << "\t" << weights[i] << endl;
    }
  }
  // Find the minimum weight to allow, want to grab the first N percent.
  sort(weights.begin(), weights.end());
  vector<double>::reverse_iterator iter;
  double sum = 0;
  double minVal = 1; // can't be any smaller than 1.
  for(iter = weights.rbegin(); iter != weights.rend(); ++iter) {
    sum += *iter;
    minVal = *iter;
    if(sum >= .9)  // get 90% of the weights, then leave the rest.
      break;
  }
  // mark probes to use as those that are above our minimum weight.
  for(int i = 0; i < W.Nrows(); i++) {
    if(W.element(i) >= minVal) {
      goodIds.insert(probeIds[i]);
    }
  }
}

void AnalysisStreamExpPcaSel::log2Matrix(Matrix &M) {
  static const double minVal = DBL_MIN;
  for(int i = 0; i < M.Nrows(); i++) {
    for(int j = 0; j < M.Ncols(); j++) {
      if(M.element(i,j) < 0) {
        Err::errAbort("log2Matrix() - Can't take log2 of negative values.");
      }
      M.element(i,j) = log2(Max(minVal, (double)M.element(i,j)));
    }
  }
}

bool AnalysisStreamExpPcaSel::fillInSelectProbes(ProbeSetGroup &selectGroup, 
						 std::vector<ChipStream *> &cStream, 
						 IntensityMart &iMart, 
						 ProbeSetGroup &psGroup, 
						 std::vector<double> &confVals) {
  int chipCount = iMart.getCelDataSetCount();
  int atomPmCount = psGroup.countPmProbes();
  int newAtomPmCount = atomPmCount;
  bool success = true;
  if(atomPmCount < m_HardMinimum) {
    success = false;
  }
  else {
    vector<bool> probesToUse(atomPmCount, false);
    vector<probeid_t> probeIds(atomPmCount, -1); // Probe ids in order as they are in psGroup
    set<probeid_t> goodIds;
    Matrix PM(atomPmCount, chipCount); // rows are probes, columns are chips
    fillInPmData(PM, probeIds, psGroup, cStream, iMart, m_DoLog);
    Matrix InfoM;
    if(m_InfoFilter != NoFilter) {
      InfoM = PM;
    }
    if(m_DoLog) {
      log2Matrix(PM);
    }
    //  cout << endl << "Original: " << psGroup.name << endl << PM << endl;
    PM = PM.t();
    ColumnVector W;
    string tempName;
    tempName = psGroup.name;
    doFeatureSelection(W, goodIds, PM, probeIds, m_DoCorrelation, &m_ProbeWeights, tempName.c_str());
    bool infoGood = true;
    vector<double> aicVals(2, 0);
    vector<double> bicVals(2, 0);
    // if we're deciding whether or not to do feature selection based
    // on information critera do so here.
    if(m_InfoFilter != NoFilter) {
      vector<int> clusters(probeIds.size(), 1);
      for(int i = 0; i < probeIds.size(); i++) {
        if(goodIds.find(probeIds[i]) != goodIds.end()) {
          clusters[i] = 0;
        }
      }
      vector<double> aicVals, bicVals;
      SpectClust::AICcRmaMedianCluster(InfoM, 2, clusters, aicVals, bicVals, true);
      //      SpectClust::AICcMedianCluster(M, 2, clusters, aicVals, bicVals);
      if(m_InfoFilter == AnalysisStreamExpPcaSel::BIC)
        confVals = bicVals;
      else if(m_InfoFilter == AnalysisStreamExpPcaSel::AIC)
        confVals = aicVals;
      else
        Err::errAbort("AnalysisStreamExpSpectSel::doFeatureSelection() - Don't recognize info filter: '" + ToStr(m_InfoFilter) + "'");
      if(confVals[0] <= confVals[1]) {
        infoGood = false;
      }
    }
    double percent = (double) goodIds.size() / atomPmCount;
    if(goodIds.size() < m_HardMinimum || percent < m_HardProportion) {
      fill(confVals.begin(), confVals.end(), 0);
    }
    if(goodIds.size() >= m_HardMinimum && infoGood) {
      makeSelectProbesetGroup(selectGroup, goodIds, psGroup);
      newAtomPmCount = selectGroup.countPmProbes();
    }
    else {
      success = false;
    }
  }
  return success;
}

 /**
 * Do the analysis for a particular group of probe sets.
 *
 * @param psGroup - Collection of probe sets to get probes from.
 * @param layout - How probes/probesets are laid out on chip.
 * @param iMart - Object containing raw data values for all chips.
 * @param doReport - Should the quantification report object be called?
 * @param alleleSummaryOnly - this is a parameter that makes sense in the base class AnalysisStream::doAnalysis, but not here.  Included here only to make the inheritance work.  Feel free to ignore.
 *
 * @return true if success, false otherwise.
 */
bool AnalysisStreamExpPcaSel::doAnalysis(ProbeSetGroup &psGroup,
                                         IntensityMart &iMart, bool doReport, bool alleleSummaryOnly) {
  if(m_Debug && !m_AllProbes.is_open()) {
    /// @todo use TsvFile
    Fs::mustOpenToWrite(m_AllProbes, Fs::join(m_OutPrefix,getName() + ".pca-select.data.txt"));
    Fs::mustOpenToWrite(m_UsedProbes,Fs::join(m_OutPrefix,getName() + ".pca-select.useddata.txt"));
    m_AllProbes << "probeset\tprobe";
    m_UsedProbes << "probeset\tprobe";

    std::vector<std::string> celFiles = iMart.getCelFileNames();
    for(int i = 0; i < celFiles.size(); i++) {
      m_AllProbes << "\t" << celFiles[i];
      m_UsedProbes << "\t" << celFiles[i];
    }
    m_AllProbes << endl;
    m_UsedProbes << endl;

    Fs::mustOpenToWrite(m_ProbeWeights,Fs::join(m_OutPrefix,getName() + ".pca-select.weights.txt"));
    m_ProbeWeights << "probeset\tprobe\tweight" << endl;
  }
  if(!m_Log.is_open()) {
    Fs::mustOpenToWrite(m_Log,Fs::join(m_OutPrefix,getName() + ".pca-select.report.txt"));
    m_Log << "probeset_id\ttotal\tused\torig_probes\tused_probes" << endl;
  }
  bool success = true;
  ProbeSetGroup psSelectGroup, *toUse = NULL;
  vector<double> confVals(2,0);
  // Use different chipstream vector if we're using our own normalization...
  vector<ChipStream *> selfQNorm(1);
  vector<ChipStream *> *stream = NULL;
  if(m_QuantNormOnly) {
    if(m_QNorm == NULL) {
      Err::errAbort("AnalysisStreamExpPcaSel::doAnalysis() - Can't have NULL SketchQuantNormTran and set m_QuantNormOnly");
    }
    selfQNorm[0] = m_QNorm;
    stream = &selfQNorm;
  }
  else {
    stream = &m_CStreams;
  }
  // Do the actual selection
  if(fillInSelectProbes(psSelectGroup, *stream, iMart, psGroup, confVals)) {
    toUse = &psSelectGroup;
  }
  else {
    toUse = &psGroup;
  }
  if(m_Log.is_open()) {
    AnalysisStreamExpPcaSel::reportProbesUsed(m_Log, *toUse, psGroup, confVals);
  }
  if(m_AllProbes.is_open()) {
    reportProbeLevelData(m_AllProbes, psGroup, iMart, *stream);
    reportProbeLevelData(m_UsedProbes, *toUse, iMart, *stream);
  }
  if(m_QMethod->setUp(*toUse, iMart, m_CStreams, *m_PmAdjust)) {
    m_QMethod->computeEstimate();
    if(doReport) {
      for(unsigned int i = 0; i < m_Reporters.size(); i++) {
        m_Reporters[i]->report(*toUse, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
      }
    }
  }
  else {
    Verbose::out(5, "Warning setup failed for name: " + ToStr(toUse->name));
    success = false;
  }
  if(!success && doReport) {
    for(unsigned int i = 0; i < m_Reporters.size(); i++) {
      m_Reporters[i]->reportFailure(*toUse, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
    }
  }
  return success;
}

void
AnalysisStreamExpPcaSel::reportProbeLevelData(std::ofstream &out,
                                              ProbeSetGroup &psGroup,
                                              IntensityMart &iMart,
                                              std::vector<ChipStream *> &cStream) {
  std::string name;
  name = psGroup.name;
  for(int psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    ProbeSet *ps = psGroup.probeSets[psIx];
    for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom *a = ps->atoms[atomIx];
      unsigned int channelIx = a->getChannelCode();
      for(int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
        Probe *p = a->probes[probeIx];
        if(name != "") {
          out << name << "\t";
        }
        else {
          out << ps->name << "\t";
        }
        out << p->id;

        for(int chipIx = 0; chipIx < iMart.getCelFileCount(); chipIx++) {
          double intensity = QuantMethod::transformPrimaryData(p->id, chipIx, iMart, cStream, channelIx);
          out << "\t" << intensity;
        }
        out << endl;
      }
    }
  }
}

void AnalysisStreamExpPcaSel::reportProbesUsed(std::ofstream &out,
                                               ProbeSetGroup &used,
                                               ProbeSetGroup &orig,
                                               std::vector<double> &confVals) {
  int atomPmCount = orig.countPmProbes();
  int newAtomPmCount = used.countPmProbes();
  out << used.name;
  out << "\t" << atomPmCount << "\t" << newAtomPmCount;
  vector<probeid_t> usedIds;
  vector<probeid_t> origIds;
//   if(confVals.size() != 0) {
//     for(int i = 0; i < confVals.size(); i++) {
//       out << "\t" << confVals[i];
//     }
//   }
  for(int psIx = 0; psIx < used.probeSets.size(); psIx++) {
    ProbeSet *ps = used.probeSets[psIx];
    for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom *a = ps->atoms[atomIx];
      for(int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
        Probe *p = a->probes[probeIx];
        usedIds.push_back(p->id);
      }
    }
  }
  for(int psIx = 0; psIx < orig.probeSets.size(); psIx++) {
    ProbeSet *ps = orig.probeSets[psIx];
    for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom *a = ps->atoms[atomIx];
      for(int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
        Probe *p = a->probes[probeIx];
        origIds.push_back(p->id);
      }
    }
  }
  std::sort(usedIds.begin(), usedIds.end());
  std::sort(origIds.begin(), origIds.end());
  out << "\t";
  for(int i = 0; i < origIds.size(); i++) {
    out << FROM_ZEROBASED(origIds[i]) << ",";
  }
  out << "\t";
  for(int i = 0; i < usedIds.size(); i++) {
    out << FROM_ZEROBASED(usedIds[i]) << ",";
  }
  out << endl;
}
