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

#include "algorithm/spectclust/SpectClust.h"
//
#include "rma/RMA.h"
#include "util/Util.h"
//

using namespace std;

/** 
 * @brief Virtual destructor for a virtual class.
 */
DistanceMetric::~DistanceMetric() {
} 

double DistanceMetric::dist(const Matrix &M, int col1Ix, int col2Ix) const {
  Numeric dist = DotProduct(M.Column(col1Ix), M.Column(col2Ix));
  return dist;
}

  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
double CorrelationMetric::dist(const Matrix &M, int col1Ix, int col2Ix) const {
    std::vector<double> c1(M.Nrows()), c2(M.Nrows());
    for(int i = 0; i < M.Nrows(); i++) {
      c1[i] = M.element(i,col1Ix-1);
      c2[i] = M.element(i,col2Ix-1);
    }
    double dist = 1 + correlation_coeff(c1.begin(), c1.end(), c2.begin());
    return dist;
  }

GuassianRadial::GuassianRadial(double sigma) {
  m_Sigma = sigma;
}

/** 
 * Compute a distance metric between two columns of a
 * matrix. <b>Note that the indexes are *1* based (not 0) as that is
 * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
 * 
 * @param M - Matrix whose columns represent individual items to be clustered.
 * @param col1Ix - Column index to be compared (1 based).
 * @param col2Ix - Column index to be compared (1 based).
 * 
 * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
 */
double GuassianRadial::dist(const Matrix &M, int col1Ix, int col2Ix) const {
  double dist = 0;
  if(col1Ix == col2Ix) 
    return 0;
  ColumnVector V = M.Column(col1Ix) - M.Column(col2Ix);
  dist = V.SumSquare() / (2 * m_Sigma * m_Sigma);
  dist = exp(-1 * dist);
  return dist;
}

AngleMetric::AngleMetric(const Matrix &M) {
  setMatrix(M);
}

/** 
 * Set the matrix. This precalculates the norms for each column so
 * they are only calculated once.
 * @param M - Matrix to be used.
 */
void AngleMetric::setMatrix(const Matrix &M) {
  m_Norms.resize(M.Ncols(),0);
  for(int i=1; i <= M.Ncols(); i++) {
    m_Norms[i-1] = sqrt(M.Column(i).SumSquare());
  }
}
  
ExpAngleMetric::ExpAngleMetric(const Matrix &M, double sigma) :  AngleMetric(M) {
  m_Sigma=sigma;
}
  
/** 
 * Compute a distance metric between two columns of a
 * matrix. <b>Note that the indexes are *1* based (not 0) as that is
 * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
 * 
 * In this case the distnance metric is the exponated 1 - cosine of the angle between the
 * two vectors, this is also sometimes known as the uncentered correlation coefficient
 * and is similar to correlation except that it requires the length of the vectors to 
 * be the same for cos(x,y) to be 1.
 * 
 * If x and y are vectors then:
 * \f$ ExpAngleDist(\vec{x},\vec{y}) = e^{-1 * (1 - cos(\theta))/2\sigma^2} \f$
 * Where:
 * \f$ cos(\theta) = \frac{ \vec{x}^t \cdot \vec{y}}{\|\vec{x}\|\|\vec{y}\|} \f$
 * Where the numerator is the dot product of <b>x</b> and <b>y</b>.
 * @param M - Matrix whose columns represent individual items to be clustered.
 * @param col1Ix - Column index to be compared (1 based).
 * @param col2Ix - Column index to be compared (1 based).
 * 
 * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
 */
double ExpAngleMetric::dist(const Matrix &M, int col1Ix, int col2Ix) const {
  double dist = 0;
  int nRow = M.Nrows();
  for(int rowIx = 0; rowIx < nRow; rowIx++) {
    dist += M.element(rowIx,col1Ix-1) * M.element(rowIx, col2Ix -1);
  }
  dist /= (m_Norms[col1Ix-1] * m_Norms[col2Ix-1]);
  dist = exp(-1 * ((1 - dist)/(2*m_Sigma*m_Sigma)));
  return dist;
}

/** 
 * Compute a distance metric between two columns of a
 * matrix. <b>Note that the indexes are *1* based (not 0) as that is
 * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
 * 
 * In this case the distnance metric is the cosine of the angle between the
 * two vectors, this is also sometimes known as the uncentered correlation coefficient
 * and is similar to correlation except that it requires the length of the vectors to 
 * be the same for cos(x,y) to be 1.
 * 
 * If x and y are vectors then:
 * \f$ angleDist(\vec{x},\vec{y}) = cos(\theta) = \frac{ \vec{x}^t \cdot \vec{y}}{\|\vec{x}\|\|\vec{y}\|} \f$
 * Where the numerator is the dot product of <b>x</b> and <b>y</b>.
 * @param M - Matrix whose columns represent individual items to be clustered.
 * @param col1Ix - Column index to be compared (1 based).
 * @param col2Ix - Column index to be compared (1 based).
 * 
 * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
 */
double AngleMetric::dist(const Matrix &M, int col1Ix, int col2Ix) const {
  double dist = 0;
  int nRow = M.Nrows();
  col1Ix--;
  col2Ix--;
  for(int rowIx = 0; rowIx < nRow; rowIx++) {
    dist += M[rowIx][col1Ix] * M[rowIx][col2Ix];
  }
  dist = (dist / (m_Norms[col1Ix] * m_Norms[col2Ix]) + 1); // cosine is always between 1 and -1 this scales to 2 0.
  return dist;
}

void SpectClust::multByMatrix(ColumnVector &N, ColumnVector &O, const Matrix &M) {
  int nRow = M.Nrows(), nCol = M.Ncols();
  if(!(N.Nrows() == O.Nrows() && M.Nrows() == M.Ncols() && M.Nrows() == N.Nrows())) {
    Err::errAbort("wrong dimensions: " + ToStr(O.Nrows()) + " " + ToStr(N.Nrows()) + " " + ToStr(M.Nrows()));
  }
  for(int rowIx = 0; rowIx < nRow; rowIx++) {
    N[rowIx] = 0;
    for(int colIx = 0; colIx < nCol; colIx++) {
      N[rowIx] += O[colIx] * M[rowIx][colIx];
    }
  }
}

double SpectClust::fast_corr(std::vector<double> &x, std::vector<double> &y) {
  assert(x.size() == y.size());
  int n = x.size();
  double xyS=0, xS=0, x2S=0, yS=0, y2S=0;
  for(int i = 0; i < n; i++) {
    xyS += x[i] * y[i];
    xS += x[i];
    yS += y[i];
    x2S += x[i] * x[i];
    y2S += y[i] * y[i];
  }
  double invN = (double)1/n;
  double cor = (xyS - (invN * xS * yS))/(sqrt(x2S - (invN * xS * xS)) * sqrt(y2S - (invN * yS * yS)));
  return cor;
}

void SpectClust::fillInDistance(SymmetricMatrix &A, const Matrix &M, const DistanceMetric &dMetric, bool expon) {
  A.ReSize(M.Ncols());
  for(int i = 1; i <= M.Ncols(); i++) {
    for(int j = i; j <= M.Ncols(); j++) {
      if(expon) 
        A(j,i) = A(i,j) = exp(-1 * dMetric.dist(M,i,j));
      else
        A(j,i) = A(i,j) = dMetric.dist(M,i,j);
    }
  }
}

double SpectClust::rowMedian(const Matrix &M, int rowIx) {
  std::vector<double> r(M.Ncols(),0);
  for(int i = 0; i < M.Ncols(); i++) {
    r[i] = M.element(rowIx,i);
  }
  double med = median_in_place(r.begin(), r.end());
  return med;
}

void SpectClust::rowMedianDivide(Matrix &M) {
  for(int rowIx = 0; rowIx < M.Nrows(); rowIx++) {
    double med = rowMedian(M, rowIx);
    for(int colIx = 0; colIx < M.Ncols(); colIx++) {
      M.element(rowIx,colIx) /= med;
    }
  }
}

void SpectClust::subColAvg(Matrix &M) {
  ColumnVector Ones(M.Nrows());
  Ones = 1.0;
  // Calculate the average for each column matrix style.
  Matrix ColAvg = Ones.t() * M;
  ColAvg = ColAvg / M.Nrows();
  // Make a matrix with row averages for each column and
  // subtract from original matrix M
  M = M - (Ones * ColAvg);
}

void SpectClust::MatrixScatter(Matrix &M, Matrix &C) {
  Matrix X(M);
  subColAvg(X); // subtract of mean of each column
  C = (X.t() * X); // expected value of (X - xmu)(Y - ymu)
}

void SpectClust::MatrixCov(Matrix &M, Matrix &C) {
  MatrixScatter(M,C);
  C = C / (M.Nrows() - 1); // expected value of (X - xmu)(Y - ymu)
}

void SpectClust::MatrixCor(Matrix &M, Matrix &C) {
  Matrix X(M);
  MatrixCov(M,C);
  DiagonalMatrix V(C.Nrows());
  V << C;
  for(int i = 0; i < V.Nrows(); i++) {
    V.element(i,i) = 1 / sqrt(V.element(i,i)); 
  }
  C = V * C * V;
}

bool SpectClust::MaxEigen(const Matrix &M, double &maxValue, ColumnVector &MaxVec, int maxIterations) {
  double maxDelta = 1e-6;
  bool converged = false;
  int i = 0;
  int nRows = M.Ncols();
  if(M.Ncols() != M.Nrows()) 
    Err::errAbort("MaxEigen() - Can't get eigen values of non square matrices.");
  if(M.Ncols() <= 0) 
    Err::errAbort("MaxEigen() - Must have positive number of rows and columns.");
  ColumnVector V(M.Ncols());
  V = 1.0 / M.Nrows(); // any vector really...
  V = V / Norm1(V);
  MaxVec.ReSize(M.Ncols());
  for(i = 0; i < maxIterations; ++i) {
    //    MaxVec = M * V;
    multByMatrix(MaxVec, V, M);
    double delta = 0;
    double norm = sqrt(SumSquare(MaxVec));
    for(int vIx = 0; vIx < nRows; vIx++) {
      MaxVec.element(vIx) = MaxVec.element(vIx) / norm; // scale so we don't get too big.
    }
    for(int rowIx = 0; rowIx < nRows; rowIx++) {
      delta += fabs((double)MaxVec.element(rowIx) - V.element(rowIx));
    }
    if(delta < maxDelta)  
      break; // we've already converged to eigen vector.
    V = MaxVec;
  }
  if(i < maxIterations) {
    converged = true;
  }
  // calculate approximate max eigen value using Rayleigh quotient (x'*M*x/x'*x).
  Matrix num = (MaxVec.t() * M * MaxVec);
  Matrix denom =  (MaxVec.t() * MaxVec);
  maxValue = num.element(0,0) / denom.element(0,0);
  return converged;
}

void SpectClust::normalizeSum(SymmetricMatrix &A) {
  RowVector Ones(A.Ncols());
  Ones = 1.0;
  // Calculate the sum of each row matrix style. Could this be a Diagonal Matrix instead?
  Matrix RowSum = Ones * (A.t());
  DiagonalMatrix D(A.Ncols());
  D = 0;
  // Take the inverse square root
  for(int i = 1; i <= A.Ncols(); i++) {
    D(i,i) = 1/sqrt(RowSum(1,i));
  }
  Matrix X = D * (A * D);
  A << X;
}

bool SpectClust::findNLargestEvals(const Matrix &M, int numLamda, std::vector<Numeric> &eVals, Matrix &EVec, int maxIterations) {
  bool converged = true;
  EVec.ReSize(M.Ncols(), numLamda);
  eVals.clear();
  eVals.reserve(numLamda);
  Matrix W = M;
    
  for(int i = 1; i <= numLamda; i++) {
    ColumnVector maxVec;
    double maxVal;
    /* Get the maximum eigen vector. */
    converged = MaxEigen(W, maxVal, maxVec, maxIterations) && converged;
    EVec.Column(i) << maxVec;
    eVals.push_back(maxVal);
     
    /* Now subtract of the largest eigen value to get the next
       largest in next iteration. */
    Matrix ToSub = maxVal * (maxVec * maxVec.t());
    W = W - ToSub;
  }
  return converged;
}

double SpectClust::calcNormalizedCut(const Matrix &D, 
                                     const std::vector<std::pair<double,int> > &indices, 
                                     int cut) {
  double assocA = 0, assocB = 0, cutAB = 0;
  // count up weight of A nodes connected to B nodes. 
  for(int aIx = 0; aIx < cut; aIx++) {
    int Aindex = indices[aIx].second;
    for(int colIx = cut; colIx < indices.size(); colIx++) {
      int Bindex = indices[colIx].second;
      if(Aindex == Bindex) {
        Err::errAbort("How can " +ToStr(Aindex) + " be in both the a and b index?");
      }
      cutAB += D.element(Aindex,Bindex);
    }
  }
  // Count up connectivity of A to entire graph.
  for(int aIx = 0; aIx < cut; aIx++) {
    int Aindex = indices[aIx].second;
    for(int colIx = 0; colIx < D.Ncols(); colIx++) {
      assocA += D.element(Aindex, colIx);
    }
  }
  // Count up connectivity of B to entire graph.
  for(int bIx = cut; bIx < D.Nrows(); bIx++) {
    int Bindex = indices[bIx].second;
    for(int colIx = 0; colIx < D.Ncols(); colIx++) {
      assocB += D.element(Bindex, colIx);
    }
  }
  double nCut = DBL_MAX;
  if(assocA != 0 && cutAB != 0 && assocB != 0)
    nCut = (cutAB/assocA + cutAB/assocB);
  return nCut;
}


void SpectClust::calcNormalizedCutValues(const Matrix &D, const std::vector<std::pair<double,int> > &indices, 
                                         std::vector<double> &cutVals) {
  cutVals.resize(indices.size());
  std::fill(cutVals.begin(), cutVals.end(), DBL_MAX);
  double assocA = 0;
  double assocB = D.Sum();
  double currentNCut = 0;
  double lastNCut = 0;
  int nCol = D.Ncols();
  for(int cut = 0; cut < indices.size() -1; cut++) {
    double rowSum = 0;
    for(int colIx = 0; colIx < nCol; colIx++) {
      rowSum += D.element(indices[cut].second, colIx);
    }
    assocA += rowSum;
    assocB -= rowSum;
    double cutChange = 0;
    for(int colIx = 0; colIx < cut; colIx++) {
      cutChange -= D.element(indices[colIx].second, indices[cut].second);
    }
    for(int colIx = cut+1; colIx < nCol; colIx++) {
      cutChange += D.element(indices[cut].second, indices[colIx].second);
    }
    currentNCut = lastNCut + cutChange;
    if(assocA == 0 || assocB == 0 || currentNCut == 0) {
      cutVals[cut] = DBL_MAX;
    }
    else {
      double nVal = (currentNCut / assocA) + (currentNCut / assocB);
      cutVals[cut+1] = nVal; // +1 is to be compatible with existing calcNormalizedCut() conventions
    }
    lastNCut = currentNCut;
  }
}


bool SpectClust::findNLargestSymEvals(const SymmetricMatrix &W, int numLamda, std::vector<Numeric> &eVals, Matrix &EVec) {
  bool converged = false;
  eVals.clear();
  eVals.reserve(numLamda);
  DiagonalMatrix D(W.Ncols());
  Matrix E;
  try {
    EigenValues(W, D, E);
    converged = true;
    EVec.ReSize(W.Ncols(), numLamda);
    int count = 0;
    for(int i = W.Ncols(); i > W.Ncols() - numLamda; i--) {
      eVals.push_back(D(i));
      EVec.Column(++count) << E.Column(i);
    }
  }
  catch(const Exception &e) {
    Err::errAbort("Exception: " + ToStr(e.what()));
  }
  catch(...) {
    Err::errAbort("Yikes couldn't calculate eigen vectors.");
  }
  return converged;
}

void SpectClust::orderIntraClusterCorr(const Matrix &M, vector<int> &clusters, int numClusters, std::vector<double> &clusterCorr) {
  std::vector<std::pair<double,int> > inCorr;
  for(int clustIx = 0; clustIx < numClusters; clustIx++) {
    vector<double> clustCorr;
    for(int col1 = 0; col1 < M.Ncols(); col1++) {
      if(clusters[col1] == clustIx) {
        vector<double> c1Dat;
        vector<double> c2Dat;
        colToVector(M, col1, c1Dat);
        for(int col2 = col1 + 1; col2 < M.Ncols(); col2++) {
          if(clusters[col2] == clustIx) {
            colToVector(M, col2, c2Dat);
            double corr = fast_corr(c1Dat,c2Dat);
            clustCorr.push_back(corr);
          } // if
        } // for col2
      } // if clusters[col1]
    } // for col1
    std::pair<double,int> p;
    if(clustCorr.size() != 0) 
      p.first = median_in_place(clustCorr.begin(), clustCorr.end());
    else 
      p.first = -2; // lower than any possible correlation
    p.second = clustIx;
    inCorr.push_back(p);
  } // for clustIx
  std::sort(inCorr.begin(), inCorr.end(), PairLess());
  clusterCorr.resize(inCorr.size());
  for(int i = 0; i < inCorr.size(); i++) {
    clusterCorr[i] = inCorr[i].first;
  }
  vector<int> newClusters(clusters.size(), -1);
  for(int icIx = 0; icIx < inCorr.size(); icIx++) {
    for(int clustIx = 0; clustIx < clusters.size(); clustIx++) {
      if(clusters[clustIx] == inCorr[icIx].second) {
        newClusters[clustIx] = icIx;
      }
    }
  }
  clusters = newClusters;
}

void SpectClust::colToVector(const Matrix &M, int colIx , std::vector<double> &v) {
  v.clear();
  v.reserve(M.Nrows());
  for(int i =0; i < M.Nrows(); i++) {
    v.push_back(M.element(i,colIx));
  }
}

void SpectClust::partitionClusters(const Matrix &D, const Matrix &E, std::vector<Numeric> &eVals, int numClusters,
                                   std::vector<int> &clusters, int hardMin, double cutVal, double margin) {
  clusters.clear();
  clusters.resize(D.Ncols(), -1);
  std::vector<std::pair<double,int> > indices;
  /* Load up all the values and their indexes. */
  for(int i = 0; i < D.Ncols(); i++) {
    std::pair<double,int> p;
    p.first = E.element(i,1);
    p.second = i;
    indices.push_back(p);
  }
  std::sort(indices.begin(), indices.end(), PairLess());
  vector<double> indexNVals;
  calcNormalizedCutValues(D, indices, indexNVals);
  double minNormCut = DBL_MAX;
  int minNCutIx = 0;
  if(cutVal == DBL_MAX) {
    for(int i = 0; i < indexNVals.size(); i++) {
      double nCutI = indexNVals[i];
      if(nCutI < minNormCut) {
        minNCutIx = i;
        minNormCut = nCutI;
      }
    }
    cutVal = minNormCut;
  }
  else {
    for(int i = 0; i < indices.size() - 1; i++) {
      if(indices[i].first <= cutVal && indices[i+1].first >= cutVal) {
        minNCutIx = i;
      }
    }
  }
  clusters.resize(indices.size());
  fill(clusters.begin(), clusters.end(), -1);
  int maxNCutIx = minNCutIx;
  int diff = indices.size() - minNCutIx;
  diff = (int) (margin*diff);
  maxNCutIx = Max((int)indices.size() - diff, maxNCutIx); 
  if(minNCutIx >= hardMin)
    minNCutIx = Max((int)(margin*minNCutIx), hardMin);
  for(int i = 0; i < indices.size(); i++) {
    if(i <= minNCutIx) {
      clusters[indices[i].second] = 0;
    }
    else if(i > maxNCutIx) {
      clusters[indices[i].second] = 1;
    }
  }
}

void SpectClust::AICcRmaMedianCluster(const Matrix &PM, int numClusters, const std::vector<int> &clusters, 
                                   std::vector<double> &AICc, std::vector<double> &BIC, bool doLog) {
  vector<vector<vector<float> > > clusterDat(numClusters); // first index is cluster, then probe, then chip.
  vector<vector<float> > allDat;
  AICc.resize(2);
  fill(AICc.begin(), AICc.end(), 0);
  BIC.resize(2);
  fill(BIC.begin(), BIC.end(), 0);
  vector<int> clusterCounts(numClusters, 0);
  int okCount = 0;
  if(PM.Nrows() != clusters.size()) 
    Err::errAbort("SpectClust::AICcRmaMedianCluster() - Number of PM rows doesn't equal cluster vector size.");
  // setup PM matrices for each cluster and one for all the data.
  for(int probeIx = 0; probeIx < clusters.size(); probeIx++) {
    if(clusters[probeIx] >= 0) {
      clusterCounts[clusters[probeIx]]++;
      okCount++;
      vector<float> dat(PM.Ncols());
      for(int chipIx = 0; chipIx < PM.Ncols(); chipIx++) {
        dat[chipIx] = PM.element(probeIx, chipIx);
      }
      clusterDat[clusters[probeIx]].push_back(dat);
      allDat.push_back(dat);
    }
  }
  if(okCount == 0) 
    Err::errAbort("SpectClust::AICcRmaMedianCluster() - No data provided.");
  // calculate the residual sum of square error from an RMA median polish fit over the clusters.
  double clustSumSquare = 0;
  vector<float> colEst;
  vector<float> rowEst;
  for(int clustIx = 0; clustIx < numClusters; clustIx++) {
    if(clusterCounts[clustIx] == 0) // no residuals if this cluster is empty.
      continue;
    RMA::medianPolishPsetFromMatrix(clusterDat[clustIx], clusterDat[clustIx].size(), clusterDat[clustIx][0].size(), colEst, rowEst, doLog);
    for(int rowIx = 0; rowIx < clusterDat[clustIx].size(); rowIx++) {
      for(int colIx = 0; colIx < clusterDat[clustIx][0].size(); colIx++) {
        clustSumSquare += clusterDat[clustIx][rowIx][colIx] * clusterDat[clustIx][rowIx][colIx];
      }
    }
  }
  double allSumSquare = 0;
  // calculate residual sum of square error from an RMA median polish over entire data set.
  RMA::medianPolishPsetFromMatrix(allDat, allDat.size(), allDat[0].size(), colEst, rowEst, doLog);
  for(int rowIx = 0; rowIx < allDat.size(); rowIx++) {
    for(int colIx = 0; colIx < allDat[0].size(); colIx++) {
      allSumSquare += allDat[rowIx][colIx] * allDat[rowIx][colIx];
    }
  }
  /* Number of parameters for each model. 

     For the single cluster model we have (NumProbes - 1) parameters
     for the probe effects (as median of probe effects is constrained
     to be zero) and additional NumChip parameters for the chip effects:
     Total = NumProbes - 1 + NumChips.

     For the two cluster model we have NumProbes_A - 1 for the first cluster
     probe effects and NumProbes_B -1 for the second cluster probe effects and additional
     NumChip parameters for the chip effects for both the A and B clusters:
     Total = NumProbes_A - 1 + NumProbes_B - 1 + (2 * NumChips) = NumProbes - 2 + (2 * NumChips)
  */
  int numParam = PM.Ncols() + okCount - 1;
  int clustNumParam = (PM.Ncols() * numClusters) + okCount - numClusters;
  int N = PM.Ncols() * okCount;
  if(okCount > numClusters+2) {
    double Izero = log(allSumSquare/N); 
    double Ione = log(clustSumSquare/N);
    double Uzero, Uone;
    Uzero = log(allSumSquare/(N -  numParam)); 
    Uone = log(clustSumSquare/(N - (clustNumParam)));
    AICc[0] = Izero + ((double)(N + numParam)/(N - numParam - 2));
    AICc[1] = Ione + ((N + (clustNumParam))/(double)(N - (clustNumParam) - 2));
    BIC[0] =  N * Izero + (numParam * log((double)N));
    BIC[1] =  N * Ione + (clustNumParam * log((double)N));
  }
}
