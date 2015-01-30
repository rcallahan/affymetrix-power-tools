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

//this is a routine to hold 'probe selection' numerics
// i.e. logistic regression, etc

#include "algorithm/selector/ProbeSelector.h"
//

double LogisticRegression(vector<double> &Summaries,
                          vector<double> &Predictor,
                          const vector<int> &TargetGenotypes,
                          const vector< vector<double> > &Covariates,
			  vector<double> PriorMeans,
			  vector<double> PriorWeights,
                          int loopMax,
                          double Converged) {
    int nobs, npred;
    int i,j;
    double convergeCheck;
    
//vector<double> priorm;
//	vector<double> weightm;
//	priorm.resize(3);
//	weightm.resize(3);
//	priorm[0] = 0.66;
//	priorm[1] = 0;
//	priorm[2] = -0.66;
//	weightm[0] = 1;
//	weightm[1] = 0.5;
//	weightm[2] = 1;

    nobs = TargetGenotypes.size();
    npred = Covariates.size()+1;
    Summaries.clear();
    // logistic regression loop
    // 0,1,2 copies of the A allele in Y
    
    ColumnVector Y(nobs);
    ColumnVector En(nobs);
    ColumnVector Beta(npred);
    Matrix Xmat(nobs,npred);

    // fill in the ones that need initialization
    double safety = 0.5/nobs;
    double ssafety = safety*safety; // tiny value for variance
    for (i=0; i<nobs; i++)
        Y.element(i) = TargetGenotypes[i]+safety*(1-TargetGenotypes[i]);
    for (i=0; i<nobs; i++)
        En.element(i) = 2; // at most 2 alleles
    for (j=0; j<npred; j++)
        Beta.element(j) = 0; // everything het
    for (i=0; i<nobs; i++) {
        Xmat.element(i,0) = 1; // intercept
        for (j=1; j<npred; j++) {
            Xmat.element(i,j) = Covariates[j-1][i];
            //cout << Covariates[j-1][i] << "\t";
        }
    }

    //cout << "initialized" << endl;
    // scratch variables used for IWLS iterations
    Matrix WXmat(nobs,npred);
    ColumnVector Weights(nobs);
    ColumnVector Nu(nobs);
    ColumnVector Mu(nobs);
    ColumnVector Vee(nobs);
    ColumnVector Zee(nobs);
    ColumnVector WZee(nobs);
    Matrix U, V;
    DiagonalMatrix D;
    ColumnVector Fitted(nobs);
    ColumnVector BetaNew(npred);
    ColumnVector Delta(npred);
    
    int loop=0;
    bool done= false;
    while (loop<loopMax && !done) {
        Nu = Xmat * Beta;
        for (i=0; i<nobs; i++) {
            Mu.element(i) = 1/(1+exp(-Nu.element(i)));
            Vee.element(i) = Mu.element(i)*(1-Mu.element(i));
            if (Vee.element(i)<ssafety) { // prevent divide by zero errors and zero weights
                Vee.element(i) = ssafety;
            }
            Zee.element(i) = Nu.element(i) + (Y.element(i)/En.element(i)-Mu.element(i))/Vee.element(i);
            Weights.element(i) = En.element(i) * Vee.element(i);
            Weights.element(i) = sqrt(Weights.element(i));
            WZee.element(i) = Zee.element(i) * Weights.element(i);
            for (j=0; j<npred; j++)
                WXmat.element(i,j) = Xmat.element(i,j) * Weights.element(i);
        }
        //cout << "SVD" << loop << "\t";
        //cout << Vee;
        //cout << WXmat;
        SVD(WXmat,D,U,V);
        
        Fitted = U.t() * WZee;
        //cout << D;
        //cout << "d.i" << "\t";
        //cout << D.i();
        BetaNew = V * ( D.i() * Fitted );
        //cout << "delta" << "\t";
        Delta = Beta-BetaNew;
        
        convergeCheck = Delta.SumSquare();
        Beta = BetaNew;
        loop++;
        if (convergeCheck<Converged)
            done=true;
        //cout << loop << endl;
        //cout << Beta << endl;
    }
    vector<double> gcor,gobs;
    gcor.resize(3);
    gobs.resize(3);
    gcor[0]=gcor[1]=gcor[2] = 0.001;
    gobs[0]=gobs[1]=gobs[2] = 0.001;
    
    double correct,bestguess, tmp;
    int best;
    double  mcor;
    double tcor;
    
    Nu = Xmat * Beta;
    for (i=0; i<nobs; i++) {
        // make "calls"
        Mu.element(i) = 1/(1+exp(-Nu.element(i)));
        Vee.element(i) = Mu.element(i)*(1-Mu.element(i));
        if (Vee.element(i)<ssafety) {// prevent div0
            Vee.element(i) = ssafety;
        }
        best = 1;
        bestguess = 1000000;
        for (j=0; j<3; j++) {
            tmp = (j/En.element(i)-Mu.element(i))/Vee.element(i);
            tmp = fabs(tmp);
            if (j==1) {
                tmp /= 2; // het stretch -logistic
            }
            if (tmp<bestguess) {
                bestguess=tmp;
                best = j;
            }
        }
        gobs[TargetGenotypes[i]] +=1;
        if (best==TargetGenotypes[i]) {
            gcor[TargetGenotypes[i]] +=1;
        }
    }


    mcor = (gcor[0]+gcor[2])/(gobs[0]+gobs[2]);
    Summaries.push_back(mcor);
    tcor = gcor[1]/gobs[1];
    Summaries.push_back(tcor);
    correct =(gcor[0]+gcor[1]+gcor[2])/(gobs[0]+gobs[1]+gobs[2]);
    Summaries.push_back(correct);
    
    // FLD computation
    
    vector<double> sq,tq,ssq,mq;
    double cvar, fld;
    
    tq.resize(3);
    mq.resize(3);
    sq.resize(3);
    ssq.resize(3);
    // sneak prior in here
    tq[0] = PriorWeights[0];
    tq[1] = PriorWeights[1];
    tq[2] = PriorWeights[2];
    sq[0] = PriorMeans[0]*PriorWeights[0];
    sq[1] = PriorMeans[1]*PriorWeights[1];
    sq[2] = PriorMeans[2]*PriorWeights[2];
    ssq[0] = PriorMeans[0]*PriorMeans[0]*PriorWeights[0];
    ssq[1] = PriorMeans[1]*PriorMeans[1]*PriorWeights[1];
    ssq[2] = PriorMeans[2]*PriorMeans[2]*PriorWeights[2];
    for (i=0; i<nobs; i++) {
        sq[TargetGenotypes[i]] += Covariates[0][i];
        tq[TargetGenotypes[i]] += 1;
        ssq[TargetGenotypes[i]] += Covariates[0][i]*Covariates[0][i];
    }
    double epsilon = 0.001;
    cvar = 0.1*epsilon; // smooth
    cvar +=  ssq[0]+ssq[1]+ssq[2];
    cvar -= sq[0]*sq[0]/tq[0];
    cvar -= sq[1]*sq[1]/tq[1];
    cvar -= sq[2]*sq[2]/tq[2];
    cvar /= tq[0]+tq[1]+tq[2]+epsilon;
    cvar = sqrt(cvar); // stddev
    mq[0] = sq[0]/tq[0];
    mq[1] = sq[1]/tq[1];
    mq[2] = sq[2]/tq[2];
    Summaries.push_back(mq[0]);
    Summaries.push_back(mq[1]);
    Summaries.push_back(mq[2]);
    Summaries.push_back(tq[0]);
    Summaries.push_back(tq[1]);
    Summaries.push_back(tq[2]);
    Summaries.push_back(cvar);
    fld = -1*(mq[1]-mq[0])/cvar;
    Summaries.push_back(fld);
    fld = -1*(mq[2]-mq[1])/cvar;
    Summaries.push_back(fld);
    fld = -1*(mq[2]-mq[0])/cvar;
    Summaries.push_back(fld);
    double entropy;
    entropy = tq[0]*log(tq[0])+tq[1]*log(tq[1])+tq[2]*log(tq[2])-(tq[0]+tq[1]+tq[2])*log(tq[0]+tq[1]+tq[2]);
    entropy = -1*entropy;
    Summaries.push_back(entropy);
    // export fitted value for logistic regression
    Predictor.resize(npred);
    
    for (i=0; i<npred; i++)
        Predictor[i]=Beta.element(i);
    
    // compute log-likelihood
    // Mu = pr(1)
    double p,q,LL;
    LL = 0;
    
    for (i=0; i<nobs; i++) {
        // safety to avoid zero
        p = (Mu.element(i)+ssafety)/(1+2*ssafety);
        q = 1-p;
        LL += Y.element(i)*log(p) + (En.element(i)-Y.element(i))*log(q);
    }
    //cout << "exited safely" << endl;
    Summaries.push_back(-2*LL+2*npred);
    //return(correct);
    return(-2*LL+2*npred);
}
