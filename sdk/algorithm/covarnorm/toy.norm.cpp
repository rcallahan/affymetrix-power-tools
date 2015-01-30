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

#include "util/AptVersionInfo.h"
#include "algorithm/covarnorm/covarnorm.h"
//

// this is a toy program to exercise covariate normalization

void executecube()
{
	// try out the cubic regression
	vector<double> C,S,W;
	vector<double> Eqn,Pred;
	float ti;
	float i;
	for (i=0; i<10; i++)
	{
		ti=i;
		S.push_back(ti);
		ti=i*i*i + 3*i*i+7*i+1;
		C.push_back(ti);
		W.push_back(1);
	}
	FitWeightedCubic(C,S,W,Eqn,Pred);
	td(C);
	td(S);
	td(W);
	td(Eqn);
	td(Pred);
}

void executenorm()
{
	NormalizationFrame Test;
	vector<double> C,S;
	float ti;
	float i;
	for (i=0; i<11; i++)
	{
		ti = (i-90)/90;
		C.push_back(ti);
		S.push_back(i);
		ti = (i+78)/90;
		C.push_back(ti);
		S.push_back(i);
		ti = (i-5.5)/90;
		C.push_back(ti);
		S.push_back(i);
	}
	td(C);
	td(S);

	Test.FillInContrast(C);
	Test.FillInCovariates(S);
	Test.FitEM();
	Test.ThisFit.TargetCenter.resize(3);
	cout << "AA:\t";
	td(Test.ThisFit.fitAA);
	cout << "AB:\t";
	td(Test.ThisFit.fitAB);
	cout << "BB:\t";
	td(Test.ThisFit.fitBB);
	cout << "sigma:\t";
	td(Test.ThisFit.sigma);

	Test.ThisFit.TargetCenter[0]=-.66;
	Test.ThisFit.TargetCenter[1] = 0;
	Test.ThisFit.TargetCenter[2] = .66;
	
	vector<double> dummy, dummypred;
	dummy.push_back(0); // dummy covariate
	dummy.push_back(0);
	dummy.push_back(0);
	dummy.push_back(0);
	for (int blah=0; blah<C.size(); blah++)
	{
		dummypred.push_back( Test.ThisFit.MakePrediction(C[blah],dummy));
	}
	/*
    cout<< "C" << "\t";
    td(C);
    cout<< "CX" << "\t";
    td(dummypred);
  */
}

int main( int argc, char *argv[])
{
	executenorm();
	// executecube();
}
