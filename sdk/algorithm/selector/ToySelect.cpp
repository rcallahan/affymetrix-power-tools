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
#include "algorithm/selector/ProbeSelector.h"
//

// this is a toy program to exercise covariate normalization

void executeSelect()
{
	cout << "running executeSelect()...\n";

	int i;

	vector<double> Output;
	vector<int> Input;
	vector< vector<double> > Covar;
	vector<double> Sums;
	vector<double> PriorM;
	vector<double> PriorW;
	PriorM.resize(3);
	PriorW.resize(3);
	PriorM[0] = 0.66;
	PriorM[1] = 0;
	PriorM[2] = -0.66;
	PriorW[0]=PriorW[1]=PriorW[2] = 0.5;

	Input.resize(100);
	Covar.resize(1);
	Covar[0].resize(100);
	for (i=0; i<33; i++)
	{
		Input[i] = 0;
		Covar[0][i] = -.5+i/1000.0;
	}
	for (;i<66; i++)
	{
		Input[i] = 1;
		Covar[0][i] = 0+i/1000.0;
	}
	for (;i<100; i++)
	{
		Input[i] = 2;
		Covar[0][i]=.5+i/1000.0;
	}
	
	LogisticRegression(Sums,Output,Input,Covar,PriorM,PriorW,100,.0001);
	for (i=0; i<Output.size(); i++)
		cout << Output[i] << endl;
  
  //
  cout << "done.\n";
}

int main( int argc, char *argv[])
{
	executeSelect();
}
