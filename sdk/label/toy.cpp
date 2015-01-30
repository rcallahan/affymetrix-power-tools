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
#include "label/snp.label.h"
//
#include "chipstream/BioTypes.h"
//
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

void five(vector<double> &z)
{
	z[0]= -.70;
	z[1]= -.75;
	z[2]= -.70;
	z[3]= -.76;
	z[4]= -.65;
	z[5]= -.73;
	z[6]= -.73;
	z[7]= .433;
	z[8] = .36;
	z[9] = -.7;
	z[10] = -.75;
	z[11] = -.72;
	z[12] = -.70;
	z[13] = -.75;
	z[14] = -.76;
	z[15] = -.73;
	z[16] = .41;
	z[17] = -.69;
	z[18] = .38;
	z[19] = -.76;
}

void seven(vector<double> &z)
{
	z[0] = -.55;
	z[1] = -.56;
	z[2] = .08;
	z[3] = -.54;
	z[4] = .30;
	z[5] = .09;
	z[6] = .20;
	z[7] = -.46;
	z[8] = -.65;
	z[9] = -.59;
	z[10] = -.35;
	z[11] = -.50;
	z[12] = .08;
	z[13] = .19;
	z[14] = .25;
	z[15] = .09;
	z[16] = -.34;
	z[17] = -.50;
	z[18] = .02;
	z[19] = -.42;
}

void bella()
{
  vector<double> z, conf,y,si;
	//vector<int> call;
  vector<affx::GType> call;
	vector<int> H;
	int length;
	snp_param sp;
	sp.Initialize();
	length = 20;
	si.assign(length,0);

	z.resize(length);
	y.resize(length);
	conf.resize(length);
	call.resize(length);
	z[1] = 0;
	z[2] = 0.66;
	/*for (int q=0; q<21; q++)
	{
		z[0] = -1.0 +q*0.1;
	}*/
	seven(z);
		bayes_label(call,conf,z,y,H,si, sp);
		printf("M: %f %f %f\n",sp.prior.aa.m,sp.prior.ab.m,sp.prior.bb.m);
		printf("P: %f %f %f\n",sp.posterior.aa.m,sp.posterior.ab.m,sp.posterior.bb.m);
		for (int i=0; i<length;i++)
		{	
			printf("Z: %f\t%d\t%f\n",z[i],call[i],conf[i]);
		}
		printf("\n");
	
}

void test_bins()
{
  vector<double> z, conf,y, si;
	vector<affx::GType> call;
	vector<int> H;
	int length;
	double tmp;
	snp_param sp;

	sp.Initialize();
	length = 20;

	z.resize(length);
	y.resize(length);
	call.resize(length);
	conf.resize(length);
	si.assign(length,0);

	for (int qi=0; qi<length; qi++)
	{
		tmp = .66+.001*qi;
		if (qi<2*length/3)
			tmp= 0 + .001*qi;
		if (qi<length/3)
			tmp = -.66+.001*qi;
		z[qi] = tmp;
	}
	sp.bins = 5;

		bayes_label(call,conf,z,y,H,si,sp);
		printf("M: %f %f %f\n",sp.prior.aa.m,sp.prior.ab.m,sp.prior.bb.m);
		printf("P: %f %f %f\n",sp.posterior.aa.m,sp.posterior.ab.m,sp.posterior.bb.m);
		for (int i=0; i<length;i++)
		{	
			printf("Z: %f\t%d\t%f\n",z[i],call[i],conf[i]);
		}
		printf("\n");
	sp.bins = 50;
	bayes_label(call,conf,z,y,H,si,sp);
		printf("M: %f %f %f\n",sp.prior.aa.m,sp.prior.ab.m,sp.prior.bb.m);
		printf("P: %f %f %f\n",sp.posterior.aa.m,sp.posterior.ab.m,sp.posterior.bb.m);
		for (int i=0; i<length;i++)
		{	
			printf("Z: %f\t%d\t%f\n",z[i],call[i],conf[i]);
		}
		printf("\n");	
}

int main( int argc, char *argv[])
{
	printf("hello, labeling!\n");
	//bella();
	test_bins();
	return(0);
}

