////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <Rmath.h>

#ifdef WIN32
#define snprintf _snprintf
#endif

#define MAX_LINE_LEN 100

int main(int argc, char *argv[]) {
  double x,val;
  char s[MAX_LINE_LEN];

  // Create table for pnorm
  ofstream pnormFile("pnorm.tab");
  double pnorm_test_min = -20;
  double pnorm_test_max = 20;
  int    pnorm_test_n   = 2000;
  double step           = (pnorm_test_max - pnorm_test_min)/pnorm_test_n;
  int i;
  for(i=0, x=pnorm_test_min; i<=pnorm_test_n; i++, x+=step) {
    val = pnorm(x,0,1,1,0);
    snprintf(s,sizeof(s),"%10.10f\t%20.20f\n",x,val);
    pnormFile << s;
  }
  pnormFile.close();

  // Create table for pwilcox
  ofstream pwilcoxFile("pwilcox.tab");
  int pwilcox_test_m_min = 1;
  int pwilcox_test_m_max = 7;
  int pwilcox_test_n_min = 1;
  int pwilcox_test_n_max = 5;
  int m,n,z;
  for(m=pwilcox_test_m_min; m<=pwilcox_test_m_max; m+=3) {
    int R_correction = (m*(m+1))/2;  // For some reason R defines the lowest rank to be 0 not 1, so this correction factor is required.
    for(n=pwilcox_test_n_min; n<=pwilcox_test_n_max; n++) {
      double xMax = (n+m)*(n+m-1)/2.0;
      for(z=0; z<=xMax; z++) {
        val = pwilcox(z-R_correction,m,n,1,0);
        snprintf(s,sizeof(s),"%d\t%d\t%d\t%20.20f\n",m,n,z,val);
        pwilcoxFile << s;
      }
    }
  }
  pwilcoxFile.close();

  // Create table for psignrank
  ofstream psignrankFile("psignrank.tab");
  int psignrank_test_n_min = 1;
  int psignrank_test_n_max = 15;
  for(int n=psignrank_test_n_min; n<=psignrank_test_n_max; n++) {
    // int R_correction = (n*(n+1))/2;  // For some reason R defines the lowest rank to be 0 not 1, so this correction factor is required.
    int tMax = (n*(n+1))/2;
    for(int t=0; t<=tMax; t++) {
      val = psignrank(t,n,1,0);
      snprintf(s,sizeof(s),"%d\t%d\t%20.20f\n",n,t,val);
      psignrankFile << s;
    }
  }
  psignrankFile.close();

  exit(EXIT_SUCCESS);
}
