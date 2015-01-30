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

// this routine holds image processing routines that handle some blemish-related
// morphological operations

//
#include "algorithm/artifact/morphology.h"
//
#include <algorithm>

using namespace std;

//////////

inline long xhash(int i,int j,int Xdim,int Ydim){
  return(j*Xdim+i);
}

void manhattan_transfer(std::vector<int> &BlemishMap, int Xdim, int Ydim){
	int MD=Xdim+Ydim;
	long tk;
    // northwest->southeast
	// array consists only of zeros and 1s to start with
	int i,j;

	//Take out the if/else in the inner loop because they're slow
	//Use the hash implicitly in address
	j=0;
	i=0;
	tk=0;
	// northwest corner is where we start, has no neighbors
	if (BlemishMap[tk]==1)
		BlemishMap[tk]=0;
	else
		BlemishMap[tk] = MD;
	// north row has no north neighbors
	for (i=1; i<Xdim; i++){
		tk=xhash(i,j,Xdim,Ydim);
		if (BlemishMap[tk]==1){
			BlemishMap[tk] =0;
		}else{
			BlemishMap[tk] = MD;
			BlemishMap[tk] = min(BlemishMap[tk],BlemishMap[tk-1]+1);
		}
	}

	// west column has no west neigbors
	i=0; 
	for (j=1; j<Ydim; j++){
		tk = xhash(i,j,Xdim,Ydim);
		if (BlemishMap[tk]==1){
			BlemishMap[tk]=0;
		} else {
			BlemishMap[tk] = MD;
			BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk-Xdim]+1);
		}
	}
	

    for (i=1; i<Xdim; i++){
        for (j=1; j<Ydim; j++){
		tk=xhash(i,j,Xdim,Ydim);
            if (BlemishMap[tk] == 1){
                BlemishMap[tk] = 0;
            } else {
                BlemishMap[tk] = MD;
		// westward neighbor
                BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk-1]+1);
		// northward neighbor
                BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk-Xdim]+1);
            }
        }
    }
    // southeast->northwest
// array consists of distances to nearest blemish in north or west

    // same unrolling to take out if/else switch in the loop
    i=Xdim-1;
    j=Ydim-1;
    // southeast corner has no neighbors, so nothing to be done
    // Proceed along south row looking at east neighbor
    for (i=Xdim-2; i>=0; i--){
	    tk = xhash(i,j,Xdim,Ydim);
	    // east neighbor
	    BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk+1]+1);
    }
    // proceed along east column looking at south neighbor
    i = Xdim-1;
    for (j=Ydim-2; j>=0; j--){
	    tk = xhash(i,j,Xdim,Ydim);
	    BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk+Xdim]+1);
	}

    for (i=Xdim-2; i>=0; i--){
        for (j=Ydim-2; j>=0; j--){
		tk = xhash(i,j,Xdim,Ydim);
		// where is the closest blemish?
            // either what we had on the first pass
            // east neighbor
		BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk+1]+1);
            // west neighbor
		BlemishMap[tk] = min(BlemishMap[tk], BlemishMap[tk+Xdim]+1);
        }
    }
}

// invert below/above threshold
void invert(std::vector<int> &BlemishMap,int k){
	for (long tk=0; tk<BlemishMap.size(); tk++){
            BlemishMap[tk] = ((BlemishMap[tk]<=k)?1:0);
        }
   }

// dilate the BlemishMap by k efficiently
void dilate(std::vector<int> &BlemishMap, int Xdim, int Ydim, int k){
    manhattan_transfer(BlemishMap,Xdim,Ydim);
    invert(BlemishMap,k);
}

// erode the BlemishMap by k by inverting the map
void erode(std::vector <int> &BlemishMap, int Xdim, int Ydim, int k){
	invert(BlemishMap,0);
	dilate(BlemishMap,Xdim,Ydim,k);
	invert(BlemishMap,0);	
}

/// put this here for now
//
void dc_block(std::vector<float> &res, int Xdim, int Ydim, int k){
	// kxk block impulse filtering
	
	vector<double> tmpc,tmpcr; 

	tmpc.assign(res.size(),0);
	tmpcr.assign(res.size(),0);
	int i,j, curpos, rpos, xb,yb,kx;
	double impulse;
	kx = 2*k+1;

	// logic does not work for k=0
	// running sum down columns
	for (i=0; i<Xdim; i++){
		curpos = i;
		tmpc[curpos] = res[curpos];
		for (j=1; j<kx; j++){
			curpos = curpos + Xdim;
			tmpc[curpos] = res[curpos] + tmpc[curpos-Xdim];
		}
		for (j=kx; j<Ydim; j++){
			curpos = curpos + Xdim;
			tmpc[curpos] = res[curpos] + tmpc[curpos-Xdim]-res[curpos-kx*Xdim];
		}
	}
	// running sum across rows
	for (j=0; j<Ydim; j++){
		curpos = Xdim*j;
		tmpc[curpos] = tmpc[curpos];
		for (i=1; i<kx; i++){
			curpos = curpos+1;
			tmpcr[curpos] = tmpc[curpos] + tmpcr[curpos-1];
		}
		for (i=kx; i<Xdim; i++){
			curpos = curpos+1;
			tmpcr[curpos] = tmpc[curpos] + tmpcr[curpos-1] - tmpc[curpos-kx];
		}
	}
	// now tmpcr has 2k+1x2k+1 block information in running sums
	int den = kx*kx; // we'll be using this a lot

	// dc_impulse = kxk sum - current item / (k*k-1)
	for (j=0; j<Ydim; j++){
		yb = min(max(j-k,0)+2*k,Ydim-1);
		for (i=0; i<Xdim; i++){
			xb = min(max(i-k,0)+2*k,Xdim-1);
			rpos = xhash(i,j,Xdim,Ydim);
			// running sum minus div item = mean residual
			impulse = (tmpcr[xhash(xb,yb,Xdim,Ydim)]) / den;
			// local residual average, to be divided out - exactly the res if k=0
			res[rpos] =  impulse;
		}
	}
}
