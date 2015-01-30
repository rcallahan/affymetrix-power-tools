////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "copynumber/Indexer.h"
//

using namespace std;

// Conversion of probe IDs to row and column indices and the inverse.
// Note that the probe ID is assumed to be from a cyto2 array because
// of the way that probe IDs are mapped to positions.  When looking at
// the image with the logo in the top left corner, rows start at the top
// and columns start at the left.  That is, a right hand coordinate
// system.  Probe IDs start at 1 in the top right corner of the image
// then go downward to the bottom and proceed into the next column
// just to the left.

// probe ID to row colunn.  row is .first, column is .second. 
pair<int,int> index_to_rc(int probe_id, int ncol)
	{
	int pos = probe_id - 1;
	return make_pair<int,int>(pos % ncol, ncol - pos/ncol - 1);
	}

// row and column to probe ID.
int rc_to_index(int row, int col, int ncol)
	{
	return ncol*(ncol - col - 1) + row + 1;
	}

pair<int,int> probeid_to_rc(int probe_id, int nrow, int ncol)
	{
	int row = (probe_id - 1) % nrow;
	int col = ncol - (probe_id - 1)/nrow - 1;
	return make_pair<int,int>(row,col);
	}

