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

#include "dm/DmPTable.h"

float dmFGetPValue(const int NumProbePairs, const int nIndex)
{
	if(nIndex < 0 || NumProbePairs < 0)
		return -1; // wrong indecies
	else if( nIndex > ((1<<(NumProbePairs + 1)) - 1))
		return - 2;  // out of bound

	long lIndex = (1 << (NumProbePairs + 1)) - 2;
	lIndex += nIndex;
	if(lIndex < DMTABLEBOUND)
		return dmPTable[lIndex];
	else return -2; // out of bound
}

