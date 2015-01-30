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

#ifndef __PTABLE_H__
#define __PTABLE_H__

extern const float MINSIGLEVEL1[];
extern const float MINSIGLEVEL1[];
extern const float tTable[];
extern const long  TABLEBOUND;
extern const float pTable[];

float  fGetPValue(const int NumProbePairs, const int nIndex);

#endif  // end of __PTABLE_H__
