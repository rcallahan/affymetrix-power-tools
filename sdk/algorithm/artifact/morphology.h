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

/**
 * @file   Morphology.h
 * @author Earl Hubbell
 * @date   June 4 2009
 * 
 */

#ifndef _MORPHOLOGY_H 
#define _MORPHOLOGY_H

//
#include <vector>

void dilate(std::vector<int> &BlemishMap, int Xdim, int Ydim, int k);
void erode(std::vector <int> &BlemishMap, int Xdim, int Ydim, int k);
void dc_block(std::vector<float> &res, int Xdim, int Ydim, int k);

#endif /* _MORPHOLOGY_H_ */

