////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#ifndef BINARYSEGMENTER_H_
#define BINARYSEGMENTER_H_

#include "copynumber/BinaryCopyNumberNode.h"

#include <vector>


void binary_segmenter(std::vector<BinaryCNode>::iterator bcn_start,
	std::vector<BinaryCNode>::iterator bcn_end,
	std::vector<BinaryCNodeParam> bcn_params,
	int level=0);


void binary_segmenter_normal(std::vector<BinaryCNode>::iterator bcn_start,
	std::vector<BinaryCNode>::iterator bcn_end,
	std::vector<BinaryCNodeParamNormal> bcn_params,
	int level=0);


std::vector<BinaryCNodeParam> make_bcn_params(std::valarray<double> mu,
		double prob0, double tprob);


std::vector<BinaryCNodeParamNormal> make_bcn_params_normal(
	std::valarray<double> mu, std::valarray<double> var,
	std::valarray<double> tprob, std::valarray<double> bias);

#endif
