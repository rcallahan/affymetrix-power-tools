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
#pragma once
#include "calvin_files/portability/src/AffymetrixBaseTypes.h"
#include <string>
#include <vector>

using namespace std;

class SnpSummaryStats
{
public:
	SnpSummaryStats(void);
	~SnpSummaryStats(void);

	void CalculateMetrics(vector<vector<u_int8_t> > &calls, int numSets, vector<vector<float> > &metrics);

private:
	float TwoDigitPrecision(float f);

};
