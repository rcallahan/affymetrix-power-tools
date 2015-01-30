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

#include "SnpSummaryStats.h"
#include "stats/statfun.h"
#include "file/CHPFileData.h"

#include <string>
#include <vector>

using namespace std;
using namespace affxchp;

SnpSummaryStats::SnpSummaryStats(void)
{
}

SnpSummaryStats::~SnpSummaryStats(void)
{
}


/**
 * Convert the flow to two digits of precision
 * @param f The float value
 * @return The float with only two digits of precision
 */
float SnpSummaryStats::TwoDigitPrecision(float f)
{
    return (float)((floor((f+0.005)*100))/100.0);
}


/*
 * Compute the metrics
 */
void SnpSummaryStats::CalculateMetrics(vector<vector<u_int8_t> > &calls, int numSets, vector<vector<float> > &metrics)
{
    // Loop over each probe set, calculate and store the results.
    for (int iset=0; iset<numSets; iset++)
    {
        int aa=0;
        int ab=0;
        int bb=0;
        int nc=0;
        int nfiles = (int)calls[iset].size();
        for (int ifile=0; ifile<nfiles; ifile++)
        {
            u_int8_t call = calls[iset][ifile];
            if (call == ALLELE_A_CALL)
                ++aa;
            else if (call == ALLELE_AB_CALL)
                ++ab;
            else if (call == ALLELE_B_CALL)
                ++bb;
            else
                ++nc;
        }
        int index=0;
        metrics[iset][index++] = TwoDigitPrecision(100.0f * (aa+ab+bb) / nfiles);
        metrics[iset][index++] = 100.0f * aa / nfiles;
        metrics[iset][index++] = 100.0f * ab / nfiles;
        metrics[iset][index++] = 100.0f * bb / nfiles;

		float freqA = (float)(2*aa + ab) / (2*aa + 2*ab + 2*bb);
		float freqB = (float)(2*bb + ab) / (2*aa + 2*ab + 2*bb);
        metrics[iset][index++] = min(freqA, freqB);
        metrics[iset][index++] = (float)affxstat::CalcHWEqPValue(aa, ab, bb);
    }
}

