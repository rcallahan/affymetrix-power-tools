////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
 * @file   TrustedProbesCheck.h
 * @brief  Extend MatrixCheck to allow nested columns
 * where the outer columns are tab delimited. 
 */
#ifndef _TRUSTED_PROBES_CHECK_H_
#define _TRUSTED_PROBES_CHECK_H_

#include "util/MatrixCheck.h"
//
#include <string>
//

struct TrustedProbesCheck {
  std::string name;
  bool posteriorMatrixCheck;
  bool matrixCheck;
  int skiplines;
  int skipcols;
  double eps;
    
};


class PosteriorMatrixCheck : public MatrixCheck
{
public:
  
  PosteriorMatrixCheck(
                const std::string &generated, 
                const std::string &gold, 
                double eps, 
                int rowSkip, 
                int colSkip, 
                bool matchNames, 
                unsigned int allowedMisMatch, 
        	double frac = 0.0 ) :
  MatrixCheck(
                generated, 
                gold, 
                eps, 
                rowSkip, 
                colSkip, 
                matchNames, 
                allowedMisMatch, 
		frac ) {}
  
  // replace the commas with tabs.
  bool check(std::string &msg);
};


#endif /* _TRUSTED_PROBES_CHECK_H_ */
