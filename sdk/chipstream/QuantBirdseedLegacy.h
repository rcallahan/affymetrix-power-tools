////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/*  
 * FILE QuantBirdseedLegacy.h
 */

#ifndef _QUANTBIRDSEEDLEGACY_H_
#define _QUANTBIRDSEEDLEGACY_H_


//
#include "chipstream/BioTypes.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantBirdseed.h"
#include "chipstream/QuantBirdseedv1.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodExprReport.h"
//
#include "birdseed-dev/PriorsReader.h"
#include "birdseed-v1/FitSNPGaussiansPriors3.h"
#include "birdseed-v1/GenotypeCaller.h"
#include "stats/stats.h"
//
#include "newmat.h"
//
#include <cfloat>
#include <cstring>
#include <string>
#include <vector>
//


/// String describing quant method
#define QUANTBIRDSEEDLEGACY "birdseed"

class QuantBirdseedLegacy : public QuantBirdseedv1 {

  public:

    /** Constructor, currently creates prior estimates from pre-calculated data from R. */
  QuantBirdseedLegacy(double confidenceThreshold, double correctionFactor);

  ~QuantBirdseedLegacy();
    
    /** Fill in our self documentation from current state. */
  static void setupSelfDoc(SelfDoc &doc);

    /** 
     * @brief What is the name of the quantification method?
     * @return name of adjuster.
     */  
  std::string getType();

  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  static SelfDoc explainSelf();

};

#endif /* _QUANTBIRDSEEDLEGACY_H_ */

/******************************************************************/
/**************************[END OF QuantBirdseed.h]****************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */

