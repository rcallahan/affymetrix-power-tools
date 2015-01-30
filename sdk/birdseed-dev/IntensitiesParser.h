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

#ifndef _INTENSITIESPARSER_H_
#define _INTENSITIESPARSER_H_

//
#include "birdseed-dev/birdseeddefs.h"
//
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
//

namespace birdseed 
{
namespace dev
{
    // Read a file of intensities, calculate the correction factor, and then
    // return the IntensityMatrix for each SNP in the input intensities file.
    class IntensitiesParser
    {
      private:
        std::string dataPath;
        FILE *fp;
        double correctionFactor;
        std::string header;
        std::string snpName;
        std::auto_ptr<IntensityMatrix> intensities;
    
    
      public:
      IntensitiesParser(const std::string& path, bool calculateCorrectionFactor);
      ~IntensitiesParser();
      double getCorrectionFactor() const
        {
            return correctionFactor;
        }
    
        bool advanceSNP();
        const std::string &getHeader() const
        {
            return header;
        }
    
        const std::string &getCurrentSNPName() const
        {
            return snpName;
        }
    
        const IntensityMatrix &getCurrentIntensities() const
        {
            return *(intensities.get());
        }
    };
 
};
};


#endif /* _INTENSITIESPARSER_H_ */

/******************************************************************/
/**************************[END OF IntensitiesParser.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
