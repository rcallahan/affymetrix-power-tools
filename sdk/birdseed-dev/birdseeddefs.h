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
 * FILE birdseeddefs.h
 */

#ifndef _BIRDSEEDDEFS_H_
#define _BIRDSEEDDEFS_H_

//
#include "birdseed-dev/Matrix.h"
//
#include "chipstream/BioTypes.h"
//
#include <cstddef>
//

namespace birdseed 
{
namespace dev
{
    
    static const size_t NUM_ALLELES = 2;
    static const size_t MAX_NUM_CLUSTERS = 3;

    static const double pi = 3.1415926535897;

    // Can now handle a single sample.
    static const size_t MIN_SAMPLES_TO_CLUSTER = 1;

    enum {A_ALLELE_INDEX = 0,
          B_ALLELE_INDEX = 1
    };

    typedef affx::Gender Gender;

    static const Gender GENDER_FEMALE = affx::Female;
    static const Gender GENDER_MALE = affx::Male;
    static const Gender GENDER_UNKNOWN = affx::UnknownGender;

    // Caution -- only works for arrays, not for pointers!!
#define ARRAY_SIZE(ary) (sizeof(ary)/sizeof(ary[0]))

    typedef FixedMatrix<double, NUM_ALLELES, NUM_ALLELES> Matrix2x2;

    typedef VarMatrix<double, NUM_ALLELES> IntensityMatrix;

    typedef VarVarMatrix<double> PXI_ZJMatrix;
 
};
};

#endif /* _BIRDSEEDDEFS_H */

/******************************************************************/
/**************************[END OF birdseeddefs.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
