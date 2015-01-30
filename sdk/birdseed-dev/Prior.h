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

#ifndef _PRIOR_H_
#define _PRIOR_H_

//
#include "birdseed-dev/Matrix.h"
#include "birdseed-dev/birdseeddefs.h"
//
#include "file/CELFileData.h" // this defines STRUCT_ALIGNMENT
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <string>
//

namespace birdseed
{
namespace dev {
// Holds a prior for a single cluster of a single SNP.
class Prior
{
public:

#ifdef _MSC_VER
#pragma pack(push, 1)
#endif
#ifdef __APPLE__
#pragma options align=packed
#endif

  // We cant have STRUCT_ALIGNMENT defined as that will cause
  // the doubles to be unaligned causing a bus error.
  // Here we undefine it so things work, but we cant read Binary format files.
#ifdef __sparc__
#undef STRUCT_ALIGNMENT
#define STRUCT_ALIGNMENT
#endif

  FixedVector<double, NUM_ALLELES> m_mean STRUCT_ALIGNMENT;

  FixedMatrix<double, NUM_ALLELES, NUM_ALLELES>  m_covarMatrix STRUCT_ALIGNMENT;

  uint32_t m_numObservations STRUCT_ALIGNMENT;

#ifdef _MSC_VER
#pragma pack(pop)
#endif
#ifdef __APPLE__
#pragma options align=reset
#endif

  Prior(const char *s);
  Prior()    {}

  std::string tostring(std::string sep = " ") const;

};

// Holds all the priors for a single SNP.
class Priors
{
public:
  enum DiploidIndex {AA_INDEX = 0,
                     AB_INDEX = 1,
                     BB_INDEX = 2
                    };

  enum HaploidIndex {A_INDEX = 0,
                     B_INDEX = 1
                    };


  Priors(uint32_t numPriors);

  size_t numObservations() const;

  uint32_t getNumPriors() const  {
    return numPriors;
  }

  bool isDiploid() const  {
    return numPriors == 3;
  }

  const Prior &getPrior(size_t i) const;

  const Prior &getDiploidPrior(DiploidIndex i) const;

  const Prior &getHaploidPrior(HaploidIndex i) const;

  const Prior &getAPrior() const;

  const Prior &getBPrior() const;

  Prior &getPrior(size_t i);

  Prior &getDiploidPrior(DiploidIndex i);

  Prior &getHaploidPrior(HaploidIndex i);

  void setPrior(const Prior &prior, size_t i);

  void setNumPriors(size_t num);

private:

#ifdef _MSC_VER
#pragma pack(push, 1)
#endif
#ifdef __APPLE__
#pragma options align=packed
#endif

  uint32_t numPriors STRUCT_ALIGNMENT;
  Prior thePriors[MAX_NUM_CLUSTERS];

#ifdef _MSC_VER
#pragma pack(pop)
#endif
#ifdef __APPLE__
#pragma options align=reset
#endif

};

};
};

#endif /* _PRIOR_H_ */

/******************************************************************/
/**************************[END OF Prior.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
