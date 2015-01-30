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

//
#include "birdseed-dev/Prior.h"
//
#include "broadutil/BroadException.h"
//
#include <cassert>
#include <cstring>
#include <sstream>
#include <string>
//

using namespace std;
using namespace birdseed::dev;


Prior::Prior(const char *s)
{
  unsigned int temp;
  int ret = sscanf(s, "%lf %lf %lf %lf %lf %u",
                   &m_mean[0], &m_mean[1], &m_covarMatrix[0][0], &m_covarMatrix[1][0], &m_covarMatrix[1][1], &temp);
  if (ret != 6) {
    throw BroadException("Error parsing prior", __FILE__, __LINE__, s);
  }
  m_covarMatrix[0][1] = m_covarMatrix[1][0];
  m_numObservations = temp;
}

std::string Prior::tostring(std::string sep) const
{
  std::stringstream strm;
  strm << m_mean[0] << sep << m_mean[1] << sep
  << m_covarMatrix[0][0] << sep << m_covarMatrix[1][0] << sep << m_covarMatrix[1][1] << sep
  << m_numObservations;
  return strm.str();
}

Priors::Priors(uint32_t numPriors):
    numPriors(numPriors)
{
  assert(numPriors == 2 || numPriors == 3);
  memset(thePriors, '\0', sizeof(thePriors));
}

size_t Priors::numObservations() const
{
  size_t tot = 0;
  for (size_t i = 0; i < numPriors; ++i) {
    tot += thePriors[i].m_numObservations;
  }
  return tot;
}

const Prior & Priors::getDiploidPrior(DiploidIndex i) const
{
  assert(isDiploid());
  return thePriors[i];
}

const Prior & Priors::getHaploidPrior(HaploidIndex i) const
{
  assert(!isDiploid());
  return thePriors[i];
}

const Prior & Priors::getAPrior() const
{
  if (isDiploid()) {
    return getDiploidPrior(AA_INDEX);
  } else return getHaploidPrior(A_INDEX);
}

const Prior & Priors::getBPrior() const
{
  if (isDiploid()) {
    return getDiploidPrior(BB_INDEX);
  } else return getHaploidPrior(B_INDEX);
}

const Prior & Priors::getPrior(size_t i) const
{
  assert(i < (sizeof(thePriors) / sizeof(thePriors[0])));
  return thePriors[i];
}

Prior & Priors::getPrior(size_t i)
{
  assert(i < (sizeof(thePriors) / sizeof(thePriors[0])));
  return thePriors[i];
}

Prior & Priors::getDiploidPrior(DiploidIndex i)
{
  assert(isDiploid());
  return thePriors[i];
}

Prior &Priors::getHaploidPrior(HaploidIndex i)
{
  assert(!isDiploid());
  assert(i < (sizeof(thePriors) / sizeof(thePriors[0])));
  return thePriors[i];
}

void Priors::setPrior(const Prior &prior, size_t i)
{
  assert(i < (sizeof(thePriors) / sizeof(thePriors[0])));
  thePriors[i] = prior;
}
void Priors::setNumPriors(size_t num)
{
  numPriors = num;
}

/******************************************************************/
/**************************[END OF Prior.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
