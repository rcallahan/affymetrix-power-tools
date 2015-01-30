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
 * @file   MultiQuantMethodListener.cpp
 * @author Alan Williams
 * @date   Fri Jan  4 16:31:26 PST 2008
 *
 * @brief  Class for reporting results from multiple quantification methods
 */

//
#include "chipstream/MultiQuantMethodListener.h"
//
#include <cstdlib>
#include <iostream>
//

using namespace std;

/**
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 *
 * @param qMethod - Quantification method to be used.
 *
 * @return true if success, false otherwise.
 */
bool MultiQuantMethodListener::prepare(QuantMethod &qMethod, const IntensityMart &iMart)
{
  return true;
}

/**
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 *
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 *
 * @return true if success, false otherwise.
 */
bool MultiQuantMethodListener::report(ProbeSetGroup &psGroup,
                                      QuantMethod &qMethod,
                                      const IntensityMart &iMart,
                                      std::vector<ChipStream *> &iTrans,
                                      PmAdjuster &pmAdjust)
{
  QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
  if(eMethod == NULL) {
    Err::errAbort("Can only use a QuantMethodExprReport with a QuantExprMethod.");
  }
  APT_ERR_ASSERT(!psGroup.probeSets.empty(),"Cant have and empty psGroup.");

  // We only report expression probe sets.
  if(psGroup.probeSets[0]->psType != ProbeSet::Expression) {
    return false;
  }
  
  // Write signal entry to buffer writer.
  vector<double> results;
  for (int chip=0; chip<eMethod->getNumTargets(); chip++) {
    results.push_back(eMethod->getSignalEstimate(chip));
  }
  m_Results.push_back(results);
  
  return true;
}

/**
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 *
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 *
 * @return true if success, false otherwise.
 */
bool MultiQuantMethodListener::reportFailure(ProbeSetGroup &psGroup,
                                             QuantMethod &qMethod,
                                             const IntensityMart &iMart,
                                             std::vector<ChipStream *> &iTrans,
                                             PmAdjuster &pmAdjust)
{
    return true;
}

/**
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * @param qMethod - Quantification method that was used.
 * @return true if success, false otherwise.
 */
bool MultiQuantMethodListener::finish(QuantMethod &qMethod)
{
    return true;
}
