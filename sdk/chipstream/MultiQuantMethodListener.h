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
 * @file   MultiQuantMethodListener.h
 * @author Alan Williams
 * @date   Fri Jan  4 16:31:26 PST 2008
 * 
 * @brief  Class for reporting results from multiple quantification methods
 */
#ifndef _MULTIQUANTMETHODLISTENER_H_
#define _MULTIQUANTMETHODLISTENER_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/AnalysisStream.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 *   Class for reporting results of quantification methods.
 */
class MultiQuantMethodListener : public QuantMethodReport {

public:
  
  /** Constructor. */
  MultiQuantMethodListener() { }

  /** Destructor. */
  ~MultiQuantMethodListener() { }

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * 
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * 
   * @return true if success, false otherwise.
   */
  bool report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
              const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust);

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * 
   * @return true if success, false otherwise.
   */
   bool reportFailure(ProbeSetGroup &psGroup, QuantMethod &qMethod, 
                      const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust);

  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * @param qMethod - Quantification method that was used.
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod);

  // Buffer the current set of results
  std::vector<std::vector<double> > m_Results;

private:
};

#endif /* _MULTIQUANTMETHODLISTENER_H_ */
