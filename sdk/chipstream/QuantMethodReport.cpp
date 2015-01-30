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
 * @file   QuantMethodReport.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 12:09:42 2005
 * 
 * @brief  Class for reporting results of quantification methods.
 */

#include "chipstream/QuantMethodReport.h"
//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/TsvReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 * Virtual destructor for a virtual class.
 */
QuantMethodReport::~QuantMethodReport() {
};

/** 
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodReport::reportFailure(ProbeSetGroup &psGroup, 
                                      QuantMethod &qMethod,
                                      const IntensityMart &iMart, 
                                      std::vector<ChipStream *> &iTrans, 
                                      PmAdjuster &pmAdjust) {
  return true;
}


void QuantMethodReport::addStdHeaders(QuantMethodReport *qReport,
                                              const std::string& execGuid, 
                                              const std::string& reportGuid,
                                              const std::string& timeStr,
                                              const std::string& commandLine,
                                              const std::string& execVersion,
                                              const AnalysisInfo& info)
{
  std::vector<std::string>::const_iterator keyIx, paramIx;
  //
  for(keyIx = info.m_ParamNames.begin(), paramIx = info.m_ParamValues.begin();
      keyIx != info.m_ParamNames.end() && paramIx != info.m_ParamValues.end();
      ++keyIx, ++paramIx) {
    // @todo should we using the AffymetrixParameterConsts.h #defined values?
    qReport->addHeader("affymetrix-algorithm-param-" + *keyIx, *paramIx);
  }
  for(keyIx = info.m_ClientInfoNames.begin(), paramIx = info.m_ClientInfoValues.begin();
      keyIx != info.m_ClientInfoNames.end() && paramIx != info.m_ClientInfoValues.end();
      ++keyIx, ++paramIx) {
    // @todo should we using the AffymetrixParameterConsts.h #defined values?
    qReport->addHeader("affymetrix-application-meta-data-info-" + *keyIx, *paramIx);
  }
}
