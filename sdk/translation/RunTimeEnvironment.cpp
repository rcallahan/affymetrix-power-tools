////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   RunTimeEnvironment.cpp
 * @author Mybrid Spalding
 * @date   Wed Apr  9 12:03:14 PDT 2008
 * @brief  A single instance class to hold the global runtime environment, similar shell enviorment variables.
 */


#include "translation/RunTimeEnvironment.h"
//
#include "util/Err.h"
//
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
//

using namespace std;

/*****************************************************************************/
/**
 * RunTimeEnvironment::initializeRunTimeEnvironment
 * Synopsis:
 * Initializes the run time enviornment Verbose levels. 
 *
 */
/*****************************************************************************/
void RunTimeEnvironment::initializeRunTimeEnvironment()
{

  APT_ERR_ASSERT(!m_adtOpts.m_progName.empty(), "");

  m_programName = m_adtOpts.m_progName;

  if (m_adtOpts.m_verbosity > ADT_VERBOSE_TMI) {
    m_adtOpts.m_verbosity = ADT_VERBOSE_TMI;
  }

  if (m_adtOpts.m_verbosity >=  0) {
    m_currentVerbosity = (ADT_VERBOSE_ENUM) m_adtOpts.m_verbosity;
  } else {
    m_currentVerbosity = ADT_VERBOSE_NORMAL;
  }


}
// end RunTimeEnvironment::initializeRunTimeEnvironment
/*****************************************************************************/
/*****************************************************************************/
/**
 * RunTimeEnvironment::setVerbosity:
 * Synopsis:
 * 
 * Should not be used during normal coding. This is a debug API
 * used to pass a different verbosity level to child functions
 * that respect such things. Although in practive most verbose levels
 * are hard-coded. 
 *
 * @param level - the ADT_VERBOSE_ENUM level to set.
 *
 */
/*****************************************************************************/
void RunTimeEnvironment::setVerbosity(ADT_VERBOSE_ENUM level)
{

  m_currentVerbosity = level;

  Verbose::setLevel(level);

}
// end RunTimeEnvironment::setVerbosity
/*****************************************************************************/
/*****************************************************************************/
/**
 * RunTimeEnvironment::profilesReport
 * Synopsis:
 *
 * UNIX ONLY: this API is not intended for Windows. There is no
 * Windows check because the Profile object is coded to be empty
 * on Windows platform. This API will just be a no-op if called
 * on the Windows platform.
 *
 * Report to standard error
 * on the various dynamic timings created in m_profiles hash.
 *
 * @param totalSeconds - the total run time. 
 *
 * @return - void, prints report to stand error. 
 */
/*****************************************************************************/
void RunTimeEnvironment::profilesReport(float totalSeconds)
{

  std::map<std::string, Profile*>::iterator itSP;

  std::vector<std::string> profiles;

  for (itSP = m_profiles.begin(); itSP != m_profiles.end(); itSP++) {
    std::stringstream profileSStr;

    profileSStr << itSP->second->getElapsedFormatedString() << ": " << itSP->first;
    profileSStr << " (";

    if (totalSeconds > 0.0) {
      int percent = (int)(100.0 * (itSP->second->getElapsedSeconds() / totalSeconds));
      profileSStr << percent << "%, ";
    }
    profileSStr << itSP->second->m_totalCalls << " calls)";
    profileSStr << endl;

    profiles.push_back(profileSStr.str());
  }

  std::sort(profiles.begin(), profiles.end());

  for (int i = 0; i < profiles.size() ; i++) {

    cerr << profiles[i];
  }

  return;

}
// end RunTimeEnvironment::profilesReport
/*****************************************************************************/
