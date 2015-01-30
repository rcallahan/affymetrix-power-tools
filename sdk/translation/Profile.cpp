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
 * @file   Profile.cpp
 * @author Mybrid Spalding
 * @date   Mon Mar 24 20:18:26 PDT 2008
 * @brief  Basic profile class. Need better means use a 3rd party pkg.
 */

#include "translation/Profile.h"
//
#include "util/Err.h" // includes Verbose.h
//
#include <sstream>
#include <stdio.h>

using namespace std;

/*****************************************************************************/
/**
 * Profile::begin
 * Synopsis: Mark the begin time to be profiled. Resets endTime to startTime.
 *
 * @return - void
 */
/*****************************************************************************/
void Profile::begin()
{

#if defined (WIN32)
#else
  m_endCount = m_startCount;
  m_startCount++;

  gettimeofday(&(m_startTimeval), NULL);
  m_endTimeval = m_startTimeval;
#endif

}
// end Profile::begin
/*****************************************************************************/
/*****************************************************************************/
/**
 * Profile::end
 * Synopsis: Mark the end time of something profiled.
 *
 * @return - void, throws exception if startTime is not set.
 */
/*****************************************************************************/
void Profile::end()
{
#if defined (WIN32)
#else


  m_endCount++;

  checkTimes();

  gettimeofday(&m_endTimeval, NULL);

  m_totalElapsedSeconds +=
    ((m_endTimeval.tv_sec - m_startTimeval.tv_sec) +
     ((m_endTimeval.tv_usec - m_startTimeval.tv_usec) / 1000000.0)
    );

  m_totalCalls++;
#endif
  return;

}
// end Profile::end
/*****************************************************************************/
/*****************************************************************************/
/**
 * Profile::checkTimes
 * Synopsis: Checks if endTime - startTime > 0
 *
 * @return - void, throws exception if endTime - startTime <= 0
 */
/*****************************************************************************/
void Profile::checkTimes()
{

#if defined (WIN32)
#else
  if ((m_startCount == 0)  || (m_endCount < m_startCount)) {
    if (m_endCount < m_startCount) {
      APT_ERR_ABORT("Profile::checkTimes: end() was not called after last begin().\n");
    } else {
      APT_ERR_ABORT("Profile::checkTimes: begin() was not called.\n");
    }
  }
#endif

}
// end Profile::checkTimes()
/*****************************************************************************/
/*****************************************************************************/
/**
 * Profile::getElapsedSeconds
 * Synopsis: returns endTime - startTime
 *
 * @return - void, throws exception if startTime and endTime are not set.
 */
/*****************************************************************************/
float Profile::getElapsedSeconds()
{

  return m_totalElapsedSeconds;

}
// end Profile::getElapsedSeconds
/*****************************************************************************/
/*****************************************************************************/
/**
 * Profile::getElapsedFormatedString
 * Synopsis: returns a formated time std::string as MMmSSs.
 *
 * Assumes time will never be more than 60 minutes, or 1 hour.
 *
 * @return - void, throws exception if startTime and endTime are not set.
 */
/*****************************************************************************/
std::string Profile::getElapsedFormatedString()
{

  std::string elapsed;

#if defined (WIN32)
#else

  std::stringstream formatedSStr;
  float elapsedSeconds = getElapsedSeconds();

  int elapsedMinutes = int(elapsedSeconds / 60.0);
  if (elapsedMinutes > 0) {
    elapsedSeconds = elapsedSeconds - (elapsedMinutes * 60.0);
  }

  char elapsed_cstr[30];
  sprintf(elapsed_cstr, "%02dm%06.3fs", elapsedMinutes, elapsedSeconds);

  elapsed = elapsed_cstr;

#endif
  return elapsed;

}
// end Profile::getElapsedFormatedString
/*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////

