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
 * @file   Profile.h
 * @author Mybrid Spalding
 * @date   Mon Mar 24 20:18:26 PDT 2008
 * @brief  Basic profile class. Need better means use a 3rd party pkg.
 */

#ifndef TRANSLATION_PROFILE_H
#define TRANSLATION_PROFILE_H

#include <cstring>
#include <string>
//
#if defined (WIN32)
#else
#include "sys/time.h"
#endif

#define PROFILE_BEGIN( _p ) \
  { \
    p.begin(); \
  }

#define PROFILE_END_VERBOSE( _p, _level, _msg )                         \
  {                                                                     \
    _p.end();                                                           \
    Verbose::out(_level, ToStr(__FILE__ ) + ToStr(":") + ToStr(__LINE__) \
                 + ToStr(": ") +  ToStr(_msg) + ToStr(": ") +           \
                 p.getElapsedFormatedString() );                        \
  }

/** Place to group program options. */
class Profile
{

public:

  float   m_totalElapsedSeconds;
  int     m_totalCalls;
#if defined (WIN32)
#else
  timeval m_startTimeval;
  timeval m_endTimeval;
#endif
  int     m_startCount;
  int     m_endCount;

  Profile() {
    m_totalCalls = 0;
    m_totalElapsedSeconds = 0.0;
    m_startCount = 0;
    m_endCount   = 0;
  }

  void          begin();
  void          checkTimes();
  void          end();
  float         getElapsedSeconds();
  std::string   getElapsedFormatedString();

};


#endif /* TRANSLATION_PROFILE_H */
