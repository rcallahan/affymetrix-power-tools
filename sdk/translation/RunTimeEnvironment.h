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
 * @file   RunTimeEnvironment.h
 * @author Mybrid Spalding
 * @date   Thu Mar 27 09:24:26 PDT 2008
 * @brief  A single instance class to hold the global runtime environment, similar shell enviorment variables.
 */


#ifndef TRANSLATION_RUN_TIME_ENVIRONMENT_H
#define TRANSLATION_RUN_TIME_ENVIRONMENT_H

#include "translation/ADTOptions.h"
#include "translation/Profile.h"
//
#include <cstring>
#include <iostream>
#include <map>
#include <string>
//

const std::string ADT_VERSION = "1.0.2";

const std::string TRANSLATION_CVS_ID =
  "$Id: RunTimeEnvironment.h,v 1.25 2009-09-30 22:33:24 csugne Exp $";

const std::string TEST_DATA_ROOT_DIR =
  "../../regression-data/data/idata/translation";

const std::string TEST_DATA_UNIT_DIR  =
  "./input/";

const std::string TEST_DATA_REGRESSION_DIR  =
  TEST_DATA_ROOT_DIR  + "/regression-data";
  
#define MARK() cerr << "***MARK***" << __FILE__ << " : " << __LINE__ << " : " << endl;

// A run-time mask is used so as to enable multipe types.
const unsigned int ADT_INVALID_OPTION_MASK_MISSING_TRANSLATION_TABLE = 0x1;
const unsigned int ADT_INVALID_OPTION_MASK_MISSING_EXPERIMENT_CHP    = 0x2;

enum ADT_VERBOSE_ENUM {
  // ADT_VERBOSE_NULL, 0: quiet, don't output anything. 
  ADT_VERBOSE_NULL,        
  // ADT_VERBOSE_NORMAL, 1: the only option respected by the Windows Console
  // This means "exceptions" have to be written to this level of verbose.
  // The "exception" verbose level is only respected by the command line tools.
  ADT_VERBOSE_NORMAL,
  // ADT_VERBOSE_EXCEPTION, 2: command line only. When used the log messages
  // will report to the log file but not the command-line console. Note
  // that this behavior is not respected by the Windows console.
  // Exceptions for the Windows console must be written ADT_VERBOSE_NORMAL.
  ADT_VERBOSE_EXCEPTION,
  // ADT_VERBOSE_INPUT_FILES, 3: Outputs information related to input files
  // that may help with understanding calling. The most common log
  // message through out the application. 
  ADT_VERBOSE_INPUT_FILES,
  // ADT_VERBOSE_CALL, 4: Specific to GeneCall and the calling algorithm.
  // Use this to see why a call was made. 
  ADT_VERBOSE_CALL,
  // ADT_VERBOSE_TMI, 5: Too much information level is used for
  // output of common problem touch points in the code base.
  // Any verbose level specified > 5 is reduced to 5. There is no level
  // greater than TMI. 
  ADT_VERBOSE_TMI,
};

// Environment
class RunTimeEnvironment
{

public:

  // Program
  std::string        m_programName;

  unsigned int       m_invalidOptionsMask;

  // Verbose.h
  std::string        m_verboseLogName;
  ADT_VERBOSE_ENUM   m_currentVerbosity;

  //Profiling
  Profile            m_programTime;

  std::map< std::string, Profile* > m_profiles;

  // Options
  ADTOptions        m_adtOpts;

  RunTimeEnvironment() {
    m_currentVerbosity = ADT_VERBOSE_NORMAL;
  };

  ~RunTimeEnvironment() {};

  void initializeRunTimeEnvironment();

  void setVerbosity(ADT_VERBOSE_ENUM level);

  void profilesReport(float totalSeconds = 0);

};


#endif /* TRANSLATION_RUN_TIME_ENVIRONMENT_H */
