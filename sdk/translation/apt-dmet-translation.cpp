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
 * @file   apt-dmet-translation.cpp
 * @author Mybrid Spalding
 * @date   Wed Jun  4 08:29:07 PDT 2008
 * @brief  Commandline controller wrapper around TranslateEngine.run().
 */

//
#include "translation/TranslationEngine.h"
//
#include "util/AptVersionInfo.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/Util.h"
#include "util/Verbose.h"

using namespace std;


/*****************************************************************************/
/**
 * translateEngineMessageHandlingCallBack
 * Synopsis:
 *
 *  There is a catch 22 where the command line message handlers need to
 * have the PGOptions information in order to set up the message handling
 * log file name, etc. Resolve this by passing in a call back to the
 * TranslateEngine constructor which initializes the PGOptions.
 *
 * @param atdOopts - filled in with all the PGOptions data.
 *
 * @return - nada, initializes message handling.
 */
/*****************************************************************************/
void translateEngineMessageHandlingCallBack(ADTOptions & adtOpts)
{

  Err::setThrowStatus(true);

  std::string verboseLogName;

  if (adtOpts.m_logFile != "") {
    verboseLogName = adtOpts.m_logFile;
  }
  else {
    verboseLogName = Fs::join(adtOpts.m_outputDir,"apt-dmet-translation.log");
  }

  std::ofstream *verboseOFStream = new std::ofstream;

  if (!Fs::isWriteableDir(adtOpts.m_outputDir)) {
    if (Fs::mkdirPath(adtOpts.m_outputDir, false) != APT_OK) {
      APT_ERR_ABORT("Can't make or write to directory: " + adtOpts.m_outputDir);
    }
  }
  Fs::mustOpenToWrite(*(verboseOFStream), verboseLogName.c_str());

  ADT_VERBOSE_ENUM logVerbosity = ADT_VERBOSE_EXCEPTION;
  if (adtOpts.m_verbosity > ADT_VERBOSE_EXCEPTION) {
    logVerbosity = (ADT_VERBOSE_ENUM) adtOpts.m_verbosity;
  }

  LogStream *verboseLogStream = new LogStream(logVerbosity,  verboseOFStream, false);

  Verbose::pushMsgHandler(verboseLogStream);
  Verbose::pushProgressHandler(verboseLogStream);
  Verbose::pushWarnHandler(verboseLogStream);
  Verbose::setLevel(adtOpts.m_verbosity);

  return;
}
// end translateEngineMessageHandlingCallBack
/*****************************************************************************/
/*****************************************************************************/
/**
 * apt-dmet-translation main
 *
 * @param argc -
 * @param argv -
 * @return - exit 0 when ok.
 */
/*****************************************************************************/
int main(int argc, char* argv[])
{
  try {
    /// @todo This isnt how an engine should be invoked or have options set.
    ///       Attemping to process the "--help" option causes BaseEngine to print:
    ///          "running TranslationEngine" before it prints the help message.
    TranslationEngine(argv, C_CMDLINE, &translateEngineMessageHandlingCallBack).run();
    exit(0);
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
// end main
/*****************************************************************************/
