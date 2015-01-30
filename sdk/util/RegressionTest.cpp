////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

#include "util/RegressionTest.h"
//
#include "util/Fs.h"
#include "util/MatrixCheck.h"
#include "util/MixedFileCheck.h"
#include "util/RowFile.h"
#include "util/SQLite.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
//
using namespace std;

RegressionTest::RegressionTest() {
  initialize();
}

RegressionTest::RegressionTest(const std::string& name,
                               const std::string& generatedFile, 
                               const std::string& goldFile,
                               double epsilon, 
                               const std::string& command,
                               int rowSkip,
                               int colSkip,
                               bool matchNames,
                               int allowedMisMatch,
                               bool negTest) {
  initialize();
  m_Name = name;
  m_Checks.push_back(new MatrixCheck(generatedFile, goldFile, epsilon,
                                     rowSkip, colSkip, matchNames, allowedMisMatch));
  m_CommandStr = command;
  m_CommandStr = Fs::convertCommandToUnc(m_CommandStr);
  m_Delete = true;
  m_NegTest = negTest;
  m_64bit = !Util::is32Bit();
}

RegressionTest::RegressionTest(	const std::string& name,
                                const std::string& command, 
                                const std::vector<RegressionCheck *>& fileChecks,
                                bool doDelete,
                                bool negTest) {
  initialize();
  m_Name = name;
  m_CommandStr = command;
  m_CommandStr = Fs::convertCommandToUnc(m_CommandStr);
  m_Checks = fileChecks;
  m_Delete = doDelete;
  m_NegTest = negTest;
  m_64bit = !Util::is32Bit();
}

void RegressionTest::initialize() {
  m_SQL = NULL;
  m_RunTime = DATA_NA;
  m_MaxMem = DATA_NA;
  m_MemOk = DATA_NA;
  m_MemLeaks = DATA_NA;
  m_NumPassed = DATA_NA;
  m_NumFailed = DATA_NA;
  m_64bit = false;
}

RegressionTest::~RegressionTest() {
  if (m_Delete) {
    for(uint32_t i = 0; i < m_Checks.size(); i++) {
      delete m_Checks[i];
    }
  }
  if (m_SQL != NULL) {
    m_SQL->close();
    Freez(m_SQL);
  }
}

bool RegressionTest::run() {
  bool success = true;
  string toRun = m_CommandStr;

  assert(m_CommandStr != "");
  if (!m_ValgrindLog.empty()) {
    toRun = "valgrind --tool=memcheck --num-callers=20 --leak-check=full --log-file=" + m_ValgrindLog + " " + m_CommandStr;
  }  
  
  printf("====\n");
  printf("==== Test: '%s'\n",m_Name.c_str());
  printf("==== Cmd:  '%s'\n",toRun.c_str());
  printf("====\n");

  time_t startTime = time(NULL);
  int retVal = system(toRun.c_str());
  time_t endTime = time(NULL);
  m_RunTime = endTime - startTime; // convert to ms
  if(retVal == -1) {
    m_Error += ToStr(" Failed to run command: ") + toRun;
      success = false;
  }
  if(retVal != 0) {
    m_Error += ToStr(" Command exited with non-zero status: ") + toRun;
    success = false;
  }
  return success;
}

void RegressionTest::setDatabase( const std::string & path ) {
  
  if ( path.empty() || !Fs::fileExists(path) ) {
    if ( !path.empty() ) {
      Verbose::warn(1, std::string("RegressionTest ignoring database file,  file not found: ") + path);
    }
    return;
  }
  if ( !Fs::isWriteable(path) ) {
    Verbose::warn(1, std::string("RegressionTest ignoring database file,  permission denied: ") + path);
    return;
  }
  if ( !Fs::isWriteableDir( Fs::dirname(path) ) ) {
    Verbose::warn(1, std::string("RegressionTest ignoring database file, parent directory write permission denied and required for sqlite3 journal file: ") + Fs::dirname(path) );
    return;
  }
  
  m_Database = path;
}

bool RegressionTest::pass() {
  bool success = true;
  std::string xmlFile = "JUnitTestResults." + m_Name + ".xml";
  std::ofstream xmlOut;
  std::vector<std::string> m_Results(m_Checks.size() + 1);

  if(Fs::fileExists(xmlFile))
    Fs::rm(xmlFile, false);

  m_NumFailed = 0;
  m_NumPassed = 0;
  int numSkipped = 0;
  m_TimeStamp = time((time_t*)NULL);
  if(run()) {
    m_NumPassed++;
    if (!m_LogToParse.empty()) {
      //parseLogFile(m_LogToParse);
    }
    if (!m_ValgrindLog.empty()) {
      parseValgrindFile(m_ValgrindLog);
    }
    m_Results[0] = "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name +".run\" time=\""+ToStr(m_RunTime)+"\"/>\n";
    for(unsigned int i = 0; i < m_Checks.size(); i++) {
      std::string msg;
      bool currentPass;
      try {
        currentPass = m_Checks[i]->check(msg);
        if(m_Checks[i]->m_NegTest) {
          currentPass = !currentPass;
          msg = "test passed, but was not expected to pass";
        }
      } 
      catch(...) {
        currentPass = false;
        msg = "unexpected exception thrown";
      }
      if(!currentPass) {
        m_Error += msg;
        success = false;
        m_NumFailed++;
        m_Results[i+1] = 
          "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name + "." + m_Checks[i]->m_Name +"\">\n" +
          "      <error message=\"check failed\">\n" + msg +
          "      </error>\n" +
          "  </testcase>\n";
      } else {
        m_NumPassed++;
        m_Results[i+1] = "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name + "." + m_Checks[i]->m_Name +"\"/>\n";
      }
    }
    if (!m_Database.empty()) {
      recordResults(success);
    }
  }
  /* failed just running the test. */
  else {
    std::string reportMsg;
    numSkipped = m_Checks.size();
    if(m_NegTest) {
      m_NumFailed = 0;
      m_NumPassed = 1;
      reportMsg = "test skipped because of expected execution failure -- this test should probably not exist";
      success = true;
      m_Results[0] = "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name +".run\"/>\n";
    } else {
      m_NumFailed = 1;
      m_NumPassed = 0;
      reportMsg = "test skipped because of unexpected execution failure";
      success = false;
      m_Results[0] = "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name +".run\">\n" +
        "      <error message=\"execution failed\">\n" + m_Error + "\n" +
        "      </error>\n" +
        "  </testcase>\n";
    }
    if (!m_Database.empty()) {
      recordResults(success);
    }
    for(size_t i = 0; i < m_Checks.size(); i++) {
      m_Results[i+1] = 
        "  <testcase classname=\"" + m_Name + "\" name=\""+ m_Name + "." + m_Checks[i]->m_Name +"\">\n" +
        "      <skipped message=\"" + reportMsg + "\"/>\n" + 
        "  </testcase>\n";
    }

  }

  // Generate xml Junit file
  xmlOut.open(xmlFile.c_str());
  xmlOut << "<?xml version=\"1.0\"?>\n";
  xmlOut << "<testsuite errors=\"0\" skipped=\"" + ToStr(numSkipped) + "\" failures=\"" + ToStr(m_NumFailed) + "\" tests=\"" + ToStr(m_Checks.size() + 1) + "\" name=\"" + m_Name + "\">\n";
  for(size_t i=0; i< m_Results.size(); i++)
    xmlOut << m_Results[i];
  xmlOut << "</testsuite>\n";
  xmlOut.close();

  return success;
}



bool RegressionTest::syncDirectory(const char *targetName, const char *localName, 
                                   const char *hostUrl, std::string &errorMsg) {
  std::string targetPath(targetName);
  std::string localPath(localName);
  targetPath = Fs::convertToUncPath(targetPath);
  localPath = Fs::convertToUncPath(localPath);
#if defined (WIN32)
  errorMsg = "Please copy data in" + targetPath + " to " + localPath;
  return false;
#else
  /* On unix, if possible try to just link the data as this is quick. */
  if(Fs::isReadableDir(targetPath)) {
    int retVal = 0;
    std::string command = ToStr("ln -s ") + targetPath + "  " + localPath;
    retVal = system(command.c_str());
    if(retVal == 0) {
      return true;
    }
    else {
      errorMsg = "Failed trying to run command: " + command;
      return false;
    }
  }
  /* Can't see that directory, try to download it. */
  else {
    std::string command = ToStr("wget -m -P ") + localPath  + 
      " " + ToStr(hostUrl);
    int retVal = 0;
    Verbose::out(1, "Retrieving remote data using wget.", false);
    retVal = system(command.c_str());
    Verbose::out(1, "\tDone.");
    if(retVal == 0) {
      return true;
    }
    else {
      errorMsg = "Failed trying to run command: " + command;
    }
  }
#endif
  return false;
}

void RegressionTest::parseLogFile(const std::string &logFile) {
  RowFile rf;
  rf.open(logFile);

  int maxMemVirtual = 0;
  vector<string> words;
  // Skip header...
  rf.nextRow(words);
  while (rf.nextRow(words)) {
    int m = Convert::toInt(words[10]);
    maxMemVirtual = Max(m, maxMemVirtual);
  }
  rf.close();
  m_MaxMem = maxMemVirtual;
}

void RegressionTest::parseValgrindFile(const std::string &vgLog) {
  ifstream in;
  string errToken =  "== ERROR SUMMARY: ";
  string leakToken = "== LEAK SUMMARY:";
  m_MemOk = -1;
  m_MemLeaks = -1;
  Fs::aptOpen(in,vgLog);
  APT_ERR_ASSERT(in.is_open() && in.good(), "Couldn't open file: " + vgLog + " to read.");
  string line;
  while (!in.eof()) {
    getline(in, line);
    if (line.find(errToken) != string::npos) {
      // Parse memory errors.
      int startPos = errToken.size() + line.find(errToken);
      int endPos = line.find(" ", startPos);
      APT_ERR_ASSERT(endPos != string::npos, "Couldn't parse error summary.");
      string countStr = line.substr(startPos, endPos-startPos);
      int count = Convert::toInt(countStr);
      if (count == 0) {
        m_MemOk = 1;
      }
      else {
        m_MemOk = 0;
      }
    }
    if (line.find(leakToken) != string::npos) {
      // Parse memory leaks...
      // Get the "definitely lost line."
      getline(in, line);
      string defLostToken = "definitely lost: ";
      int foundToken = line.find(defLostToken);
      int startPos = defLostToken.length() + foundToken;;
      int endPos = line.find(" ", startPos);
      APT_ERR_ASSERT(endPos != string::npos, "Couldn't parse memory leak.");
      string countStr = line.substr(startPos, endPos-startPos);
      Util::removeChar(countStr,',');
      int count = Convert::toInt(countStr);
      m_MemLeaks = count;
      // Parse leaks
    }
  }
  in.close();
}

void RegressionTest::recordResults(bool success) {
  APT_ERR_ASSERT(!m_Database.empty(), "Database not set.");
  if (m_SQL == NULL) {
    m_SQL = new SQLiteDatabase();
    m_SQL->open(m_Database);
  }
  
  stringstream sql;
  sql << "insert into regression_log values ("
      << "'" << m_SuiteName << "', "
      << "'" << m_Name << "', " 
      <<  m_TimeStamp << ", " 
      << (success ? 1 : 0) << ", " 
      << (m_NegTest ? 1 : 0) << ", " 
      << m_64bit << ", "
      << "'" << m_Build << "', "
      << "'" << m_OpSys << "', "
      << (m_ValgrindLog.empty() ? 0 : 1) << ", "
      << m_MaxMem << ", "
      << m_RunTime << ", "
      << m_MemOk << ", "
      << m_MemLeaks << ", "
      << m_NumPassed << ", "
      << m_NumFailed 
      << ");";
  string toExecute = sql.str();
  Verbose::out(1, toExecute);
  try {
    m_SQL->execute(toExecute, true, false);
  }
  catch ( SQLiteException e) {
    Verbose::warn(1, std::string("RegressionTest::recordResults FAILED: ") + e.getMessage()); 
  }
  catch ( Except &e) {
    Verbose::warn(1, std::string("RegressionTest::recordResults FAILED exception: ") + e.what()); 
  }
  catch (...) {
    Verbose::warn(1, std::string("RegressionTest::recordResults FAILED with unknown exception."));
  }
}

void RegressionTest::setSuite(const RegressionSuite &suite, 
                              const std::string &outDir,
                              const std::string &logDir, 
                              const std::string &valgrindLog) {

  if ( !Fs::dirExists(outDir) ) {
    Fs::mkdirPath(outDir, false);
  }
  setLogFile(logDir);
  setSuiteName(suite.getSuiteName());
  setDatabase(suite.getDatabase());
  if (suite.doValgrind()) {
    setValgrindLog(valgrindLog);
  }
}
