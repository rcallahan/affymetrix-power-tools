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

/**
 * @file   RegressionTest.h
 * @author Chuck Sugnet
 * @date   Sun Dec  4 23:26:39 2005
 * 
 * @brief  Utility class for helping to do regression testing.
 */
#ifndef REGRESSION_TEST_H
#define REGRESSION_TEST_H

#include "util/RegressionSuite.h"
#include "util/RegressionCheck.h"
#include "util/MatrixCheck.h"
#include "util/MixedFileCheck.h"
//
#include <fstream>
#include <iostream>
#include <string>
//

class SQLiteDatabase;

/**
 * Utility class for helping to do regression testing.
 */
class RegressionTest {

public:

  /** Constructor. */
  RegressionTest(const std::string& name,
                 const std::string& generatedFile, 
                 const std::string& goldFile,
                 double epsilon, 
                 const std::string& command,
                 int rowSkip,
                 int colSkip,
                 bool matchNames,
                 int allowedMisMatch,
                 bool negTest = false);
  
  RegressionTest(const std::string& name,
                 const std::string& command, 
                 const std::vector<RegressionCheck *>& fileChecks,
                 bool doDelete=true,
                 bool negTest = false);

  ~RegressionTest();

  /** 
   * Set the log from LogStream to parse memory and other stats from.
   * 
   * @param logFile - path to log file
   */
  void setLogFile(const std::string &logFile) { m_LogToParse = logFile; }

  /** 
   * Set the valgrind --log-file to read the memory leak and checks from.
   * 
   * @param vgLog - path to valgrind log file
   */
  void setValgrindLog(const std::string &vgLog) { m_ValgrindLog = vgLog; }

  void setDatabase(const std::string &path);

  void setSuiteName(const std::string &suiteName) { m_SuiteName = suiteName; }

  int getNumPassed() const { return m_NumPassed; }

  int getNumFailed() const { return m_NumFailed; }
  
  /** 
   * Run the command to check the test. 
   * @return true if executes sucessfully, false otherwise.
   */
  bool run();

  /** 
   * Do we pass our post-run tests??? 
   * @return true if we pass our test
   */
  bool pass();

  /** 
   * Get whatever message was reported. 
   *
   * @return the concatanation of any error messages that have been
   * sent.
   */
  std::string getErrorMsg() { return m_Error; }

  /** 
   * Try to link and if necessary copy a directory from another place.
   * @param targetName - path to directory to be synced
   * @param localName - where we want data synced to locally.
   * @param hostUrl - url for data if copying necessary.
   * @param errorMsg - fill in error message here if something goes wrong.
   * 
   * @return - true if sucessful, false otherwise.
   */  
  static bool syncDirectory(const char *targetName, const char *localName, 
                            const char *hostUrl, std::string &errorMsg);


  void setSuite(const RegressionSuite &suite, 
                const std::string &outDir,
                const std::string &logDir, 
                const std::string &valgrindLog);

public:
  std::string m_Name;
  bool m_NegTest;

protected:

  /** Constructor */
  RegressionTest();

  /** 
   * Set appropriate data members to no data values
   */
  void initialize();

  /** 
   * Read a log file to populate things like memory used, etc.
   * 
   * @param logFile - File produced by LogStream
   */
  void parseLogFile(const std::string &logFile);

  /** 
   * Read through valgrind output to determine memory leaks and errors. Note
   * that according to the valgrind manual it is best to compile the program
   * with -O0 (optimization level zero) to avoid false positives.
   * http://valgrind.org/docs/manual/quick-start.html#quick-start.prepare
   * 
   * @param vgLog - valgrind output with the --log-file option
   */
  void parseValgrindFile(const std::string &vgLog);

  void recordResults(bool success);

  static const int DATA_NA = -1;

  std::vector<RegressionCheck *> m_Checks; 
  std::string m_CommandStr;     ///< Command that generates the generated file, platform adjusted.
  std::string m_Error;          ///< Any errors generated are stored here.
  std::string m_SuiteName;
  bool m_Delete;
  int64_t m_RunTime; ///< How long did the test take to run in seconds
  int m_MaxMem; ///< Maximum resident virtual memory during processing
  int m_MemOk;  ///< 0 false, 1 true, -1 not set
  int m_MemLeaks; ///< Bytes definitely leaked according to valgrind
  int m_TimeStamp;
  bool m_DoMemCheck;
  std::string m_ValgrindLog;
  std::string m_LogToParse;
  std::string m_Database;
  std::string m_OpSys;
  std::string m_Build;
  bool m_64bit;
  int m_NumPassed;
  int m_NumFailed;
  SQLiteDatabase *m_SQL;
};

#endif /* REGRESSIONTEST_H */
