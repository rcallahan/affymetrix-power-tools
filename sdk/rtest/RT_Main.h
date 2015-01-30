////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
 * @file   RT_Main.h
 * @author Alexander Wong
 * @date   Sun Oct  4 11:49:37 2009
 * @brief  Class used to contain and execute a suite of regression tests
 */

#ifndef RT_MAIN_H
#define RT_MAIN_H

#include "rtest/RT_Test.h"


/// @class RT_Main
/// @brief Definition of options and the values.
class RT_Main
{
 public:

  ///@brief Creates and initializes a new RT_Main
  RT_Main(const std::string &name);
  ~RT_Main();

  /// @brief Retrieves the library of environment variables
  /// @return Returns a pointer to m_mainEnv
  RT_Env* getMainEnv();

  /// @brief Sets or creates a new variable
  /// @param varKey=variable key or tag  
  /// @param varVal=value associated with that key
  void setVar(const std::string &varKey, const std::string &varVal);  

  /// @brief Sets or creates a new variable, parses the input string
  /// @param varString=string containing variable key and value in the form "varKey=varValue"  
  void setVar(const std::string &varString);  

  /// @brief Retrieves the variable value associated with the input key
  /// @param varKey=Variable Key
  /// @return Returns the Variable Value in the form of a String
  std::string getVarVal(const std::string &varKey);

  /// @brief Stores a list of default test names associated with this suite of regression tests
  /// @param testNameList a string vector containing test names
  void setTestNameList(std::vector<std::string> testNameList);

  /// @brief Retrieves a list of default test names associated with this suite of regression tests
  /// @param Returns a string vector containing test names
  std::vector<std::string> getTestNameList();
  
  /// @brief Add or create a new test and store it in m_tests
  /// @param test = pointer to an instance of an RT_Test
  void addTest(RT_Test* test);

  /// @brief Runs all tests stored in m_tests
  void runTests();

  /// @brief Defines acceptable Vars for Regression Testing
  /// @param opts=pointer to an instance of PgOptions
  void defineOptions(PgOptions* opts);
  /// @brief Parses the command line arguments passed in and stores them in m_vars
  /// @param argc=count of number of elements in the argv array
  /// @param argv=array of characters terminated by NULL
  void parseArgv(int argc, char* argv[]);
  /// @brief Parse a simple bash file into our RT Framework
  void parseFile();

  /// @brief Prints out values stored in all class variables
  void dump();
  /// @brief Generates a quick reference file of test results
  /// @param passedTests = a List of all regression tests passed during this run 
  /// @param failedTests = a List of all regression tests failed during this run
  void printQuickRefFile(std::vector<std::string> passedTests, std::vector<std::string> failedTests);
  
  void printTests();

 private:
  std::string m_name;                              ///< Name of this suite of tests
  std::vector<RT_Test*> m_tests;                   ///< Contains all Tests to be run
  RT_Env* m_mainEnv;                               ///< Object holding environment variables
  int m_numFailed;                                 ///< Number of tests failed
  int m_numPassed;                                 ///< Number of tests passed
  std::vector<std::string> m_defaultTestNameList;  ///< Contains a list of all the possible default tests
};

#endif
