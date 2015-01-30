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
#ifndef RT_TEST_H
#define RT_TEST_H

#include "rtest/RT_Cmd.h"
#include "rtest/RT_Check.h"
#include "rtest/RT_Util.h"


class RT_Test
{
 public:
  RT_Test();
  RT_Test(const std::string &name);
  RT_Test(const std::string &name, RT_Env* env);
  ~RT_Test(); 

  ///@brief Creates a new RT_Cmd to store in the m_cmd list of commands
  RT_Cmd* newCmd();

  ///@brief Creates a new RT_Check to store in the m_check list of Checks
  RT_Check* newCheck();

  /// @brief Sets the environment to point to the env argument
  /// @param env = pointer to a instance of RT_Env
  void setEnv(RT_Env* env);

  /// @brief sets a new environment variable
  /// @param varKey=variable key that acts as an identifier
  /// @param varValue=value associated with the varKey
  void setVar(const std::string &varKey, const std::string &varVal);

  /// @brief Overloaded method to add or create a RT_Cmd object to the test 
  void addCmd(RT_Cmd* cmd);
  void addCmd(const std::string &exe,std::vector<std::string> args);
  
  /// @brief Overloaded method to add or create a RT_Check object to the test 
  void addCheck(RT_Check* check);
  void addCheck(const std::string &exe, std::vector<std::string> args);
  void addCheck(const std::string &exe, std::vector<std::string> args, bool negTest);

  /// @brief Generates multiple new checks and stores them in this test 
  /// @param Takes in the Check wrapper executable, and two vector of file names (gen and gold) which must be the same size
  void addCheckMult(const std::string &exe, std::vector<std::string> gen, std::vector<std::string> gold, std::vector<std::string> args=std::vector<std::string>());

  /// @brief Runs all commands and Checks within this test
  /// @return Returns 1 for failed test, 0 for passed test
  int run();

  /// @brief Adds a subtest to this test
  void addSubtest(RT_Test* test);
  
  /// @brief Retrives class variable m_testName
  std::string getName();

  /// @brief Set the test as positive or negative test
  void setExpectedRetVal(int retVal);

  /// @brief Retrives m_expectedRetVal which signifies a positive or negative test
  int getExpectedRetVal();

  /// @brief Prints out values stored in all class variables
  void dump();

 private:
  std::vector<RT_Cmd*> m_cmd;              ///< Contains all apt executables to be run for this test
  std::vector<RT_Check*> m_check;          ///< Contains all check wrapper executables to be run for this test
  std::string m_testName;                  ///< Name of this test
  int m_expectedRetVal;                    ///< Flag for a negative test; 0 = positive, 1 = negative test
  std::vector<RT_Test*> m_subtests;        ///< All subtests which comprise this test.  Subtests are also RT_Test's
  RT_Env* m_testEnv;                       ///< Environment which holds our dictionary of variables for this test
};

#endif
