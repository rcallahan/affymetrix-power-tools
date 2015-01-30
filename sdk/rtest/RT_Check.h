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
#ifndef RT_CHECK_H
#define RT_CHECK_H

#include "rtest/RT_Args.h"
#include "rtest/RT_Cmd.h"
#include <iostream>


class RT_Check : public RT_Cmd
{
 public:
  //Standard Constructor
  RT_Check();
  RT_Check(const std::string &exe, std::vector<RT_Args*> args);
  //Constructor for negative test/check
  RT_Check(const std::string &exe, std::vector<RT_Args*> args, int negTest);
  //RT_Check(const std::string &exe, std::vector<RT_Args*> args, std::vector<std::string> gen, std::vector<<std::string> gold, RT_Test* rt);
  ~RT_Check();

  /// @brief sets the m_exe class variable
  //void setExe(const std::string &exe);
  /// @brief Retrieves m_exe class variable
  //std::string getExe();

  /// @brief adds an argument to the m_args class vector variable
  //void addArg(const std::string &argName, const std::string &argVal);
  //void addArg(const std::string &argString);
  /// @brief Retrieves m_args class variable
  //std::vector<RT_Args*> getArgs();

  void setNegTest();

  void setGenFiles(std::vector<std::string> genFiles);
  std::vector<std::string> getGenFiles();

  void setGoldFiles(std::vector<std::string> goldFiles);
  std::vector<std::string> getGoldFiles();

  bool isGenGoldEmpty();
  void setGenGoldFileIsArg(bool check);
  bool getGenGoldFileIsArg();
  std::vector<std::string> getCheckMult();

  /// @brief Returns a string representing a call to a check wrapper executable (i.e. apt-check-matrix, apt-check-calvinchp). It is built up by combining m_exe and m_args
  std::string getCheck();

  /// @brief Prints out values stored in all class variables
  void dump();


 private:
  //std::vector<RT_Args*> m_args;          ///< Contains all Arguments to be passed to the Check Wrapper
  //std::string m_exe;                     ///< Contains the name and path to the Check Wrapper
  std::vector<std::string> m_genFiles;   ///< Contains a list of generated Files to be used in the creation of multiple checks from this single RT_Check
  std::vector<std::string> m_goldFiles;  ///< Contains a list of gold Files to be used in the creation of multiple checks from this single RT_Check
  int m_negTest;                         ///< Flag for a negative test; 0 = positive, 1 = negative test
  bool m_genGoldFileIsArg;             ///< Flag to see if gen and gold vectors should be added as arguments instead of with flags. Depends on the check wrapper being used.
};

#endif
