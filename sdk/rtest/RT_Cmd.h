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
#ifndef RT_CMD_H
#define RT_CMD_H

#include "rtest/RT_Args.h"


class RT_Cmd
{
 public:
  
  RT_Cmd();
  RT_Cmd(const std::string &exe, std::vector<RT_Args*> args);
  ~RT_Cmd();

  /// @brief sets the m_exe class variable
  void setExe(const std::string &exe);
  /// @brief Retrieves m_exe class variable
  std::string getExe();

  /// @brief adds an argument to the m_args class vector variable
  void addArg(const std::string &argName, const std::string &argVal);
  void addArg(const std::string &argString);
  /// @brief Retrieves m_args class variable
  std::vector<RT_Args*> getArgs();

  /// @brief Returns a string representing a call to a apt executable (i.e. apt-probeset-genotype, apt-cel-extract). It is built up by combining m_exe and m_args
  std::string getCmd();

  /// @brief Prints out values stored in all class variables
  void dump();

 protected:
  std::vector<RT_Args*> m_args;        ///< Contains all Arguments to be passed to the apt executable on the command line
  std::string m_exe;                   ///< Contains the name and path to the apt executable
};

#endif
