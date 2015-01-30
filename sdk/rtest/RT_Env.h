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
#ifndef RT_ENV_H
#define RT_ENV_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include "util/Util.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"

//Define new type from std::map to hold environment variables
typedef std::map<std::string, std::string> maptype;

class RT_Env
{
 public:
  RT_Env();
  ~RT_Env();
  
  /// @brief Retrieves the m_vars class variable
  maptype getVars();

  /// @brief Define all default variables and assign them a default value
  void defineDefaultVars();

  /// @brief Returns Variable value with name that matches the queryString
  std::string getVarVal(const std::string &varKey);

  /// @brief Create new variable or edit existing variable to have passed in Values
  /// @param varKey=variable key or tag  
  /// @param varVal=value associated with that key
  void setVar(const std::string &varKey, const std::string &varVal);
  /// @brief Create new variable or edit existing variable to have passed in Values
  /// @param varString=string containing variable key and value in the form "varKey=varValue"  
  void setVar(const std::string &varString);

  /// @brief Takes in a input string and replaces all instances of variable key values with the variable values 
  std::string expandVarString(const std::string &strArg);

  /// @brief Retrieves the cpucomsys operating system environment variable
  std::string getCpuComSys();

  /// @brief Prints out values stored in all class variables
  void dump();

 private:
  maptype m_vars;
};


#endif
