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
#ifndef RT_ARGS
#define RT_ARGS

#include "rtest/RT_Env.h"

class RT_Args
{
 public:
  RT_Args();
  RT_Args(const std::string &name, const std::string &value);
  ~RT_Args();

  //Assign value to m_Name or m_Value
  void setName(const std::string &name);
  void setValue(const std::string &value);

  //Retrieve m_Name or m_Value
  std::string getName();
  std::string getValue();

  //Print out values stored in the RT_Args variables
  void dump();

 private:
  std::string m_Name;      //denotes the name of the command, the flag
  std::string m_Value;     //denotes the value or argument to be used in conjunction with the name/flag
};
#endif
