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
#include "rtest/RT_Cmd.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Cmd::RT_Cmd()
{
  this->m_exe = "";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Cmd::RT_Cmd(const std::string &exe, std::vector<RT_Args*> args)
{
  this->m_exe = exe;
  this->m_args.assign(args.begin(), args.end());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Cmd::~RT_Cmd()
{
  //Delete Stuff Here
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//sets the  m_exe class variable
void RT_Cmd::setExe(const std::string &exe)
{
  this->m_exe = exe;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns/retrieves m_exe class variable
std::string RT_Cmd::getExe()
{
  return this->m_exe;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///Adds a arg element to m_args vector
void RT_Cmd::addArg(const std::string &argName, const std::string &argVal)
{
  RT_Args* arg = new RT_Args(argName, argVal);
  this->m_args.push_back(arg);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///Adds a arg element to m_args vector
void RT_Cmd::addArg(const std::string &argString)
{
  size_t pos = 0;
  
  pos = argString.find("=");
      
  //If no equal sign found, then this is a flag or variable without a value
  if(pos == string::npos)
    {
      RT_Args* arg = new RT_Args(argString, argString);
      this->m_args.push_back(arg);
    }
  else
    {
      //Otherwise, if a equal sign is found, then grab the name and value separately
      std::string argName = argString.substr(0, pos);
      std::string argVal = argString.substr(pos+1, argString.length()-pos);
      
      RT_Args* arg = new RT_Args(argName, argVal);
      this->m_args.push_back(arg);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns/retrieves m_args class variable
std::vector<RT_Args*> RT_Cmd::getArgs()
{
  return this->m_args;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a string representing a call to a apt executable (i.e. apt-probeset-genotype, apt-cel-extract)
//Is built up from m_exe and m_args
std::string RT_Cmd::getCmd()
{
  std::string cmd = this->m_exe;
 
  for(std::vector<RT_Args*>::iterator i = this->m_args.begin();i < m_args.end();i++)
    {
      RT_Args* arg = *i;
      //If we have a name-value pair where name is equal to value, only grab it once
      if(arg->getName() == arg->getValue())
	{
	  cmd += " " + arg->getValue();
	}
      //Otherwise we want both the name and its value separately passed in
      else
	{
	  cmd += " " + arg->getName() + " "  + arg->getValue();
	}
    }
  return cmd;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
void RT_Cmd::dump()
{
  Verbose::out(2, "***RT_Cmd DUMP***");
  Verbose::out(2, "Exe: "+ this->m_exe);
  for(std::vector<RT_Args*>::iterator i = this->m_args.begin();i < this->m_args.end();i++)
    {
      RT_Args* arg;
      arg = *i;
      arg->dump();
    }
}
