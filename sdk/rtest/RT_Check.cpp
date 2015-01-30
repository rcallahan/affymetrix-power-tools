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
#include "rtest/RT_Check.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Standard Constructor
RT_Check::RT_Check()
{
  this->m_exe = "";
  this->m_genGoldFileIsArg = false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Standard Constructor
RT_Check::RT_Check(const std::string &exe, std::vector<RT_Args*> args)
{
  this->m_exe = exe;
  this->m_args.assign(args.begin(), args.end());
  this->m_genGoldFileIsArg = false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for negative test/check
RT_Check::RT_Check(const std::string &exe, std::vector<RT_Args*> args, int negTest)
{
  this->m_exe = exe;
  this->m_args.assign(args.begin(), args.end());
  this->m_genGoldFileIsArg = false;
  
  //Create new arg to flag for negative test
  //This will be passed to and picked up by the check wrapper 
  if(negTest == 1)
    {
      RT_Args* arg;
      arg = new RT_Args("--expectedRetVal", "1");
      this->m_args.push_back(arg);
    }
  else
    {
        RT_Args* arg;
	arg = new RT_Args("--expectedRetVal", "0");
	this->m_args.push_back(arg);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Check::~RT_Check()
{
  //Delete Stuff Here
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Check::setNegTest()
{
  //Create new arg to flag for negative test
  //This will be passed to and picked up by the check wrapper 
  RT_Args* arg;
  arg = new RT_Args("--expectedRetVal", "1");
  this->m_args.push_back(arg);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
void RT_Check::dump()
{
  Verbose::out(2, "***RT_CHECK DUMP***");
  Verbose::out(2, "EXE: " + this->m_exe);
  for(std::vector<RT_Args*>::iterator i = this->m_args.begin();i < this->m_args.end();i++)
    {
      RT_Args* myArg;
      myArg = *i;
      
      myArg->dump();
    }
}
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////
//sets the  m_exe class variable
void RT_Check::setExe(const std::string &exe)
{
  this->m_exe = exe;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns/retrieves m_exe class variable
std::string RT_Check::getExe()
{
  return this->m_exe;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///Adds a arg element to m_args vector
void RT_Check::addArg(const std::string &argName, const std::string &argVal)
{
  RT_Args* arg = new RT_Args(argName, argVal);
  this->m_args.push_back(arg);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///Adds a arg element to m_args vector
void RT_Check::addArg(const std::string &argString)
{
  size_t pos = 0;
  
  pos = argString.find("=");
      
  //If no equal sign found, then this is a flag or variable without a value
  if(pos == string::npos)
    {
      RT_Args* arg = new RT_Args(argString, argString);
      this->m_args.push_back(arg);
    }

  //Otherwise, if a equal sign is found, then grab the name and value separately
  std::string argName = argString.substr(0, pos);
  std::string argVal = argString.substr(pos+1, argString.length()-pos);
  
  RT_Args* arg = new RT_Args(argName, argVal);
  this->m_args.push_back(arg);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns/retrieves m_args class variable
std::vector<RT_Args*> RT_Check::getArgs()
{
  return this->m_args;
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Check::setGenFiles(std::vector<std::string> genFiles)
{
  this->m_genFiles = genFiles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> RT_Check::getGenFiles()
{
  return this->m_genFiles;
} 

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Check::setGoldFiles(std::vector<std::string> goldFiles)
{
  this->m_goldFiles = goldFiles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> RT_Check::getGoldFiles()
{
  return this->m_goldFiles;
} 

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool RT_Check::isGenGoldEmpty()
{
   if(!this->m_genFiles.empty() && !this->m_goldFiles.empty())
    {
      int genSize = this->m_genFiles.size();
      int goldSize = this->m_goldFiles.size();

      if(genSize != goldSize)
	{
	  Verbose::out(0, "Error in the size of the gen and gold vectors for this check.  Exiting.");
	  exit(1);
	}
      return false;
    }
   else
     return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Check::setGenGoldFileIsArg(bool check)
{
  this->m_genGoldFileIsArg = check;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool RT_Check::getGenGoldFileIsArg()
{
  return this->m_genGoldFileIsArg;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> RT_Check::getCheckMult()
{
  if(this->isGenGoldEmpty())
    {
      Verbose::out(0, "Gen and Gold Vectors are empty! Exiting.");
      exit(1);
    }
  else
    {
      //If we want the gold and gen files to be added as an arg at the end of hte line
      if(this->getGenGoldFileIsArg() == true)
	{
	  std::vector<std::string> checkMult;
	  int genSize = this->m_genFiles.size();

	  for(int i = 0; i < genSize;i++)
	    {
	      std::string check = this->m_exe;

	      //Add any additional args
	      for(std::vector<RT_Args*>::iterator j = this->m_args.begin();j < m_args.end();j++)
		{
		  RT_Args* arg = *j;

		  //If we have a name-value pair where name is equal to value, only grab it once
		  if(arg->getName() == arg->getValue())
		    {
		      check += " " + arg->getValue();
		    }
		  //Otherwise we want both the name and its value separately passed in
		  else
		    {
		      check += " " + arg->getName() + " "  + arg->getValue();
		    }
		}
	      std::string gen = this->m_genFiles[i];
	      std::string gold = this->m_goldFiles[i];
	      //Add gold and gen files
	      check += " " + gen;
	      check += " " + gold; 

	      checkMult.push_back(check);
	    }
	  return checkMult;
	}
      else
	{
	  std::vector<std::string> checkMult;
	  int genSize = this->m_genFiles.size();

	  for(int i = 0; i < genSize;i++)
	    {
      
	      std::string check = this->m_exe;
	      std::string gen = "--gen " + this->m_genFiles[i];
	      std::string gold = "--gold " + this->m_goldFiles[i];
	      //Add gold and gen files
	      check += " " + gen;
	      check += " " + gold; 

	      //Add any additional args
	      for(std::vector<RT_Args*>::iterator j = this->m_args.begin();j < m_args.end();j++)
		{
		  RT_Args* arg = *j;

		  //If we have a name-value pair where name is equal to value, only grab it once
		  if(arg->getName() == arg->getValue())
		    {
		      check += " " + arg->getValue();
		    }
		  //Otherwise we want both the name and its value separately passed in
		  else
		    {
		      check += " " + arg->getName() + " "  + arg->getValue();
		    }
		}
	      checkMult.push_back(check);
	    }
	  return checkMult;
	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a string representing a call to a check executable
//Is built up from m_exe and m_args
std::string RT_Check::getCheck()
{
  std::string check = this->m_exe;

  for(std::vector<RT_Args*>::iterator i = this->m_args.begin();i < m_args.end();i++)
    {
      RT_Args* arg = *i;

      //If we have a name-value pair where name is equal to value, only grab it once
      if(arg->getName() == arg->getValue())
	{
	  check += " " + arg->getValue();
	}
      //Otherwise we want both the name and its value separately passed in
      else
	{
	  check += " " + arg->getName() + " "  + arg->getValue();
	}
    }

  return check;
}
