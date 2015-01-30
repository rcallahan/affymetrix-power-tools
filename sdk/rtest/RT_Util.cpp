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
#include "RT_Util.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Takes in a vector of string arguments and returns a vector of RT_Args
std::vector<RT_Args*> RT_Util::makeArgVec(std::vector<std::string> stringVec)
{
  //Take string vectore and make Name-Value argument pairs
  std::vector<RT_Args*> args;
  std::vector<RT_Args*>::iterator argIt = args.begin();
  bool nameFlag = false;
  RT_Args* tempArg;
	      
  //Assumption: Each option flag or "NAME" has only a single "Value" corresponding to it
  for(std::vector<std::string>::iterator i = stringVec.begin();i < stringVec.end();i++)
    {
      std::string str = *i;
      //If it is a option flag, then put it in "Name"
      if(str[0] == '-' && str[1] == '-')
	{
	  //If a flag was found previously but no value for it was given
	  if(nameFlag == true)
	    {
	      //Since this is just a "Name" or flag
	      //Set the Value equal to the "Name"
	      tempArg->setValue(tempArg->getName());
	      args.push_back(tempArg);
	      argIt++;
	      nameFlag = false;
	    }
	  //If we are at the end, we need to make sure we store the last value in the vector
	  if(i == stringVec.end()-1)
	    {
	      tempArg = new RT_Args;
	      std::string tempstr = "";
	      tempstr = str;
	      tempArg->setName(tempstr);
	      tempArg->setValue(tempstr);
	      args.push_back(tempArg);
	      argIt++;
	      nameFlag = false;
	    }

	  //Add new element as Name
	  nameFlag = true;
	  tempArg = new RT_Args;
	  std::string name;
	  name = str;
	  tempArg->setName(name);
	}
      //Only single dash
      else if(str[0] == '-')
	{
	  //If a flag was found previously but no value for it was given
	  if(nameFlag == true)
	    {
	      //Since this is just a "Name" or flag
	      //Set the Value equal to the "Name"
	      tempArg->setValue(tempArg->getName());
	      args.push_back(tempArg);
	      argIt++;
	      nameFlag = false;
	    }
	  //If we are at the end, we need to make sure we store the last value in the vector
	  if(i == stringVec.end()-1)
	    {
	      tempArg = new RT_Args;
	      std::string tempstr = "";
	      tempstr = str;
	      tempArg->setName(tempstr);
	      tempArg->setValue(tempstr);
	      args.push_back(tempArg);
	      argIt++;
	      nameFlag = false;
	    }
	  
	  //Add new element as Name
	  nameFlag = true;
	  tempArg = new RT_Args;
	  std::string name;
	  name = str;
	  tempArg->setName(name);
	}
      //If an option flag was found previously
      else if(nameFlag == true)
	{
	  nameFlag = false;
	  tempArg->setValue(str);
	  args.push_back(tempArg);
	  argIt++;
	  }
      //We have an argument without a flag...likely and input file name
      else if(nameFlag == false)
	{
	  tempArg = new RT_Args;
	  tempArg->setName(str);
	  tempArg->setValue(str);		      
	  args.push_back(tempArg);
	  argIt++;
	}
    }
  return args;	      
}
