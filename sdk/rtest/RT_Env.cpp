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
#include "rtest/RT_Env.h"
#include <stdlib.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Env::RT_Env()
{
  this->defineDefaultVars();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Env::~RT_Env()
{
    //Delete Stuff Here
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
maptype RT_Env::getVars()
{
  return this->m_vars;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns Variable value with name that matches the varName
//Returns empty string if no match found
std::string RT_Env::getVarVal(const std::string &varKey)
{
  std::string retStr = "";
  maptype::const_iterator iter = this->m_vars.find(varKey);
  if(iter != this->m_vars.end())
    {
      //Return the value in the map
      retStr = iter->second;
      return retStr;
    }
  else
    {
      retStr = "";
      
      return retStr;
    }
  return retStr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create new variable or edit existing variable to have passed in Values
void RT_Env::setVar(const std::string &varKey, const std::string &varVal)
{
  this->m_vars[varKey] = varVal;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create new variable or edit existing variable to have passed in Values
void RT_Env::setVar(const std::string &varString)
{
  size_t pos = 0;
  
  pos = varString.find("=");
      
  if(pos == string::npos)
    {
      Verbose::out(0, "Error: Could not set variable");
      exit(1);
    }

  std::string varKey = varString.substr(0, pos);
  std::string varVal = varString.substr(pos+1, varString.length()-pos);
  
  this->m_vars[varKey] = varVal;
  

  /*
  int varLength = varString.length();
  const char* varChar = varString.c_str();
  bool isKey = true;

  std::string varKey = "";
  std::string varVal = "";

  for(int i = 0;i < varLength;i++)
    {
      if(varChar[i] == '=')
	{
	  isKey = false;
	}
      else
	{
	  if(isKey == true)
	    {
	      varKey = varKey + varChar[i];
	    }
	  else
	    {
	      varVal = varVal + varChar[i];
	    }
	}
    }
  this->m_vars[varKey] = varVal;
  */
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string RT_Env::getCpuComSys()
{
  ///@TODO Handle Windows case
  char* str = getenv("CPUCOMSYS");
  if(str==NULL)
    {
      return std::string("unknown");
    }

  std::string returnString = str;
  return returnString;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Define all variables and assign them a default value
void RT_Env::defineDefaultVars()
{
  this->setVar("sdktop", "../../..");
  this->setVar("alex", getCpuComSys());
  this->setVar("testname", "");

  this->setVar("valgrind", "valgrind --leak-check=yes --log-file=${regout_dir}/valgrind/${testname}");
  this->setVar("doValgrind", "false");

  this->setVar("regout_dir", "${sdktop}/rtest/testOutput");
  this->setVar("exe_dir", "${sdktop}/output/${alex}/bin");
  this->setVar("out_dir", "test-generated");
  this->setVar("gold_dir", "${sdktop}/regression-data/data");
  this->setVar("out_testname", "${out_dir}/${testname}");

  this->setVar("apt_probeset_genotype", "${exe_dir}/apt-probeset-genotype");
  this->setVar("apt_probeset_summarize", "${exe_dir}/apt-probeset-summarize");
  this->setVar("apt_cel_extract", "${exe_dir}/apt-cel-extract");
  this->setVar("apt_cel_transformer", "${exe_dir}/apt-cel-transformer");
  this->setVar("apt_geno_qc", "${exe_dir}/apt-geno-qc");
  this->setVar("apt_summary_genotype", "${exe_dir}/apt-summary-genotype");
  this->setVar("apt_dmet_genotype", "${exe_dir}/apt-dmet-genotype");
  this->setVar("apt_canary", "${exe_dir}/apt-canary");
  this->setVar("apt_cel_convert", "${exe_dir}/apt-cel-convert");
  this->setVar("apt_chp_to_txt", "${exe_dir}/apt-chp-to-txt");
  this->setVar("apt_matrix_diff", "${exe_dir}/apt-matrix-diff");
  this->setVar("apt_tsv_join", "${exe_dir}/apt-tsv-join");
  this->setVar("apt_cdf_export", "${exe_dir}/apt-cdf-export");

  this->setVar("apt_check_matrix", "${exe_dir}/apt-check-matrix");
  this->setVar("apt_check_mixedfile", "${exe_dir}/apt-check-mixedfile");
  this->setVar("apt_check_cel", "${exe_dir}/apt-check-cel");
  this->setVar("apt_check_chp", "${exe_dir}/apt-check-chp");
  this->setVar("apt_check_calvinchp", "${exe_dir}/apt-check-calvinchp");
  this->setVar("apt_check_textfile", "${exe_dir}/apt-check-textfile");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
std::string RT_Env::expandVarString(const std::string &strArg)
{
  size_t openPos = 0;
  std::string str = strArg;
  while(1)
    {
      //Find the first variable found starting from where we left off
      openPos=str.find("${", openPos);

      //If we have parsed the whole string, then just return it
      if(openPos==string::npos)
	{
	  return str;
	}
      else
	{
	  //Start where the current opening bracket was found
	  size_t closePos;
	  closePos = str.find("}", openPos);
	  
	  std::string subString = str.substr(openPos+2, closePos-2-openPos);
	  std::string varVal = getVarVal(subString);
	  
	  varVal = expandVarString(varVal);
    
	  str.erase(openPos, closePos+1-openPos);
	  str.insert(openPos, varVal);

	  openPos = openPos + varVal.length();
	}
    }
}
  /*
  int strLength = strArg.length();
  const char* charArray = strArg.c_str();
  bool isVar = false;
  bool foundDollar = false;
  bool isVarExpanded = false;

  std::string returnStr = "";
  std::string varKey = "";

  //Checking two characters at a time for varval delimiter
  //If found, replaces variable keys with their values
   for(int i = 0;i < strLength;i++)
    {
      //Is this the start of a variable?
      //Found a dollar sign...possible variable initiation
      if(charArray[i] == '$')
	{
	  foundDollar = true;
	}
      //Ok, we found a variable
      else if(charArray[i] == '{')
	{
	  if(foundDollar == true)
	    {
	      isVar = true;
	      foundDollar = false;
	    }
	}
      //Variable name done.
      //Find value and add to output string
      else if(charArray[i] == '}' && isVar == true)
	{
	  returnStr = returnStr + getVarVal(varKey);
	  //Reset the varKey to empty
	  varKey = "";
	  isVar = false;
	  isVarExpanded = true;
	}
      else
	{
	  //Did not find a bracket immediately following the dollar sign
	  //Was a false variable declaration
	  if(foundDollar = true)
	    {
	      foundDollar = false;
	    }
	  
	  //We found a variable
	  //Start building up the variable name/key
	  if(isVar == true)
	    {
	      varKey = varKey + charArray[i];
	    }
	  //If we aren't in a variable just save contents from original string
	  else
	    {
	      returnStr = returnStr + charArray[i];
	    }
	}
    }
   
   //If a variable was expanded, there is a possibility it was replaced with a value containing another variable
   //Therefore we need to recursively call this method again
   if(isVarExpanded == true)
     {
       returnStr = expandVarString(returnStr);
     }
   
   return returnStr;
  
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
void RT_Env::dump()
{
  Verbose::out(2, "***RT_VARS DUMP***");
  
  for(maptype::iterator i=this->m_vars.begin();i != this->m_vars.end();i++)
    {
      std::string str1 = "Key: "+ i->first + "\n";
      std::string str2 = "Value: " + i->second + "\n";
      Verbose::out(2, str1);
      Verbose::out(2, str2);
    }
}
