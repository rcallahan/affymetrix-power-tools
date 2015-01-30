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
#include "rtest/RT_Test.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test::RT_Test()
{
  this->m_testName = "";
  this->m_expectedRetVal = 0;
  this->m_testEnv = new RT_Env;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test::RT_Test(const std::string &name)
{
  this->m_testName = name;
  this->m_expectedRetVal = 0;
  this->m_testEnv = new RT_Env;
  this->m_testEnv->setVar("test_name", name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test::RT_Test(const std::string &name, RT_Env* env)
{
  this->m_testName = name;
  this->m_expectedRetVal = 0;
  
  //Create new instance of env, a copy of the contents of the env argument pointer
  RT_Env tempEnv = *env;
  this->m_testEnv = &tempEnv;
  this->m_testEnv->setVar("test_name", name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test::~RT_Test()
{
  //Delete Stuff Here
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set the Env Variable
//Takes in a pointer to another RT_Env instance
//Creates a new instance from this and stores it for this test
void RT_Test::setEnv(RT_Env* env)
{
  //Create new instance of RT_Env and make it a copy of the contents of the env pointer argument
  RT_Env* tempEnv = new RT_Env;
  *tempEnv = *env;
  this->m_testEnv = tempEnv;

  //Set the testname variable to the name of current test
  std::string name = this->getName();
  this->m_testEnv->setVar("testname", name);

  //For any existing subtests, copy the test's new env and change the name to the name of the subtest
  for(std::vector<RT_Test*>::iterator i=this->m_subtests.begin();i<this->m_subtests.end();i++)
    {
      //grab the iterator pointer contents
      RT_Test* subtest = *i;

      //Create new instance of RT_Env and make it a copy of the contents of the main test's
      RT_Env* subEnv = new RT_Env;
      *subEnv = *tempEnv;
      subtest->setEnv(subEnv);

      //Set the testname variable for this subtest
      std::string subName = subtest->getName();
      subtest->setVar("testname", subName);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create or modify a env Variable
void RT_Test::setVar(const std::string &varKey, const std::string &varVal)
{
  this->m_testEnv->setVar(varKey, varVal);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns m_testName
std::string RT_Test::getName()
{
  return this->m_testName;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Cmd* RT_Test::newCmd()
{
  RT_Cmd* cmd = new RT_Cmd();
  this->m_cmd.push_back(cmd);

  return cmd;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Check* RT_Test::newCheck()
{
  RT_Check* check = new RT_Check();
  this->m_check.push_back(check);

  return check;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//to add or create a RT_Cmd object to the test 
void RT_Test::addCmd(RT_Cmd* cmd)
{
  this->m_cmd.push_back(cmd);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//to add or create a RT_Cmd object to the test 
void RT_Test::addCmd(const std::string &exe, std::vector<std::string> args)
{

  std::vector<RT_Args*> cmdArgs;
  cmdArgs = RT_Util::makeArgVec(args);

  RT_Cmd* cmd;
  cmd = new RT_Cmd(exe, cmdArgs);

  this->m_cmd.push_back(cmd);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//to add or create a RT_Check object to the test 
void RT_Test::addCheck(RT_Check* check)
{
  this->m_check.push_back(check);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//to add or create a RT_Check object to the test 
void RT_Test::addCheck(const std::string &exe, std::vector<std::string> args)
{
  std::vector<RT_Args*> checkArgs;
  checkArgs = RT_Util::makeArgVec(args);

  RT_Check* check;
  check = new RT_Check(exe, checkArgs);

  this->m_check.push_back(check);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//to add or create a RT_Check object to the test 
void RT_Test::addCheck(const std::string &exe, std::vector<std::string> args, bool negTest)
{
  std::vector<RT_Args*> checkArgs;
  checkArgs = RT_Util::makeArgVec(args);

  RT_Check* check;
  check = new RT_Check(exe, checkArgs, negTest);

  this->m_check.push_back(check);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generates multiple Check objects from the gen and gold vectors (Breaks them down into file pair entries for new Checks) 
//Any Arguments other than gen and gold should be pushed on to an rt_args vec and passed in.
void RT_Test::addCheckMult(const std::string &exe, std::vector<std::string> gen, std::vector<std::string> gold, std::vector<std::string> args)
{
  //Must be equal since we are creating checks with file pairs from gen and gold.
  if(gen.size() == gold.size())
    {
      for(int i = 0;i < gen.size(); i++)
	{
	  //Store gen and gold file
	  std::vector<std::string> argVec;
	  argVec.push_back("--gen");
	  argVec.push_back(gen[i]);
	  argVec.push_back("--gold");
	  argVec.push_back(gold[i]);
  
	  //Add any additional arguments
	  for(std::vector<std::string>::iterator j = args.begin();j < args.end();j++)
	    {
	      argVec.push_back(*j);
	    }

	  this->addCheck(exe, argVec);
	}
    }
  else
    {
      Verbose::out(0, "Error in RT_Test variables for generateMultiChecks.  Incorrect number of arguments given in gen or gold vector");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Adds a subtest to this test
void RT_Test::addSubtest(RT_Test* test)
{
  test->setEnv(this->m_testEnv);
  this->m_subtests.push_back(test);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Runs all Cmds and Checks within this RT_Test
int RT_Test::run()
{
  int passed = 0;

  Verbose::out(1, "#TEST\n" + this->getName());
  
  //Iterate through and run all subtests
  for(std::vector<RT_Test*>::iterator m = this->m_subtests.begin();m < m_subtests.end();m++)
    {
      RT_Test* test = *m;
      int subRetVal = 0;
      subRetVal = test->run();

      //Check Return Value for pass or fail
      if(subRetVal != 0)
	{
	  passed = 1;
	}
    }
  
  //Iterate through and run all Commands
  for(std::vector<RT_Cmd*>::iterator i = this->m_cmd.begin();i < m_cmd.end();i++)
    {
      //Get the command string and expand any variables within it to there values
      RT_Cmd* cmd = *i;
      std::string cmdString = cmd->getCmd();
      if(this->m_testEnv->getVarVal("doValgrind") == "true")
	{
	  cmdString = this->m_testEnv->getVarVal("valgrind") + " " + cmdString;
	}
      cmdString = this->m_testEnv->expandVarString(cmdString);

      Verbose::out(1, "#CMD\n" + cmdString);
      ///@TODO This should be a fork/exec
      int retVal = system(cmdString.c_str());

      //If Positive Test
      if(this->m_expectedRetVal == 0)
	{
	  //Something went wrong in executing the command line
	  if(retVal == -1)
	    {
	      Verbose::out(0, "Error: Command did not pass or run in " + this->m_testName);
	      //std::string errStr = "Error: Command did not pass or run in " + this->m_testName;
	      //Verbose::out(0, errStr);
	      passed = 1;
	    }
	  //The test Failed
	  else if(retVal != 0)
	    {
	      passed = 1;
	    }
	}
      //If Negative Test
      else
	{
	  if(retVal == -1)
	    {
	      Verbose::out(0, "Error: Command did not pass or run in " + this->m_testName);
	      //std::string errStr = "Error: Command did not pass or run in " + this->m_testName;
	      //Verbose::out(0, errStr);
	      passed = 1;
	    }
	  //The Test Failed
	  if(retVal == 0)
	    {
	      passed = 1;
	    }
	}
    }
  
  //Iterate through and run all Checks
  for(std::vector<RT_Check*>::iterator j = this->m_check.begin();j < m_check.end();j++)
    {
      //Get the check string and expand any variables within it to there values
      RT_Check* check = *j;

      //If the gold and gen vectors are not empty, then we want to create multiple checks
      if(!check->isGenGoldEmpty())
	{
	  std::vector<std::string> multiCheckString = check->getCheckMult();
	  for(std::vector<std::string>::iterator m = multiCheckString.begin();m < multiCheckString.end();m++)
	    {
	      std::string checkString = *m;
	      checkString + "--testName" + this->m_testName;
	      checkString = this->m_testEnv->expandVarString(checkString);

	      Verbose::out(1, "#CHECK\n" + checkString);
	      int retVal = system(checkString.c_str());

	      //Check Return Value
	      //Something went wrong in executing the command line
	      if(retVal == -1)
		{
		  Verbose::out(0, "Error: In " + this->m_testName + " a Check did not pass or run");
		  passed = 1;
		}
	      //Test Failed
	      else if(retVal != 0)
		{
		  passed = 1;
		}
	    }
	}
      //Otherwise just get a single check
      else
	{
	  std::string checkString = check->getCheck();
	  checkString = checkString + " --testName " + this->m_testName;
	  checkString = this->m_testEnv->expandVarString(checkString);

	  Verbose::out(1, "#CHECK\n" + checkString);
	  int retVal = system(checkString.c_str());

	  //Check Return Value
	  if(retVal == -1)
	    {
	      Verbose::out(0, "Error: In " + this->m_testName + " a Check did not pass or run");
	      passed = 1;
	    }
	  else if(retVal != 0)
	    {
	      passed = 1;
	    }
	}
    }
  return passed;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
void RT_Test::dump()
{
  Verbose::out(2, "***RT_Test Dump***");
  Verbose::out(2, "Test Name: " + this->m_testName);
  for(std::vector<RT_Cmd*>::iterator i = this->m_cmd.begin();i < this->m_cmd.end();i++)
    {
      RT_Cmd* cmd;
      cmd = *i;
      
      cmd->dump();
    }
  
  for(std::vector<RT_Check*>::iterator j = this->m_check.begin();j < this->m_check.end();j++)
    {
      RT_Check* check;
      check = *j;
      
      check->dump();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set the test as positive or negative test
//Changes m_expectedRetVal class variable
void RT_Test::setExpectedRetVal(int retVal)
{
  this->m_expectedRetVal = retVal;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns m_expectedRetVal
int RT_Test::getExpectedRetVal()
{
  return this->m_expectedRetVal;
}

