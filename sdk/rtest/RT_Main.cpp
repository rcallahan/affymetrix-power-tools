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
#include "rtest/RT_Main.h"
#include "util/Fs.h"
//
#include <cstdlib>
#include <sstream>
//
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Main::RT_Main(const std::string &name)
{
  this->m_name = name;
  this->m_numPassed = 0;
  this->m_numFailed = 0;
  this->m_mainEnv = new RT_Env();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
RT_Main::~RT_Main()
{
  //Delete Stuff Here
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out values stored in all class variables
void RT_Main::dump()
{
  Verbose::out(2, "***RT_Main DUMP***");
  for(std::vector<RT_Test*>::iterator i = this->m_tests.begin();i < this->m_tests.end();i++)
    {
      RT_Test* test;
      test = *i;
      test->dump();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Main::printTests()
{
  Verbose::out(1, "***List of Default Tests***");
  for(std::vector<std::string>::iterator i=this->m_defaultTestNameList.begin();i < this->m_defaultTestNameList.end();i++)
    {
      std::string testName = *i;
      Verbose::out(1, testName);
    }
  Verbose::out(1, "\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Store test name list for help/tutorial reasons
void RT_Main::setTestNameList(std::vector<std::string> testNameList)
{
  this->m_defaultTestNameList = testNameList;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> RT_Main::getTestNameList()
{
  return this->m_defaultTestNameList;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add/create a new test and store it in m_tests
void RT_Main::addTest(RT_Test* test)
{
  test->setEnv(this->m_mainEnv);
  this->m_tests.push_back(test);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Runs all tests stored in m_tests
void RT_Main::runTests()
{
  std::vector<std::string> passedTests;
  std::vector<std::string> failedTests;

  for(std::vector<RT_Test*>::iterator i = m_tests.begin();i < m_tests.end();i++)
    {
      RT_Test* test = *i;
      int passed = 0;
      passed = test->run();
      //If test Passed
      if(passed == 0)
	{
	  this->m_numPassed++;
	  passedTests.push_back(test->getName());
	  Verbose::out(1, "Passed " + test->getName() + "\n\n");
	}
      //If test failed
      else
	{
	  this->m_numFailed++;
	  failedTests.push_back(test->getName());
	  Verbose::out(1, "Failed " + test->getName() + "\n\n");
	}
    }

  Verbose::out(1, "\n\nTest Results: ");
  std::string sPassed;
  std::stringstream passout;
  passout<<this->m_numPassed;
  sPassed = passout.str();

  std::string sFailed;
  std::stringstream failout;
  failout<<this->m_numFailed;
  sFailed = failout.str();
  
  std::string numPassedString = "NumPassed: " + sPassed;
  std::string numFailedString = "NumFailed: " + sFailed;
  Verbose::out(1, numPassedString);
  Verbose::out(1, numFailedString);
  //Verbose::out(1,  "NumPassed: " + this->m_numPassed + ", NumFailed: " + this->m_numFailed + "\n");
  this->printQuickRefFile(passedTests, failedTests);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generates a quick reference file of test results
void RT_Main::printQuickRefFile(std::vector<std::string> passedTests, std::vector<std::string> failedTests)
{
  ofstream fileout;
  std::string filename = "";
  //Check to set outfile to the proper directory depending on the test
  if(this->m_name == "apt-rt-dmet-genotype" || this->m_name == "apt-rt-canary")
    {
      filename = "../../rtest/testOutput/" + this->m_name + "-QResults";
    }
  else
    {
      filename = "../../../rtest/testOutput/" + this->m_name + "-QResults";
    }
  fileout.open(filename.c_str(), ios::out); 
  fileout<<"****Test Results "<<this->m_name<<" *****"<<endl;
  fileout<<"====================================="<<endl;
  fileout<<"NUMPASSED: "<<this->m_numPassed<<endl;
  fileout<<"NUMFAILED: "<<this->m_numFailed<<endl;
  fileout<<endl;
  
  if(this->m_numPassed != 0)
    {
      fileout<<"Passed Tests: "<<endl;
      for(std::vector<std::string>::iterator j = passedTests.begin();j < passedTests.end();j++)
	{
	  fileout<<*j<<endl;
	}
    }

  cout<<endl;

  if(this->m_numFailed != 0)
    {
      fileout<<"Failed Tests: "<<endl;
      for(std::vector<std::string>::iterator m = failedTests.begin();m < failedTests.end();m++)
	{
	  fileout<<*m<<endl;
	}
    }
  
  fileout.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//REturns a pointer to m_mainEnv
RT_Env* RT_Main::getMainEnv()
{
  return this->m_mainEnv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Main::setVar(const std::string &varKey, const std::string &varVal)
{
  this->m_mainEnv->setVar(varKey, varVal);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void RT_Main::setVar(const std::string &varString)
{
  this->m_mainEnv->setVar(varString);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns true if the variable exists in m_vars, false otherwise
std::string RT_Main::getVarVal(const std::string &varKey)
{
  std::string varVal =  this->m_mainEnv->getVarVal(varKey);
  return varVal;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Defines acceptable Vars for Regression Testing
void RT_Main::defineOptions(PgOptions* opts)
{
  opts->defineOption("s", "simple", PgOpt::STRING_OPT,
		     "specifies a simple bash command script to create a regression test.  This is an alternative to the usual command line regressiontest calls.",
		     "");
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
		     "Produces a help page/tutorial for the Regression Testing Framework",
		     "false");
  opts->defineOption("l", "listtests", PgOpt::BOOL_OPT,
		     "Lists the tests",
		     "false");
  opts->defOptMult("D", "varval", PgOpt::STRING_OPT,
		   "set var=val for this test.  This overrides the default value for this run or creates a new variable completely",
		   "");
  opts->defineOption("v", "valgrind", PgOpt::BOOL_OPT,
		     "Run Valgrind on all commands, will print to the specified log file",
		     "false");
  opts->defineOption("e", "verbose", PgOpt::INT_OPT,
		     "Set the Verbosity Level",
		     "2");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parses the command line arguments passed in and stores them in m_vars
void RT_Main::parseArgv(int argc, char* argv[])
{
  PgOptions* opts = new PgOptions();
      
  //Set Usage message for Matrix Check
  opts->setUsage("The command line format for running a regression test:\napt-rt-program --simple simplebashtest.txt --D var=val test1 test2\n\n");

  this->defineOptions(opts);

  ///Parse argv and store it within PgOptions class variables
  ///Match the command line arguments to the proper options
  opts->parseArgv(argv);

  //If the User wanted Help
  if(opts->getBool("help"))
    {
      int choice = -1;
      std::string waitingForEnter = "";

      while(choice != 0)
	{
	  cout<<endl;
	  cout<<endl;
	  cout<<"****Regression Test Help Menu****"<<endl;
	  cout<<"1. How-to run regression tests"<<endl;
	  cout<<"2. Default tests"<<endl;
	  cout<<"3. Command line Options"<<endl;
	  cout<<"4. Default Variables"<<endl;
	  cout<<"(Enter 0 to Exit)"<<endl;
	  cout<<endl;
	  cout<<"Please enter your number choice from the list above:";
	  cin>>choice;
	  cout<<endl;

	  //How to run examples
	  if(choice == 1)
	    {
	      cout<<"The command line format for running a regression test:\napt-rt-program --simple simplebashtest.txt --D var=val test1 test2\n\n"<<endl;
	      cout<<endl;
	      cout<<"Press any letter key then Enter to continue"<<endl;
	      cin>>waitingForEnter;
	    }
	  //List of Built in tests
	  else if(choice == 2)
	    {
	      int testNum = 0;
	      cout<<"List of Default Tests"<<endl;
	      for(std::vector<std::string>::iterator i=this->m_defaultTestNameList.begin();i < this->m_defaultTestNameList.end();i++)
		{
		  std::string testName = *i;
		  testNum++;
		  cout<<testNum<<". "<<testName<<endl;
		}
	      cout<<endl;
	      cout<<"Press any letter key then Enter to continue"<<endl;
	      cin>>waitingForEnter;
	    }
	  //List of command line options
	  else if(choice == 3)
	    {
	      opts->usage();
	      cout<<endl;
	      cout<<"Press any letter key then Enter to continue"<<endl;
	      cin>>waitingForEnter;
	    }
	  //List of default Var=Val environment variables 
	  else if(choice == 4)
	    {
	      cout<<"List of Default Variables"<<endl;
	      maptype vars = this->m_mainEnv->getVars();
	      for ( maptype::iterator it=vars.begin(); it != vars.end(); it++)
		{
		  cout << (*it).first << " => " << (*it).second << endl;
		}

	      cout<<endl;
	      cout<<"Press any letter key then Enter to continue"<<endl;
	      cin>>waitingForEnter;
	    }
	  else if(choice == 0)
	    {
	      continue;
	    }
	  else
	    {
	      cout<<"Invalid Choice"<<endl;
	      choice = 0;
	    }
	}
      exit(0);	  
    }

  //If one or more variables were defined
  if(opts->get("varval") != "")
    {
      PgOpt* oneOpt;
      oneOpt = opts->mustFindOpt("varval");

      std::vector<std::string> myVec;
      
      myVec = oneOpt->getValueVector();

      for(std::vector<std::string>::iterator j=myVec.begin();j<myVec.end();j++)
	{
	  std::string str = *j;
	  this->setVar(str);
	}
    }

  if(opts->getBool("valgrind"))
    {
      this->setVar("doValgrind", "true");
      if( !Fs::dirExists("${regout_dir}/valgrind") ) {
        Fs::mkdirPath(this->m_mainEnv->expandVarString("${regout_dir}/valgrind"), false);
      }
    }

  if(opts->getBool("listtests"))
    {
      this->printTests();
      exit(0);
    }

  //If a simple Bash File was specified
  if(opts->get("simple") != "")
    {
      std::string value = opts->get("simple");
      this->setVar("simple", value);
      this->parseFile();
    }

  if(opts->getInt("verbose") != 2)
    {
      int level = opts->getInt("verbose");
      Verbose::setLevel(level);
    }
  else
    {
      Verbose::setLevel(2);
    }
  
  //Get name of default tests to run
  std::vector<std::string> args = opts->getArgVector();
  for(std::vector<std::string>::iterator i = args.begin();i < args.end();i++)
    {
      std::string value = *i;
      this->setVar(value, value);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parse a simple bash file into our RT Framework
void RT_Main::parseFile()
{
  RT_Test* rt = NULL;

  std::string fileName;
  fileName =  this->getVarVal("simple");
  std::ifstream inFile;
  inFile.open(fileName.c_str());
  
  std::string line;
  std::string fileString;

  bool testFlag = false;
  bool cmdFlag = false;
  bool checkFlag = false;
  bool varFlag = false;

  while(getline(inFile, line))
    {
      //Skip blank lines
      if(line == "")
	{
	  continue;
	}
      else if(line == "#test")
	{
	  testFlag = true;
	  cmdFlag = false;
	  checkFlag = false;
	  varFlag = false;
	}
      else if(line == "#cmd")
	{
	  cmdFlag = true;
	  testFlag = false;
	  checkFlag = false;
	  varFlag = false;
	}
      else if(line == "#check")
	{
	  checkFlag = true;
	  testFlag = false;
	  varFlag = false;
	  cmdFlag = false;
	}
      else if(line == "#var")
	{
	  varFlag = true;
	  testFlag = false;
	  checkFlag = false;
	  cmdFlag = false;
	}
      else
	{
	  //If Test flag found, we should find the test name and create a new RT_Test instance
	  if(testFlag == true)
	    {
	      //Store previously built up test
	      if(rt != NULL)
		{
		  this->addTest(rt);
		}
	      rt = new RT_Test(line);
	    }
	  //If Var Flag found, we should save and store each successive variable declared
	  if(varFlag == true)
	    {
	      //If no Test has yet been declared, create a new test instance with default name
	      if(rt == NULL)
		{
		  rt = new RT_Test("The Test With No Name");
		}
	      std::vector<std::string> lineArgs;
	      char delimiter = '=';
	      std::string exe;
	      //Remove Dollar Sign preceding variable name
	      int removeDollarSign = 1;

	      //Break the line into a series of strings based upon white space as delimiter
	      Util::chopString(line, delimiter, lineArgs);
	      
	      std::string varKey = "";
	      std::string varVal = "";
	      
	      //Process variable name/key
	      std::string tempString = lineArgs[0];
	      int strLength = tempString.length();
	      const char* charArray = tempString.c_str();
	      for(int i = removeDollarSign;i < strLength;i++)
		{
		  varKey = varKey + charArray[i];
		}
	      
	      //Get variable value
	      varVal = lineArgs[1];

	      this->m_mainEnv->setVar(varKey,varVal);	 
	    }
	  //If this is a Command...
	  else if(cmdFlag == true)
	    {
	      //If no Test has yet been declared, create a new test instance with default name
	      if(rt == NULL)
		{
		  rt = new RT_Test("The Test With No Name");
		}
	      std::vector<std::string> lineArgs;
	      char delimiter = ' ';
	      std::string exe;

	      //Parse the given string and expand any known variables within it
	      //line = this->m_mainEnv->expandVarString(line);
	      //Break the line into a series of strings based upon white space as delimiter
	      Util::chopString(line, delimiter, lineArgs);

	      //Iterate through each word from the line
	      for (std::vector<std::string>::iterator i = lineArgs.begin(); i < lineArgs.end();i++)
		{
		  //If it's the exe, store it separately
		  if(i == lineArgs.begin())
		    {
		      exe = *i;
		      lineArgs.erase(i);
		    }
		}
	      
	      //Add the Command to the test
	      rt->addCmd(exe, lineArgs);
	    }
	  //If this is a Check...
	  else if(checkFlag == true)
	    {
	      //If no Test has yet been declared, create a new test instance with default name
	      if(rt == NULL)
		{
		  rt = new RT_Test("The Test With No Name");
		}
	      
	      std::vector<std::string> lineArgs;
	      char delimiter = ' ';
	      std::string exe;

	      //Parse the given string and expand any known variables within it
	      //line = this->m_mainEnv->expandVarString(line);

	      //Break the line into a series of strings based  upon white space as delimiter
	      Util::chopString(line, delimiter, lineArgs);

	      //Iterate through each word from the line
	      for (std::vector<std::string>::iterator i = lineArgs.begin(); i < lineArgs.end();i++)
		{
		  //If it's the exe pull it out
		  if(i == lineArgs.begin())
		    {
		      exe = *i;
		      lineArgs.erase(i);
		    }
		}

	      //Add the Check to this test
	      rt->addCheck(exe, lineArgs);
	    }
	}
    }
  //Add this built up test to this RT_Main instance
  this->addTest(rt);
}
