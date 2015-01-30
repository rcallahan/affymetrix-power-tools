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

/**
 * @file   apt-snp-model-converter-test
 * @author Mybrid Spalding
 * @date   Tue Jan 25 15:25:43 PST 2011
 *
 * @brief  Program for doing regression tests on models files.
 *
 *
 */

#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/RegressionTest.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <string>
#include <vector>
//

class ModelFileValidCheck : public RegressionCheck
{

private:
  std::string  m_modelFile;
  std::string  m_outFile;
  bool  m_invalidCheck;
  
public:
  ModelFileValidCheck(const std::string &modelFile, const std::string outFile, bool testInvalid = false) :
    m_modelFile(modelFile), m_outFile(outFile), m_invalidCheck(testInvalid) {
    m_Name = modelFile;
  }

  virtual ~ModelFileValidCheck() {}

  bool check(std::string & msg) {
    //AxiomGT1.snp-posteriors.txt valid model file detected.
    Verbose::out(1, m_modelFile + " check...", false);
    std::ifstream results;
    std::string validTest = "valid model file detected.";
    std::string invalidTest = "not a recognized model file.";
    Fs::aptOpen(results, m_outFile);
    std::string line;
    char buf[BUFSIZ];
    results.getline(buf, BUFSIZ);
    line = ToStr(buf);
    bool valid = true;
    while(valid && !line.empty()) {
      if(line.find(m_modelFile) != std::string::npos) {
        if(m_invalidCheck) {
          if(line.find(invalidTest) == std::string::npos) {
            valid = false;
          }
        } else if(line.find(validTest) == std::string::npos) {
          valid = false;
        }
      }
      results.getline(buf, BUFSIZ);
      line = ToStr(buf);
    }
    results.close();

    if(valid) {
      Verbose::out(1, "ok");
    } else {
      Verbose::out(1, "failed");
    }
    return valid;
  }

};

class ModelFileHeaderCheck : public RegressionCheck
{

  std::string  m_outDir;
  std::string  m_logFile;
  std::string  m_libFile;
  std::string  m_goldFile;
  std::string  m_generatedFile;
  
public:

  ModelFileHeaderCheck(const std::string outDir, const std::string logFile,  const std::string &libFile, const std::string gold, std::string gen) :
    m_outDir(outDir), m_logFile( logFile ), m_libFile(libFile), m_goldFile(gold), m_generatedFile(gen) {
    m_Name = Fs::basename(libFile);
  }

  virtual ~ModelFileHeaderCheck() {}

  bool check(std::string & msg) {

    std::string ignoreFile = m_outDir + "/ignore.txt";
    Fs::rm(m_logFile, false);
    Fs::convertToUncPathInPlace(m_generatedFile,10);

    std::string cmd =
      "./apt-snp-model-converter --input-file " +  m_libFile +
      " --log-file " + m_logFile + 
      " --dump-headers GTC4.1 > \"" + m_generatedFile + "\" 2>" + ignoreFile;

    int retVal = system(Fs::convertCommandToUnc(cmd).c_str());
    if(retVal != 0) {
      Verbose::warn(1, ToStr(" ") + cmd + " Failed to run command.", "\nERROR: ");
      return false;
    }

    TextFileCheck textCheck(m_generatedFile, m_goldFile, 0);
    std::string checkMsg = m_libFile + " validated GTC4.1 headers.";
    bool checkRv = textCheck.check(checkMsg);
    if ( checkRv ) {
      Verbose::out(1, checkMsg);
    }
    else {
      Verbose::warn(1, checkMsg, "\nERROR: ");
    }
    return checkRv;
  }

};

class SnpModelConverterTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  SnpModelConverterTest() {
    numPassed = 0;
    numFailed = 0;
  }

  void qt_doModelInvalid(const std::string & outDir);
  void qt_doModelValid(const std::string & outDir);
  void qt_doModelHeaders(const std::string & outDir);
};

void SnpModelConverterTest::qt_doModelValid(const std::string & outDir)
{

  std::string name = "qt-doModelValid";
  std::string testDir = Fs::join(outDir, name);
  Fs::mkdir(testDir, false);
  std::string logFile = Fs::join(testDir, name + ".log");
  Fs::rm(logFile, false);
  std::string ignoreFile = Fs::join(testDir, "ignore.txt");
  
  std::vector< std::string > modelFiles;
  const std::string modelFilesListing = "../../../regression-data/data/idata/snp-model-converter/qt-doModelValid_test.txt";
  APT_ERR_ASSERT(affx::TsvFile::extractColToVec(modelFilesListing, "model-files",
                 &modelFiles) == affx::TSV_OK,
                 modelFilesListing + std::string(" invalid file column header 'model-files' is required."));
  std::string command = std::string("./apt-snp-model-converter --model-files ") +
                        modelFilesListing  + " --out-dir " + testDir + " --log-file " + logFile +
                        " > " + ignoreFile + " 2>&1 ";

  std::vector<RegressionCheck *> checks;

  for(int i = 0; i < modelFiles.size(); i++) {
    if(Fs::fileExists(modelFiles[i])) {
      checks.push_back(new ModelFileValidCheck(modelFiles[i], logFile));
    }
  }

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpModelConverterTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }
}

void SnpModelConverterTest::qt_doModelInvalid(const std::string & outDir)
{

  std::string name = "qt-doModelInvalid";
  std::string testDir = Fs::join(outDir, name);
  Fs::mkdir(testDir, false);
  std::string logFile = Fs::join(testDir, name + ".log");
  Fs::rm(logFile, false);
  std::string ignoreFile = Fs::join(testDir, "ignore.txt");

  std::vector< std::string > invalidFiles;
  const std::string invalidFilesListing = "../../../regression-data/data/idata/snp-model-converter/qt-doModelInvalid_test.txt";
  APT_ERR_ASSERT(affx::TsvFile::extractColToVec(invalidFilesListing, "model-files",
                 &invalidFiles) == affx::TSV_OK,
                 invalidFilesListing + std::string(" invalid file column header 'model-files' is required."));
  std::string command = std::string("./apt-snp-model-converter --model-files ") +
                        invalidFilesListing  + " --out-dir " + testDir + " --log-file " + logFile +
                        " > " + ignoreFile + " 2>&1 ";

  std::vector<RegressionCheck *> checks;

  for(int i = 0; i < invalidFiles.size(); i++) {
    if(Fs::fileExists(invalidFiles[i])) {
      checks.push_back(new ModelFileValidCheck(invalidFiles[i], logFile, true));
    }
  }

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpModelConverterTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }
}

void SnpModelConverterTest::qt_doModelHeaders(const std::string & outDir)
{

  std::string name = "qt-doModelHeaders";
  std::string testDir = Fs::join(outDir, name);
  Fs::mkdir(testDir, false);
  std::string logFile = Fs::join(outDir, name + ".log");


  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;
  std::string testFiles = "../../../regression-data/data/idata/snp-model-converter/qt-doModelHeaders_test.tsv";
  APT_ERR_ASSERT(tsv.open(testFiles) == affx::TSV_OK, testFiles + " unable to open");

  std::string libPath, goldPath, generatedPath;
  tsv.bind(0, "lib_path", &libPath, affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "results_path", &goldPath, affx::TSV_BIND_REQUIRED);

  std::vector<RegressionCheck *> checks;
  while(tsv.nextLevel(0) == affx::TSV_OK) {
    checks.push_back(new ModelFileHeaderCheck(testDir, logFile, libPath, goldPath,
                     testDir + "/" + Fs::basename(goldPath)));
  }

  tsv.clear();
  
  std::string command = "echo Running apt-snp-model-converter per lib file...";
  RegressionTest test(name.c_str(), command.c_str(), checks);

  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpModelConverterTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }
}

/** Everybody's favorite function. */
int main(int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("chipstream/snp-model-converter-qt", true);
    
    SnpModelConverterTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.qt_doModelInvalid(test.testDir);
    test.qt_doModelValid(test.testDir);
    test.qt_doModelHeaders(test.testDir);
    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  } catch(...) {
    Verbose::out(1, "Unexpected Error: uncaught exception.");
    return 1;
  }
  return 1;
}

