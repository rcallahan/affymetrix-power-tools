////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file   SummaryGenotypeEngine.h
 * @author Luis Jevons
 * @date   Wed April 17 16:01:09 2009
 * 
 * @brief Core routines for summary-genotype binaries. By separating the
 * command line parsing form the computation we allow a GUI application to share
 * the core computation once setting up the option class. 
 */
#ifndef _SUMMARYGENOTYPEENGINE_H_
#define _SUMMARYGENOTYPEENGINE_H_

#include "file/TsvFile/TsvFile.h"
#include "util/BaseEngine.h"
#include "file5/File5_Tsv.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * Class to wrap up the complexity of the different possible input
 * files to SummaryGenotypeEngine
 */
class SgIoHelper {
  
public: 
  
  /** 
   * Constructor. Take a number of files and figure out which one(s) has
   * been set to read from.
   * 
   * @param summariesFile   Text summary file with results interleved.
   * @param summariesFileA  Text summary file for A alleles
   * @param summariesFileB  Text summary file for B alleles
   * @param file5           File5 version of the summary file above
   * @param file5Group      Which group the summary file5 is in with name "summaries"
   */
  SgIoHelper(const std::string &summariesFile,
           const std::string &summariesFileA, const std::string &summariesFileB,
           const std::string &file5,
           const std::string &file5Group);  


  /** Destructor */
  ~SgIoHelper();

  /** 
   * Read in a single line from out summaries file and fill label and values.
   * 
   * @param label      Probesetname to be filled in
   * @param vals       Values to be filled in
   * @param summaries  File to read next row from
   * @param colNames   Names of the column in summaries file
   */
  void getLine(std::string &label, std::vector<double> &vals, 
               affx::TsvFile &summaries, std::vector<std::string> &colNames);

  /** 
   * Read in a single line from out summaries file5 and fill label and values.
   * 
   * @param label      Probesetname to be filled in
   * @param vals       Values to be filled in
   * @param summaries  File5 to read next row from
   * @param colNames   Names of the column in summaries file
   */
  void getFile5Line(std::string &label, std::vector<double> &vals, 
                    affx::File5_Tsv *summaries, std::vector<std::string> &colNames);

  /** 
   * Read the next snp marker from the single interleved file5 summary file
   *  
   * @param colNames Names of the columns in summary file.
   * @param aVals Summary values for the A allele
   * @param bVals Summary values for the B allele
   * @param probesetName Name of probeset
   * 
   * @return true if successful, false otherwise
   */
  bool getNextMarkerTsvSingle(std::vector<std::string> &colNames, std::vector<double> &aVals, 
                              std::vector<double> &bVals, std::string &probesetName);

  /** 
   * Read the next snp marker from the file5 
   *  
   * @param colNames Names of the columns in summary file.
   * @param aVals Summary values for the A allele
   * @param bVals Summary values for the B allele
   * @param probesetName Name of probeset
   * 
   * @return true if successful, false otherwise
   */
  bool getNextMarkerFile5(std::vector<std::string> &colNames, std::vector<double> &aVals, 
                          std::vector<double> &bVals, std::string &probesetName);
  /** 
   * Read the next snp marker from the paired summary A & B files
   *  
   * @param colNames Names of the columns in summary file.
   * @param aVals Summary values for the A allele
   * @param bVals Summary values for the B allele
   * @param probesetName Name of probeset
   * 
   * @return true if successful, false otherwise
   */
  bool getNextMarkerTwoTsv(std::vector<std::string> &colNames, std::vector<double> &aVals, 
                           std::vector<double> &bVals, std::string &probesetName);
  /** 
   * Read the next snp marker.
   *  
   * @param colNames Names of the columns in summary file.
   * @param aVals Summary values for the A allele
   * @param bVals Summary values for the B allele
   * @param probesetName Name of probeset
   * 
   * @return true if successful, false otherwise
   */
  bool getNextMarker(std::vector<std::string> &colNames, std::vector<double> &aVals, 
                     std::vector<double> &bVals, std::string &probesetName);

  /** 
   * Get the names of the columns in the summary files.
   * 
   * @return names in order of the columns
   */
  std::vector<std::string> getColNames() {
    return m_ColNames;
  }

private:
  
  /// Column names from the summary files not including the probeset_id column
  std::vector<std::string> m_ColNames;
  /// All the summaries as a single interleved TSV file
  affx::TsvFile summaries;
  /// All the A allele summaries as a TSV File
  affx::TsvFile summariesA;
  /// All the B allele summaries as a TSV File
  affx::TsvFile summariesB;
  /// Summaries in a file5 format
  affx::File5_File f5Summaries;
  /// Within the file5 format, the group that the intensities are in.
  affx::File5_Group *f5SumGroup;
  /// The tsv file within the File5 that contains the "summaries"
  affx::File5_Tsv *f5SumTsv;
};

class SummaryGenotypeEngine : public BaseEngine {

public:

  virtual std::string getEngineName() { return SummaryGenotypeEngine::EngineName(); }
  static const std::string EngineName() { return "SummaryGenotypeEngine"; }

  /**
   * Constructor
   */
  SummaryGenotypeEngine();

  /**
   * Destructor
   */
  ~SummaryGenotypeEngine();

  /*! A class to register the summary engine. */
  class Reg : public EngineReg
  {
  public:
    /*! Constructor - register the summary engine. */
    Reg() : EngineReg(SummaryGenotypeEngine::EngineName())
    {
    }

    /*! Creates an object.
     * @return The object.
     */
    BaseEngine *MakeObject() { return new SummaryGenotypeEngine; }
  };

  /*! The one and only registration object. */
  static Reg reg;

  /*! Converts the type to the summary engine type.
   * @param chip The pointer to the base engine object.
   * @return The summary engine type or NULL if not compatible.
   */
  static SummaryGenotypeEngine * FromBase(BaseEngine *engine);

private:

  virtual void defineOptions();
  virtual void checkOptionsImp();
  virtual void defineStates();
  virtual void runImp();
};

#endif
