////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
 * @file   EngineUtil.h
 * @author Chuck Sugnet
 * @date   Tue Jul 11 11:49:05 2006
 * 
 * @brief  Some common functions that application executables make use of.
 */

#ifndef _ENGINEUTIL_H_
#define _ENGINEUTIL_H_

//
#include "chipstream/AnalysisStream.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/ProbeSet.h"
//
#include "file5/File5_File.h"
#include "file5/File5_Group.h"
#include "portability/affy-base-types.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <string>
#include <vector>
//


class EngineUtil {

public:

  /**
   * Read in a file of probe IDs
   * @param probeIds - vector of probeid_t to fill in, these are zero base
   * @param probeFileName - the file to load, file has 1 base probe IDs
   */
  static void readProbeFile(std::vector<probeid_t> &probeIds, const std::string& probeFileName);
  static void readProbeFile(std::vector< std::vector< probeid_t > > &probeIds, 
			    const std::string& probeFileName);
  static void readProbeFiles(std::vector<probeid_t> &probeIds, 
                             const std::vector<std::string>& probeFileNames);

  /**
   * Read in a file of probeset names or IDs
   * @param probeSets - vector of probeset names or IDs to fill in
   * @param probeSetFileName - the file to load, file has 1 base probe IDs
   */
  static void readProbeSetFile(std::vector<std::string> &probeSets, const std::string& probeSetFileName);
  static void readProbeSetFiles(std::vector<std::string> &probeSets,
                                const std::vector<std::string>& probeSetFileNames);

  /**
   * Read in the probe class file
   *
   * @param file   - The file to read
   * @param masks  - Vector of probe masks to fill in
   */
  static void readProbeClassFile(const std::string &file, 
          const int arraySize, 
          std::map<std::string, std::vector<bool> > &classMap);

  /**
   * Return true if the qType of chip is found in the list of possible
   * types that are ok. Note that this is a case sensitive operation and
   * is meant to be so.
   * 
   */
  static bool checkChipTypeOk(const std::string &qType, 
                              const std::vector<std::string> &possibleTypes);

  /**
   * Check chip types in TSV file against those in vector. Call Err::errAbort
   * if no match is found.
   */
  static void checkTsvFileChipType(const std::string &path, 
                                   const std::vector<std::string> &chipTypes,
                                   bool mustSee = false);

  /**
   * Check to make sure that there is at least one match between the two
   * vectors of chip types. Call Err::errAbort() if not.
   */
  static void checkChipTypeVectors(std::vector<std::string> &chipTypesOk,
                                   std::vector<std::string> &chipTypesCheck);

  /** 
   * Check to make sure all the cel files are compatible with the
   * chiptypes supported and number or probes on them. Treats the empty
   * string "" as a wildcard, so if chipTypes contains "" then any chip
   * type is accepted, but only if number of probes is correct. Calls
   * Err::errAbort() if something is wrong.
   * 
   * @param chipTypes - Names of chiptypes that are compatible.
   * @param probeCount - How many probes (aka features, ada cels) should be on array.
   * @param celFiles - Vector of cel file names to check.
   */
  static void checkCelChipTypes(std::vector<std::string> &chipTypes, 
                                uint32_t probeCount,
                                std::vector<std::string> &celFiles, 
                                uint32_t rows, 
                                uint32_t cols);

  /** 
   * @brief Get a list from a text file to select a subset of probe sets for
   * analysis.
   * 
   * @param fileName - Path to file for reading probe subset from.
   * @param psGroups - Groups to be filled in, must be deleted later
   * @param psNameLoadMap - Map of probesets names to be loaded which is
   * filled in.
   */
  static void makePSetNameSubset(const std::string &fileName, 
                                 std::vector<ProbeSetGroup *> &psGroups, 
                                 std::map<std::string, bool> &psNameLoadMap);

  /** 
   * Open up a pgf and clf file pair and read the valid chiptypes
   * (i.e. HuEx-1_0-st-v2, HuEx-1_0-st-v1, HuEx-1_0-st-ta1) and number
   * of probes expected from them.
   * 
   * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
   * @param probeCount - Number of probes that are on chip according to clf file.
   * @param pgfFile - Name of pgf file to be opened and read from.
   * @param clfFile - Name of clf file to be opened and read from.
   */
  static void getPgfChipType(std::vector<std::string>& chipTypes, 
                             colrow_t& rowCount, 
                             colrow_t& colCount,
                             int& probeCount,
                             // @todo: which isnt this signature like the others?
                             // int& probeSetCount,
                             const std::string& pgfFile, 
                             const std::string& clfFile);

  /** 
   * Open up a spf (simple probe format) file and read the valid chiptypes 
   * (i.e. HuEx-1_0-st-v2, HuEx-1_0-st-v1, HuEx-1_0-st-ta1) and number
   * of probes expected.
   * 
   * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
   * @param probeCount - Number of probes that are on chip according to spf file.
   * @param spfFile - Name of spf file to be opened and read from.
   */
  static void getSpfChipType(std::vector<std::string>& chipTypes, 
                             colrow_t& rowCount, 
                             colrow_t& colCount,
                             int& probeCount,
                             int& probeSetCount,
                             const std::string& spfFile);

  /** 
   * Open up a cdf file and get the chip type (i.e. 'Mapping250K_Sty') and number
   * of probes on the chip.
   * 
   * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
   * @param rowCount - Number of rows on the chip.
   * @param colCount - Number of columns on chip.
   * @param probeCount - Number of probes that are on chip according to clf file.
   * @param cdfFile - Name of cdf file to be opened and read from.
   */
  static void getCdfChipType(std::vector<std::string>& chipTypes, 
                             colrow_t& rowCount, 
                             colrow_t& colCount,
                             int& probeCount,
                             int& probeSetCount,
                             const std::string& cdfFile);


  /**
   * @brief Print out a header and then a series of self documenting
   * objects.
   * @param header - Descriptive header for this collection of self documenting objects.
   * @param doc - Collection of self documenting objects.
   */
  void static printSelfDocs(const char *header, std::vector<SelfDoc> doc);

  void static explainParameter(const std::string& query);

  void static getCelFiles(std::vector<std::string> &celFiles, BaseEngine *engine);

  void static getChpFiles(std::vector<std::string> &celFiles, BaseEngine *engine);

};

#endif /* _ENGINEUTIL_H_ */
