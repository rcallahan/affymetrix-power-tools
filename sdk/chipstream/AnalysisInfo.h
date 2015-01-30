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
 * @file   AnalysisInfo.h
 * @author Chuck Sugnet
 * @date   Wed Aug  2 14:05:13 2006
 *
 * @brief  Information about a particular analysis run.
 */

#ifndef _ANALYSISINFO_H_
#define _ANALYSISINFO_H_

//
#include "file/CDFFileData.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cstring>
#include <map>
#include <string>

/**
 * Information about a particular analysis run.
 */
class AnalysisInfo {
public:

  /** Constructor. */
  AnalysisInfo() {
    m_MaxPsNameLength = 0;
    m_NumExpression = 0;
    m_NumGenotyping = 0;
    m_NumReporting = 0;
    m_NumRows = 0;
    m_NumCols = 0;
    m_NumProbeSets = 0;
    m_ProbeSetType = affxcdf::UnknownProbeSetType;
    m_ChipType = "";
    m_ProgID = "";
    m_ProgramName = "";
    m_ProgramVersion = "";
    m_ProgramCompany = "";
    m_AlgName = "";
    m_AlgVersion = "";
    m_AnalysisGuid = "";
    m_ExecGuid = "";
  }

  /**
   * Add a key/value parameter pair to our inforomation. Keys must be unique strings.
   *
   * @param key - name of key to be stored.
   * @param value - matching value to be associated with key.
   */
  void addParam(const std::string &key, const std::string &value) {
    std::map<std::string, uint32_t>::iterator mapIx = m_KeysIndex.find(key);
    if(mapIx != m_KeysIndex.end()) {
      Err::errAbort("AnalysisInfo::addParam() - Key '" + key + "' has already been used.");
    }
    m_KeysIndex[key] = (uint32_t) m_ParamNames.size();
    m_ParamNames.push_back(key);
    m_ParamValues.push_back(value);
  }

  /**
   * Add client supplied key/value parameter pair to our information. Keys must be unique strings.
   *
   * @param key - name of key to be stored.
   * @param value - matching value to be associated with key.
   */
  void addClientMetaInfo(const std::string &key, const std::string &value) {
    std::map<std::string, uint32_t>::iterator mapIx = m_ClientInfoKeysIndex.find(key);
    if(mapIx != m_ClientInfoKeysIndex.end()) {
      Err::errAbort("AnalysisInfo::addClientMetaInfo() - Key '" + key + "' has already been used.");
    }
    m_ClientInfoKeysIndex[key] = (uint32_t) m_ClientInfoNames.size();
    m_ClientInfoNames.push_back(key);
    m_ClientInfoValues.push_back(value);
  }


  /**
    * Turn a Info set of parameters into a text representation of:
    *    key1=value1
    *    key2=value2
    *    key3=value3
    *
    * @param sep - optional separatar string, default newline
    *
    * @return - string representation of Info set.
    */
  std::string toString(std::string sep = "\n") {
      std::string str;
     str += "apt-info-chiptype=" + m_ChipType + sep;
     str += "apt-info-program-id=" + m_ProgID + sep;
     str += "apt-info-program-name=" + m_ProgramName + sep;
     str += "apt-info-program-version=" + m_ProgramVersion + sep;
     str += "apt-info-alg-name=" + m_AlgName + sep;
     str += "apt-info-alg-version=" + m_AlgVersion + sep;
     str += "apt-info-num-rows=" + ToStr(m_NumRows) + sep;
     str += "apt-info-num-cols=" + ToStr(m_NumCols) + sep;
     str += "apt-info-num-probesets=" + ToStr(m_NumProbeSets) + sep;
     for(uint32_t i = 0; i < m_ParamNames.size(); i++) {
        str += ("apt-info-" + m_ParamNames[i] + "=" + m_ParamValues[i] + sep);
     }
     return str;
  }

  int getNumCols() {
    return m_NumCols;
  }

  /// Length of the longest probeset name
  int m_MaxPsNameLength;
  /// Name of the analysis
  std::string m_AnalysisName;
  /// Number of expression probesets
  int m_NumExpression;
  /// Number of genotyping expression
  int m_NumGenotyping;
  /// How many snps to report in CHP file.
  int m_NumReporting;
  /// How many rows are on the grid of the chip.
  int m_NumRows;
  /// How many columns are the grid of the chip.
  int m_NumCols;
  /// How many probesets are on the chip and in cdf file. This is
  /// critical as probesets must be the exact same in number and order
  /// as in the cdf file.
  int m_NumProbeSets;
  /// What type of chp file is this? For our purposes here it is genotyping
  affxcdf::GeneChipProbeSetType m_ProbeSetType;
  /// What chip is this?
  std::string m_ChipType;
  /// What is our id for GCOS?
  std::string m_ProgID;
  /// What is the name of the program
  std::string m_ProgramName;
  /// What program created this chp file.
  std::string m_ProgramVersion;
  /// What company created this chp file
  std::string m_ProgramCompany;
  /// What algorithm are we reporting
  std::string m_AlgName;
  /// What version of the algorithm is being used.
  std::string m_AlgVersion;
  /// GUID (global unique ID) associated with a specific analysis
  std::string m_AnalysisGuid;
  /// GUID associated with a specific process
  std::string m_ExecGuid;
  /// These next two vectors are matched pairs, one contains the
  /// keys and the other the values. Would use a map, but want to be
  /// guaranteed of the order.
  std::vector<std::string> m_ParamNames;
  std::vector<std::string> m_ParamValues;

  /// A mapping of key names to their index in m_ParamNames.
  std::map<std::string, uint32_t> m_KeysIndex;

  /// Two vectors for client supplied meta-data. Matched pairs of vectors, one contains the
  /// keys and the other the values. Would use a map, but want to be
  /// guaranteed of the order.
  std::vector<std::string> m_ClientInfoNames;
  std::vector<std::string> m_ClientInfoValues;

  /// A mapping of key names to their index in m_ParamNames.
  std::map<std::string, uint32_t> m_ClientInfoKeysIndex;
  
  /// Vector of probeset names and groups as they should appear in the chp file.
  /// this memory is owned elsewhere as it is relatively large in size.
  /// Only populated if there is CHP output
  std::vector<const char *> m_ProbesetNames;

  /// Vector of probeset names and groups as they should appear in the chp file.
  /// this memory is owned elsewhere as it is relatively large in size.
  /// Only populated if there is CHP output
  std::vector<const char *> m_ProbesetDisplayNames;

};

#endif /* _ANALYSISINFO_H_ */
