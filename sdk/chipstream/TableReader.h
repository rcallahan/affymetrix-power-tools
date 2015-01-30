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
 * @file   TableReader.h
 * @author Chuck Sugnet
 * @date   Mon Nov  7 11:15:39 2005
 * 
 * @brief IntensityMart object for reading intensity files for text
 * files.
 */

#ifndef _TABLEREADER_H_
#define _TABLEREADER_H_

//
#include "chipstream/IntensityReader.h"
//
#include "util/Convert.h"
#include "util/Err.h"
//

/**
 * @brief IntensityMart object for reading intensity files for text,
 * useful for testing via text files.
 */
class TableReader : public IntensityReader {

public:

  /** 
   * @brief Constructor
   * @param clfProbeCount - Number of probes in file.
   * @param columnHeaderNames - Does the table have column header names to be read.
   */
  TableReader(int clfProbeCount, bool columnHeaderNames=false) : 
    m_ClfProbeCount(clfProbeCount), m_RowColNames(columnHeaderNames), m_Read(false), m_Passed(false) {}

  /** 
   * @brief Read the file that contains the intensity data.
   * @return true if opened sucessfully, false otherwise.
   */
  bool readFiles();

private: 
  /// How many probes 
  int m_ClfProbeCount;
  /// Are ther row and column names?
  bool m_RowColNames;
  /// Are we done reading
  bool m_Read;
  /// Are we done initializing stream
  bool m_Passed;
};

#endif /* _TABLEREADER_H_ */
