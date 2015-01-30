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
 * @file   IntensityReader.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:02:16 2005
 * 
 * @brief Interface for an object that reads data from disk into memory.
 */
#ifndef _INTENSITYREADER_H_
#define _INTENSITYREADER_H_

//
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 *  Abstract base class for a class that read files and fills in memory
 * representation and chip stream objects that need initialization.
 * 
 */
class IntensityReader {

public:

  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~IntensityReader() {}

  /** 
   * Set the names of the files to read data from.
   * @param fileNames - vector of file names to be opened.
   */
  void setFiles(const std::vector<std::string> &fileNames) {
    m_FileNames = fileNames;
  }

  /** 
   * return list of file names registered with the reader
   */
  std::vector<std::string> getFiles() {
    return m_FileNames;
  }

  /** 
   * @brief Push this chip stream onto the collection that will
   * receive data.
   * @param stream - Stream to be registered.
   */
  void registerStream(ChipStream *stream) {
    m_Streams.push_back(stream);
  }

  /** 
   * @brief Push an intensity mart onto the list of those that should
   * receive data.
   * @param iMart - intensity mart to be put on the list.
   */
  void registerIntensityMart(IntensityMart *iMart) {
    m_IntenMarts.push_back(iMart);
  }

  /** 
   * @brief Do the heavy lifting of reading data from the cel files
   * and passing to all the streams and intensity marts.
   * @return true on success or false on error.
   */
  virtual bool readFiles() = 0;

protected:
  /// Files that are being used.
  std::vector<std::string> m_FileNames;
  /// Streams we are initializing
  std::vector<ChipStream *> m_Streams;
  /// Intensity Marts we are initializing
  std::vector<IntensityMart *> m_IntenMarts;

  std::vector<double> m_Saturation;

};

#endif /* _INTENSITYREADER_H_ */
