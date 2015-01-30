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
 * @file   CelReader.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:11:45 2005
 * 
 * @brief Class for reading cel files and passing the data to memory
 * representation and analysis objects.
 */

#ifndef _CELREADER_H_
#define _CELREADER_H_

//
#include "chipstream/CelListener.h"
#include "chipstream/IdxGroup.h"
#include "chipstream/IntensityReader.h"
//
#include "util/Err.h"
#include "util/Verbose.h"
//
#include <vector>
//

/**
 *  Class for reading cel files and passing the data to memory
 * representation and analysis objects.
 */
class CelReader : public IntensityReader {

public:

  CelReader() : m_Size(0) {}

  /** 
   * @brief Do the heavy lifting of reading data from the cel files
   * and passing to all the streams and intensity marts.
   * @return true on success or false on error.
   */
  bool readFiles();

  /**
   * @brief Register a cel listener to see the cel file when it is
   * read
   */
  void registerCelListener(CelListener *reader) {
    m_CelListeners.push_back(reader);
  }

  /** 
   * @brief keep track how CEL channels are grouped
   *
   * @param cel_channels - lookup object that will be populated when
   * CELs are read
   */
  void setCelChannels(IdxGroup* cel_channels) {
    m_CelChannels = *cel_channels;
  }

  inline double GetSaturationValue (int idx)
  {
	  if (idx < 0 || idx >= m_Saturation.size ())
		  return -1;

	  return m_Saturation[idx];
  }

private:

  /// Vector of CelListeners to pass the cel file to
  std::vector<CelListener *> m_CelListeners;
  /// Number of features in the cel file
  int m_Size;
  /// object to track multi-CEL channels
  IdxGroup m_CelChannels;

};

#endif /* _CELREADER_H_ */
