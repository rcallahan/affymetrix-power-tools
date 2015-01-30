////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License
// (version 2.1) as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/*! file CELFileConversionOptions.cpp This file contains a class for options when converting CEL files. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"

namespace affymetrix_cel_converter
{

  /*! Constructor */
CELFileConversionOptions::CELFileConversionOptions()
{
  m_ChipType = NULL;
  m_DATFileName = NULL;
}

std::string CELFileConversionOptions::newDatName(const std::string &datHeader, const std::string &newName)
{
  std::string s;

  int start = datHeader.find("]  ", 0);
  if (start == std::string::npos) {
    std::string empty;
    return empty;
  } else {
    start += 3;
  }

  int stop = datHeader.find(":CLS", 0);
  if (stop == std::string::npos) {
    std::string empty;
    return empty;
  } else {
    stop -= 1;
  }

  s = datHeader.substr(0, start) + newName + datHeader.substr(stop + 1);

  return s;
}

};

