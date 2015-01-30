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


#ifndef _CELFileConversionOptions_HEADER_
#define _CELFileConversionOptions_HEADER_

#include <cstring>
#include <string>

/*! file CELFileConversionOptions.h This file contains a class for options when converting CEL files. */

namespace affymetrix_cel_converter
{

/*! This class holds options for converting CEL files. */
class CELFileConversionOptions
{
public:

  /*! Constructor */
  CELFileConversionOptions();

  /*! The chip type to write out */
  const char *m_ChipType; // we do not own this memory, null == no change

  /*! The file name to put in the header */
  const char *m_DATFileName; // we do not own this memory, null == no change

  static std::string newDatName(const std::string &datHeader, const std::string &newName);

};

}

#endif
