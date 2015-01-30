////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

//
// affy/sdk/file5/File5_File.h ---
//
// $Id: File5_File.h,v 1.14 2009-09-18 03:37:27 mspald Exp $
//

#ifndef _FILE5_FILE_H_
#define _FILE5_FILE_H_

//
#include "file5/File5_Group.h"
#include "file5/File5_Object.h"
#include "file5/File5_types.h"
//
#include <map>
#include <set>
#include <vector>
//

//
class affx::File5_File : public affx::File5_Group {
public:
  std::string m_file_name;

  //
  File5_File();
  virtual ~File5_File();

  //
  void setFilename(const std::string& file_name);

  //
  bool is_open();

  //
  File5_return_t open(const std::string& file_name,int flags);
  File5_return_t init();
  File5_return_t flush();
  File5_return_t dump();
  File5_return_t dump1();
  File5_return_t close();
  //
  virtual char file5_kind_char();

  static bool isHdf5file(const std::string& file_name);

  static bool equivalent(	const std::string& strFileName1, 
                                const std::string& strFileName2, 
                                const std::string& strGroupName, 
                                const std::string& strTsvName, 
                                std::set<std::string>& setIgnore, 
                                double dEpsilon = 0.0001, 
                                double dCorrelationCutoff = 1.0,
                                bool bAllowNegation=false,
                                bool flagNaNNumDiff=false);

};

#endif // _FILE5_FILE_H_
