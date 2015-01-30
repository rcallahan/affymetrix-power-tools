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
//
// sdk/bboard/dao/DaoUtil.h ---
//
// $Id: DaoUtil.h,v 1.2 2009-11-05 20:41:39 harley Exp $
//

#ifndef _DAOUTIL_H_
#define _DAOUTIL_H_

//
#include "bboard/Apt2Types.h"
#include "bboard/BboardTypes.h"
#include "bboard/dao/DaoTypes.h"

// io obj defs
#include "bboard/Bboard.h"
#include "util/PgOptions.h"

// these are collected here to avoid littering the Dao_Table interface.

#define DAOUTIL_DEF_OBJ(_type)                                      \
  static AptErr_t readFromFile(const FsPath& dao_p,_type* obj);    \
  static AptErr_t readFromFile(const std::string& path,_type* obj); \
  static AptErr_t readFromTable(Dao_Table* dao_t,_type* bb);        \
  static AptErr_t writeToFile(const FsPath& dao_p,_type* obj);     \
  static AptErr_t writeToFile(const std::string& path,_type* obj);  \
  static AptErr_t writeToTable(Dao_Table* dao_t,_type* bb);

class DaoUtil {
public:
  //
  DAOUTIL_DEF_OBJ(Bboard);
  DAOUTIL_DEF_OBJ(PgOptions);
};

#endif // _DAOUTIL_H_
