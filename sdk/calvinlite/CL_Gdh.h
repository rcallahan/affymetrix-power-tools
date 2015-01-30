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
// CalvinLite/CL_Gdh.h ---
//
// $Id: CL_Gdh.h,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#ifndef _CALVINLITE_GENERICDATAHEADER_H_
#define _CALVINLITE_GENERICDATAHEADER_H_

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_ObjectWithParams.h"
#include "calvinlite/CL_Param.h"
#include "calvinlite/CL_string.h"
//
#include <set>
#include <string>
#include <vector>

class CL_Gdh : public CL_ObjectWithParams {
public:

  // @todo Perhaps this should be written not to have "m_"
  //       but store the values as params with special names?

  // in gdh file order
  CL_string m_datatype_str;
  CL_string m_uuid_str;
  CL_string m_datetime_str;
  CL_string m_locale_str;

  // we use the word "leaves" as we want to use the words
  // "parent and child" for containership in CalvinLite
  // since the Calvin terminlogy is mixed up.

  std::vector<CL_Gdh*> m_gdh_leaves;

  //
  CL_Gdh();
  ~CL_Gdh();
  //
  void clear();
  void defaultValues();
  //
  void setDataType(const std::string& val);
  void setUuid(const std::string& val);
  void setDateTime(const std::string& val);
  void setLocale(const std::string& val);
  //
  CL_Gdh* newLeafGdh();
  //
  int getLeafCount();
  CL_Gdh* getLeafGdh(int i);

  //
  CL_Gdh* getGdhByPath(std::vector<std::string>& path_vec);

  //
  void dump();
  void dump(const std::string& prefix);
  //
  void print(const std::string& prefix);
  void printMatchingKeys(const std::vector<std::string>& key_vec,const std::string& prefix);
  void printMatchingKeys(const std::set<std::string>& key_set,const std::string& prefix);
};

#endif
