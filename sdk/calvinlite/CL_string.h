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
// CalvinLite/CL_string.h ---
// 
//  $Id: CL_string.h,v 1.1 2009-10-27 16:53:53 harley Exp $
// 

#ifndef _CALVINLITE_STRING_H_
#define _CALVINLITE_STRING_H_

//
#include "calvinlite/CalvinLite.h"
//
#include <string>

enum CL_string_fmt_t {
  CL_STRING_NULL,
  CL_STRING_BLOB,
  CL_STRING_NSTR,
  CL_STRING_WSTR
};

class CL_string {
public:
  // the narrow 8b string
  mutable std::string m_nstr;
  mutable bool m_nstr_valid;
  // the wide 16b string
  mutable std::string m_wstr;
  mutable bool m_wstr_valid;
  
  // what the format should be in the calvin file
  CL_string_fmt_t m_calvin_fmt;

  //
  CL_string();
  void clear();
  void ensureNstrValid();
  void ensureWstrValid();

  //
  static void NstrToWstr(const std::string& nin,std::string& wout);
  static void WstrToNstr(const std::string& win,std::string& nout);

  //
  std::string Nstr();
  std::string Wstr();
  const char* c_str();
  const char* blob_str();

  void setNstr(const std::string& val);
  void setWstr(const std::string& val);

  char* resizeNstr(int size);
  char* resizeWstr(int size);
  char* resizeBlob(int size);

  // in bytes not size.
  int byteLen();

  void dump() const;
};

#endif
