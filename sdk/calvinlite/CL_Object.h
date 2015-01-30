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
// CalvinLite/CL_Object.h ---
// 
// $Id: CL_Object.h,v 1.1 2009-10-27 16:53:53 harley Exp $
// 

#ifndef _CALVINLITE_OBJECT_H_
#define _CALVINLITE_OBJECT_H_

//
#include "calvinlite/CalvinLite.h"
//
#include <string>

class CL_Object {
public:
  /// All calvin objects have a parent file. (It might be themselves.)
  CL_File* m_parent_file;
  //
  int m_fpos_start;
  int m_fpos_size;
  int m_fpos_end;
  //
  int m_orig_fpos_start;
  int m_orig_fpos_size;
  int m_orig_fpos_end;
  //
  int m_data_dirty;
  //
  int m_bytesize;
  int m_bytesize_pad;
  int m_bytesize_dirty;

  /// the numeric error message.
  CL_Err_t m_err_num;
  /// the text of the error message
  std::string m_err_msg;
  /// the path of the file which cause the error.
  //  Included to be symmetric with CalvinLiteException
  std::string m_err_path;
  
  bool m_err_should_throw;

  //
  CL_Object();
  void init();
  void clear();
  void dump();

  //
  // should have "raiseErr" as well to cause a setErr and throw.
  CL_Err_t setErr(CL_Err_t err_num,
                  const std::string& err_msg,
                  const std::string& err_path);
  CL_Err_t errNum();
  std::string errMsg();
  std::string errPath();
  void printErrMsg();
  void clearErr();
  /// setting this to false will prevent CalvinLiteException from being thrown on error.
  /// defaults to true.
  bool throwOnErr(bool val);

  //
  void setParentFile(CL_File* cl_file);
  CL_File* getParentFile();

  //
  void setDataDirty();
  void setDataDirty(bool val);
  bool isDataDirty();

  //
  int getBytesize();
  int setBytesize(int val);
  bool isBytesizeDirty();

  //
  void clearFpos();
  // inclusive start
  int getFposStart() const;
  int setFposStart(int val);
  int getFposSize() const;
  int setFposSize(int val);
  // exclusive end
  int getFposEnd() const;
  int setFposEnd(int val);

  //
  void clearOrigFpos();
  //
  int getOrigFposStart() const;
  int setOrigFposStart(int val);
  int getOrigFposSize() const;
  int setOrigFposSize(int val);
  int getOrigFposEnd() const;
  int setOrigFposEnd(int val);
  //
  bool changedFpos() const;
};

#endif
