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
// affy/sdk/file5/File5_Tsv.h ---
// 
// $Id: File5_Tsv.h,v 1.33 2009-10-28 19:06:41 csugne Exp $
// 
//
// @file   File5_Tsv.h
// @brief  Headers for the File5_Tsv classes.

/*
  We would like to store TsvFile like data in an HDF5
  formated file.  To accommodate the features of TsvFile in
  an HDF5 file, we keep each column in its own vector and
  use the curent "line" to figure out which level is being
  looked at and the offset into that set of tables.

  The column isnt in "tsv_lvl_idx", maps the current line to
  the correct level table and index into that table.

| LINE   LVL-IDX       | COL-0-0  COL-1-0    COL-2-0
+-----------------------+----------------------------
| 0       0-0          |foo
| 1       0-1          |bar
| 2       0-2          |baz
| 3       1-0          |        level1
| 4       1-1          |        level1-again
| 5       1-2          |        level1-once more
| 6       2-0          |                      lvl-2
| 7       2-1          |                      lvl-2

col-${level_}-${column_idx}

NOTES:

* The tsvFile shouldnt be a "File" it should be a "group" or "set"
For organizing File5TsvColumns

*/

#ifndef _FILE5_TSV_H_
#define _FILE5_TSV_H_

//
#include "file5/File5_File.h"
#include "file5/File5_Vector.h"
#include "file5/File5_types.h"
//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

#define TSV5_PREFIX "tsv-"

//class affx::File5_TsvColumn : public affx::File5_Vector {
//}

class affx::File5_TsvColumn : public affx::File5_Vector {
public:

  std::string m_col_name;
  int m_clvl;
  int m_cidx;

  // Current value for the column.
  int m_val_i;
  float m_val_f;
  double m_val_d;
  std::string m_val_string;

  //
  File5_TsvColumn();
  virtual ~File5_TsvColumn();

  //
  int init();
  int flush();
  int close();
  virtual char file5_kind_char();

  //
  std::string getColumnName();

  int set_c(char val);
  int set_s(short val);
  int set_i(int val);
  int set_f(float val);
  int set_d(double val);
  int set_string(const std::string& val);

  //
  int push_back_val();

  //
  int open(const std::string& name);
  int open(const std::string& name,
           const std::string& cname,
           int clvl,
           int cidx,
           affx::File5_dtype_t ctype,
           int flags);
  //
  int dump();
};

//////////

class affx::File5_Tsv : public affx::File5_Group {
  //
  typedef uint32_t linenum_t;
  
  //
public:
  // 
  std::vector<std::vector<affx::File5_TsvColumn*> > m_columns;
  std::vector<std::map<std::string,int> > m_column_name_map;

  //
  int m_lineNum; ///< The current line number

  affx::File5_Vector* m_line_clvl_idx;
  affx::File5_Vector* m_line_clvl;
  std::vector<int>    m_clvl_last_idx;
  
public:
  File5_Tsv();
  File5_Tsv(affx::File5_File* file5_ptr);
  virtual ~File5_Tsv();

  //
  int create(const std::string& name,int flags);
  int open(const std::string& name,int flags);

  int init();
  int flush();
  int close();
  virtual char file5_kind_char();

  /// This is an internal pointer - do not delete
  affx::File5_TsvColumn* getColumnPtr(int clvl,int cidx);
  affx::File5_TsvColumn* getColumnPtr(int clvl,const std::string& cname);
  //
  int getColumnIdx(int clvl,const std::string& cname);

  //
  affx::File5_dtype_t getColumnDtype(int clvl,int cidx);
  affx::File5_dtype_t getColumnDtype(int clvl,const std::string& cname);

  //
  int writeLevel(int clvl);
  int writeLevel(int clvl, int cidx);

  /// @brief     Copy the "format" of the file
  /// @param     tsv5     the File5_Tsv to copy the format from
  void copyFormat(affx::File5_Tsv& tsv5);

  //
  void openLineIndexes();

  // rewinds back to the start for reads.  (AKA: gotoLine(0))
  int rewind();
  // 
  int gotoLine(int line);
  //
  int nextLine();
  int nextLevel(int clvl);
  //
  int lineLevel();
  int lineNum();
  int getLineCount();

  // set all the columns to this size. (only works for single level Tsvs)
  void resize(size_t new_size);

  //
  int findFirst(int clvl,const std::string& cname,const std::string& val);
  int findFirst(int clvl,int cidx,const std::string& val);

  //
  void printfTabs(int cnt);
  int printfHeader();
  int printfRow();

  //
  int getLevelCount();
  int getColumnCount(int clvl);
  int getColumnName(int clvl,int cidx,std::string* cname);

  //
  int getAsStr(int clvl,int cidx,std::string* cname);

  //
  int get(int clvl,int cidx,char* val);
  int get(int clvl,int cidx,int* val);
  int get(int clvl,int cidx,float* val);
  int get(int clvl,int cidx,double* val);
  int get(int clvl,int cidx,std::string* val);

  int getLine(int clvl,int cidx, int linxIx, float* val);
  //
  int get(int clvl,const std::string& cidx,char* val);
  int get(int clvl,const std::string& cidx,int* val);
  int get(int clvl,const std::string& cidx,float* val);
  int get(int clvl,const std::string& cidx,double* val);
  int get(int clvl,const std::string& cidx,std::string* val);
  //
  int set_c(int clvl,int cidx,char val);
  int set_i(int clvl,int cidx,int val);
  int set_f(int clvl,int cidx,float val);
  int set_d(int clvl,int cidx,double val);
  int set_string(int clvl,int cidx,const std::string& val);
  //
  int set_c(int clvl,const std::string& cidx,char val);
  int set_i(int clvl,const std::string& cidx,int val);
  int set_f(int clvl,const std::string& cidx,float val);
  int set_d(int clvl,const std::string& cidx,double val);
  int set_string(int clvl,const std::string& cidx,const std::string& val);

  // @todo make a function that copies from one column to another
  // When defining a string must specify the maximum possible string size or -1 for unlimited string length
  int defineColumn(const int clvl,
                   const int cidx,
                   const std::string& cname,
                   const File5_dtype_t ctype,
                   int str_size);
  // without the string size.
  int defineColumn(const int clvl,
                   const int cidx,
                   const std::string& cname,
                   const File5_dtype_t ctype);
  // no need for the dtype
  int defineStringColumn(const int clvl,
                         const int cidx,
                         const std::string& cname,
                         int str_size);
  //
  int registerColumn(affx::File5_TsvColumn* col,
                     const int clvl,
                     const int cidx,
                     const std::string& cname);
  //
  int dump();

  //
  static int getFileTsvLineCount(const std::string& file_name, 
                                 const std::string& tsv_name);

  static int getFileTsvLineCount(const std::string& file_name,
								 const std::string& group_name,
                                 const std::string& tsv_name);
};

#endif // _FILE5_TSV_H_
