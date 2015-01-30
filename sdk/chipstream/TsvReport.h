////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
// affy/sdk/chipstream/TsvReport.h ---
// 
// $Id: TsvReport.h,v 1.12 2009-05-28 22:02:52 awilli Exp $
// 

/// @file   chipstream/TsvReport.h
/// @brief  Header for the TsvReports.

/// *sigh* So this isnt a "report" at all. 
/// Reports should have the signature "prepare/report/finish" which this does not.
/// However, it is used by the reports, hence the name.
/// It is really an aggeration or union of the TsvFile and A5/Tsv.
/// Once the format is set, it directs the calls to the underling file IO.
/// Maybe "ReportIO" would be the better name?

#ifndef _TSVREPORT_H_
#define _TSVREPORT_H_

//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "util/Fs.h"
//
#include <utility>

// suggested sizes for report columns.
// TsvFile ignores them, A5-Tsv uses them. 
// (-1 for unlimited size, but this will chew up disk space.)
#define TSVREPORT_PROBESET_STRLEN  30
#define TSVREPORT_CELFILE_STRLEN  200
#define TSVREPORT_METRIC_STRLEN   200

/// The string to use for per-report guids.
#define TSVREPORT_GUID_HDR "guid"

namespace affx {
  class TsvReport;
}

/// @brief     Unifies the presentation of file/TsvFile and file
/// @return    
class affx::TsvReport {

public:
  enum TsvReportFmt_t {
    FMT_UNSET = 0,
    FMT_TSV,
    FMT_A5,
  };
  typedef std::pair<std::string,std::string> keyval_t;

  /// 
  bool m_is_open;

  /// Has the guid been set for this report?
  int m_has_guid_header;

  /// selected output format for this report.
  TsvReportFmt_t m_format;
  /// Where the report will go.
protected:
  std::string m_dir_path;
  std::string m_file_prefix;
  std::string m_file_name;
  std::string m_group_name;
  int m_use_default_suffix;

  ///
  affx::TsvFile* m_tsv;
  /// How many digits after the decimal place do we write.
  int m_precision;
  
  //
  affx::File5_Group* m_a5_shared_group;
  affx::File5_Group* m_a5_group;
  affx::File5_File*  m_a5_file;
  affx::File5_Tsv*   m_a5_tsv;

  //
  int m_is_header_buffer;
  std::vector<keyval_t> m_header_buffer;

public:

  //
  TsvReport();
  ~TsvReport();

  //
  void init();
  void close();

  //
  bool is_open();

  //
  void setDirPath(const std::string& dirpath);
  std::string getFileprefix();
  void setFileprefix(const std::string& fileprefix);
  void setFilename(const std::string& fileName);
  void setFormat(affx::TsvReport::TsvReportFmt_t format);
  void setGroupname(const std::string& groupname);
  void setIsHeaderBuffer(int val);
  void setUseDefaultSuffix(int val);

  //
  void setPrecision(int precision);
  void setPrecision(int clvl,int cidx,int precision);

  //
  std::string getFilePath();
  std::string getFmtSuffix();

  //
  void setA5SharedGroup(affx::File5_Group* shared_group);
  static affx::File5_File* openA5File(const std::string& filename, bool replace=true);
  static void closeA5File(affx::File5_File*& a5_file);

  //
  int get(int clvl,int cidx,int* val);
  int get(int clvl,int cidx,float* val);
  int get(int clvl,int cidx,double* val);
  int get(int clvl,int cidx,std::string* val);
  //
  int set_c(int clvl,int cidx,char val);
  int set_i(int clvl,int cidx,int val);
  int set_f(int clvl,int cidx,float val);
  int set_d(int clvl,int cidx,double val);
  int set_string(int clvl,int cidx,const std::string& val);

  //
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
  int getColumnCount(int clvl);
  //
  void defineColumns(const std::vector<std::string>& colNames,File5_dtype_t ctype,int str_size);

  //
  int addHeader(const std::string& key, const std::string& val);
  int addHeaderComment(const std::string& comment);
  int addHeaderComments(const std::vector<std::string>& comments);

  //
  void addHeadersFrom(const std::vector<keyval_t>& header_vec);
  //
  int writeLevel(int lvl);

  // 
  int writeTsv_v1(const std::string& filename,bool raise_on_err);
  int writeTsv_v1(const std::string& filename) {
    return writeTsv_v1(filename,true);
  };
  int writeTsv_v1() {
    return writeTsv_v1(getFilePath(),true); 
  };
  int writeTsv_v2(const std::string& filename,bool raise_on_err);
  int writeTsv_v2(const std::string& filename) {
    return writeTsv_v2(filename,true);
  }
  int writeTsv_v2() {
    return writeTsv_v2(getFilePath(),true); 
  };
  //
  void copyOptionsTo(affx::TsvReport& tsv);
  //
  void flushHeaderBuffer();
  //
  void ensureGuidHeader();
};

//
#endif
