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
// CalvinLite/CL_File.h ---
//
// $Id: CL_File.h,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#ifndef _CALVINLITE_FILE_H_
#define _CALVINLITE_FILE_H_

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_Object.h"
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_DataGroup.h"
//
#if CL_WITH_TSVFILE==1
#include "file/TsvFile/TsvFile.h"
#endif
//
#include <fstream>
#include <string>
#include <utility>
#include <vector>

extern int CL_File_debug_flags;

typedef std::vector<std::pair<std::string,std::string> >  CL_changevec_t;

class CL_File : public CL_Object {
public:
  //
  std::string m_filename;
  std::fstream m_fstream;
  //
  int m_opt_debug;
  int m_opt_verbose;
  //
  int m_file_end_fpos;
  //
  int m_hdr_magic;
  int m_hdr_version;
  //
  int m_hdr_datagroup_cnt;
  int m_hdr_datagroup_pos;
  CL_Gdh* m_hdr_gdh;
  //
  CL_DataGroup* m_bad_datagroup;
  CL_DataSet*   m_bad_dataset;
  CL_Gdh* m_bad_gdh;

  //
  std::vector<std::string> m_tsvexport_filenames;

  //
  std::vector<CL_DataGroup*> m_datagroups;

  //
  static CL_Err_t isCalvinFormat(const std::string& pathname);

  //
  CL_File();
  ~CL_File();
  //
  void init();
  void clear();
  CL_Err_t dump();
  CL_Err_t dumpSegs();

  /// shortcut error setter, supplies our filename as the third arg.
  CL_Err_t setErr(CL_Err_t err_num,const std::string& err_msg);

  std::string getFilename() const;

  //
  bool f_open_r();
  bool f_open_w();
  bool f_seekg(int pos);
  bool f_seekp(int pos);
  int f_tellg();
  int f_tellp();

  //
  CL_Err_t read_1B(int& val);
  CL_Err_t read_4B(int& val);

  //
  CL_Err_t read_4B_Nstring(CL_string& str);
  CL_Err_t read_4B_Wstring(CL_string& str);

  //
  CL_Err_t read_4B_Blob(CL_string& str);
  CL_Err_t read_Nstring_len(CL_string& str,int byte_len);
  CL_Err_t read_Wstring_len(CL_string& str,int byte_len);

  //
  CL_Err_t open(const std::string& filename);
  CL_Err_t close();

  //
  static CL_File* openFile(const std::string& filename);
  static CL_File* openFileNoErr(const std::string& filename);

  //
  CL_Gdh* getGdh();
  CL_Gdh* getGdhByNumber(const std::string& path);
  CL_Gdh* getGdhByPath(std::vector<std::string>& path_vec);
  // @todo rename to "defineDataGroup"
  CL_DataGroup* newDataGroup(const std::string& dg_name);

  //
  CL_Err_t read_File(const std::string& filename,CL_ReadOpt_t readopt=CL_READOPT_ALL);
  CL_Err_t read_Gdh(CL_Gdh* gdh);
  CL_Err_t read_Param(CL_Param* param);
  CL_Err_t read_DataGroup(CL_DataGroup* dgroup);
  CL_Err_t read_DataSet(CL_DataSet* dset);
  CL_Err_t read_DataSetCol(CL_DataSetCol* dcol);

  //
  int bytesize_Nstr(CL_string& str);
  int bytesize_Wstr(CL_string& str);
  int bytesize_ParamVec(std::vector<CL_Param>& params);
  int bytesize_Param(CL_Param* param);
  int bytesize_File(CL_File* dfile);
  int bytesize_Gdh(CL_Gdh* gdh);
  int bytesize_DataGroup(CL_DataGroup* dgroup);
  int bytesize_DataSet(CL_DataSet* dset);
  int bytesize_DataSetCol(CL_DataSetCol* dcol);

  //
  int getDataGroupCount() const;
  int getDataSetCount(int dg_idx) const;
  void chainDataGroups();

  //
  CL_DataGroup* getDataGroup(int dg_idx);
  CL_DataGroup* getDataGroup(const std::string& dg_name);
  //
  CL_DataSet* getDataSet(int dg_idx,int ds_idx);
  CL_DataSet* getDataSet(const std::string& dg_name,const std::string& ds_name);

  //
  int computeSizesAndOffsets();
  int roundToPage(int val);
  static int roundTo1KB(int);

  //
  CL_Err_t write_Padding(int cnt,int val);
  CL_Err_t write_Zeros(int cnt);
  CL_Err_t write_1B(int val);
  CL_Err_t write_4B(int val);
  CL_Err_t write_4B_Nstring(const std::string& val);
  CL_Err_t write_4B_Nstring(CL_string& val);
  CL_Err_t write_4B_Wstring(const std::string& val);
  CL_Err_t write_4B_Wstring(CL_string& val);
  CL_Err_t write_4B_Blob(CL_string& val);
  //
  CL_Err_t write(const std::string& filename);
  CL_Err_t write_File(const std::string& filename);
  CL_Err_t write_File(const std::string& filename,bool forceWrite);
  CL_Err_t write_FileHdr(bool forceWrite);
  CL_Err_t write_Gdh(CL_Gdh* gdh,bool forceWrite);
  CL_Err_t write_Param(CL_Param* param);
  CL_Err_t write_DataGroup(CL_DataGroup* dg,bool forceWrite);
  CL_Err_t write_DataSet(CL_DataSet* ds,bool forceWrite);
  CL_Err_t write_DataSetColDef(CL_DataSetCol* dcol);

  //
  bool requiresRewrite();

  //
  static CL_Err_t copyFile(const std::string& from, const std::string& to);

  //
  static CL_Err_t dumpFile(const std::string& filename);
  static CL_Err_t dumpFileSegs(const std::string& filename);
  //
  static CL_Err_t changeFile(const std::string& filename,const CL_changevec_t& change_vec,bool substr);
  static CL_Err_t changeFileGdh(CL_Gdh* gdh,const std::string& from,const std::string& to,bool substr,int* cnt);
  static CL_Err_t changeFileGdhValue(const std::string& key,
                                     CL_string* clstr,
                                     const std::string& from,
                                     const std::string& to,
                                     bool substr,int* cnt);
  //
  static CL_Err_t printFileParamsAll(const std::string& filename);
  static CL_Err_t printFileParamsAll(const CL_File* dfile);
  static CL_Err_t printFileParams(const std::string& filename,
                                  const std::vector<std::string>& key_vec);
  static CL_Err_t printFileParams(const CL_File*,
                                  const std::vector<std::string>& key_vec);

#if CL_WITH_TSVFILE==1
  //
  static CL_Err_t ExportToTsv(const std::string& filename,int dg_idx,int ds_idx,
                              const std::string& outdir);
  //
  static CL_Err_t ExportGdhToTsv(const std::string& file_name,const std::string& outdir);

  //
  CL_Err_t tsvExport(int dg_idx,int ds_idx,
                     const std::string& outdir);
  CL_Err_t tsvExportDataSet(CL_DataSet* ds,int dg_idx,int ds_idx,
                            const std::string& tsv_name);

  //
  CL_Err_t tsvExportGdh(CL_Gdh* gdh,const std::string& tsv_name);
  CL_Err_t tsvExportGdh(CL_Gdh* gdh,const std::string& gdh_name,affx::TsvFile* tsv);
  //
  CL_Err_t tsvExportGdh_Param(affx::TsvFile* tsv,
                         CL_Param* param);
  CL_Err_t tsvExportGdh_Param(affx::TsvFile* tsv,
                         const std::string& key_str,
                         const std::string& val_str);
  CL_Err_t tsvExportGdh_Param(affx::TsvFile* tsv,
                         const std::string& key_str,
                         CL_TypeCode_t val_type_code,
                         int byte_len,
                         const std::string& val_str,
                         const std::string& val_raw);
  //
  static CL_Err_t tsvImport(const std::string& file_in,
                            const std::string& file_out,
                            const std::vector<std::string>& tsvfiles);
  static CL_Err_t tsvImportScript(const std::string& script_file);
  
  //
  CL_Err_t tsvImport(const std::string& tsv_name);
  //
  CL_Err_t tsvImportData(const std::string& tsv_name);
  CL_Err_t tsvImportDataCopy(CL_DataSet* dset,affx::TsvFile* tsv);
  //
  CL_Err_t tsvImportGdh(const std::string& tsv_name);

#endif

};

#endif
