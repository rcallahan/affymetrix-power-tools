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
// CalvinLite/CL_File.cpp ---
//
// $Id: CL_File.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

//
#include "calvinlite/CL_File.h"
//
#include "calvinlite/CL_private.h"
#include "calvinlite/CL_util.h"
// for snprintf
#include "portability/affy-system-api.h"
#include "util/Util.h"

//
#include <assert.h>
#include <ios>
#include <stdio.h>
#include <string.h>
#include <stdio.h>

//
int CL_File_debug_flags=0;

//
#define CL_SIZEOF_SCALARBUF 8

//
//#define CL_PRINT_POS(_lbl) { int pos=m_fstream.tellg(); printf("POS: %s: %08x \n",_lbl,pos); fflush(NULL); }
//#define CL_PRINT_DAT_D(_dat) { printf("    => %d      %08x\n",_dat,_dat); }
//#define CL_PRINT_DAT_S(_dat) { printf("    => '%s'\n",_dat.c_str()); }

#ifndef CL_PRINT_POS
#define CL_PRINT_POS(_lbl)
#endif
#ifndef CL_PRINT_DAT_D
#define CL_PRINT_DAT_D(_x)
#endif
#ifndef CL_PRINT_DAT_S
#define CL_PRINT_DAT_S(_x)
#endif

// for when this is optional
//#define CL_DEBUG_BYTESIZE_ON 1
//#ifdef CL_DEBUG_BYTESIZE_ON

#define CL_DEBUG_BYTESIZE_START(_lbl)           \
  int cl_debug_bytesize_fstart=f_tellg();       \
  int cl_debug_bytesize_fend=0;                 \
  int cl_debug_bytesize_len;                    \
  if (CL_File_debug_flags==1) {                 \
    printf("=== %s: start\n",_lbl);             \
  }

#define CL_DEBUG_BYTESIZE(_lbl,_expected_len)                           \
  if (CL_File_debug_flags==1) {                                         \
    cl_debug_bytesize_fend=f_tellg();                                   \
    cl_debug_bytesize_len=cl_debug_bytesize_fend-cl_debug_bytesize_fstart; \
    printf("=== %s => S=%6d  E=%6d  Len:  Expect=%6d  Actual=%6d (%+d)\n", \
           _lbl,                                                        \
           cl_debug_bytesize_fstart,cl_debug_bytesize_fend,             \
           _expected_len,                                               \
           cl_debug_bytesize_len,                                       \
           cl_debug_bytesize_len-_expected_len);                        \
    fflush(NULL);                                                       \
    m_fstream.flush();                                                  \
    assert(cl_debug_bytesize_len==_expected_len);                       \
  }

#define CL_DEBUG_BYTESIZE_FPOSEND(_lbl,_end)            \
  if (CL_File_debug_flags==1) {                         \
    cl_debug_bytesize_fend=f_tellg();                   \
    printf("== %s == expected= %d actual= %d  (%+d)\n", \
           _lbl,                                        \
           _end,cl_debug_bytesize_fend,                 \
           cl_debug_bytesize_fend-_end);                \
    assert(cl_debug_bytesize_fend==_end);               \
  }

// #else
// #define CL_DEBUG_BYTESIZE_START(_lbl)
// #define CL_DEBUG_BYTESIZE(_lbl,_expected_len)
// #define CL_DEBUG_BYTESIZE_FPOSEND(_end)
// #endif

CL_Err_t CL_File::isCalvinFormat(const std::string& pathname)
{
  CL_File dfile;
  dfile.throwOnErr(false);
  return dfile.read_File(pathname,CL_READOPT_HEADONLY);
}

//////////

CL_File::CL_File()
{
  init();
  clear();
}

CL_File::~CL_File()
{
  clear();
}

void CL_File::init()
{
  // refer ourselves to ourself
  setParentFile(this);
  m_bad_datagroup=NULL;
  m_bad_dataset=NULL;
  m_bad_gdh=NULL;
  m_hdr_gdh=NULL;
  //
  m_opt_debug=0;
  m_opt_verbose=0;
}

void CL_File::clear()
{
  delete_ptr_vec(m_datagroups);
  //
  m_filename="";
  m_hdr_magic=-1;
  m_hdr_version=-1;
  //
  delete m_hdr_gdh;
  m_hdr_gdh=NULL;
  delete m_bad_datagroup;
  m_bad_datagroup=NULL;
  delete m_bad_dataset;
  m_bad_dataset=NULL;
  delete m_bad_gdh;
  m_bad_gdh=NULL;
  //
  m_fstream.close();
  m_fstream.clear();
  //
  CL_Object::clear();
}

CL_Err_t CL_File::setErr(CL_Err_t err_num,const std::string& err_msg)
{
  // supply our filename as the last arg to what will be the error message.
  return CL_Object::setErr(err_num,err_msg,m_filename);
}

//////////

std::string CL_File::getFilename() const
{
  return m_filename;
}

//////////

int CL_File::getDataGroupCount() const
{
  return m_datagroups.size();
}

int CL_File::getDataSetCount(int dg_idx) const
{
  if ((dg_idx<0)||(dg_idx>=m_datagroups.size())) {
    return -1;
  }
  return m_datagroups[dg_idx]->getDataSetCount();
}

CL_DataGroup* CL_File::getDataGroup(int dg_idx)
{
  if ((dg_idx<0)||(dg_idx>=m_datagroups.size())) {
    return NULL;
  }
  return m_datagroups[dg_idx];
}

CL_DataGroup* CL_File::getDataGroup(const std::string& dg_name)
{
  for (int dg_idx=0;dg_idx<m_datagroups.size();dg_idx++) {
    if (m_datagroups[dg_idx]->m_name.Nstr()==dg_name) {
      return m_datagroups[dg_idx];
    }
  }
  return NULL;
}

CL_DataSet* CL_File::getDataSet(int dg_idx,int ds_idx)
{
  CL_DataGroup* dg=getDataGroup(dg_idx);
  if (dg==NULL) {
    return NULL;
  }
  return dg->getDataSet(ds_idx);
}

CL_DataSet* CL_File::getDataSet(const std::string& dg_name,const std::string& ds_name)
{
  CL_DataGroup* dg=getDataGroup(dg_name);
  if (dg==NULL) {
    return NULL;
  }
  return dg->getDataSet(ds_name);
}
////////////////////

#define PRINT_STARTSIZE(_var) {                       \
    printf("start=%8d  0x%08x   size=%8d 0x%08x",     \
           _var->getFposStart(),_var->getFposStart(), \
           _var->getFposSize(),_var->getFposSize());  \
  }

#define CL_FPOS_CMP(_lbl,_var,_to) {                                    \
    printf("%s   %d => %d  (= %d)\n",_lbl,_var->getFposStart(),_to,_to-_var->getFposStart()); \
  }

void CL_File::chainDataGroups()
{
  for (int dgi=0;dgi<m_datagroups.size();dgi++) {
    m_datagroups[dgi]->m_next_datagroup=NULL;
    if (dgi>0) {
      m_datagroups[dgi-1]->m_next_datagroup=m_datagroups[dgi];
    }
  }
}

int CL_File::computeSizesAndOffsets()
{
  chainDataGroups();
  //
  int fpos=0;
  //
  setFposStart(fpos);
  int fh_size=bytesize_File(this);

  //
  fpos+=fh_size;
  //fpos=roundToPage(fpos);

  //
  for (int dgi=0;dgi<m_datagroups.size();dgi++) {
    CL_DataGroup* dg=getDataGroup(dgi);
    // chain
    dg->chainDataSets();
    //
    dg->setFposStart(fpos);
    bytesize_DataGroup(dg);
    //
    fpos=dg->getFposEnd();
    fpos=roundToPage(fpos);
    //
    for (int dsi=0;dsi<dg->m_datasets.size();dsi++) {
      CL_DataSet* ds=dg->getDataSet(dsi);
      ds->setFposStart(fpos);
      bytesize_DataSet(ds);
      //
      fpos=ds->getFposEnd();
      fpos=roundToPage(fpos);
    }
  }
  //
  //dumpSegs();
  //
  return CL_OK;
}

////////////////////

bool CL_File::f_open_r()
{
  m_fstream.clear();
  m_fstream.open(m_filename.c_str(),std::ios::in|std::ios::binary);
  bool rv=m_fstream;
  return rv;
}
bool CL_File::f_open_w()
{
  m_fstream.clear();
  m_fstream.open(m_filename.c_str(),std::ios::out|std::ios::binary);
  bool rv=m_fstream;
  return rv;
}

bool CL_File::f_seekg(int pos)
{
  // printf("fseekg: 0x%08x\n",pos); fflush(NULL);
  m_fstream.seekg(pos,std::ios_base::beg);
  return true;
}
bool CL_File::f_seekp(int pos)
{
  // printf("fseekp: 0x%08x\n",pos); fflush(NULL);
  m_fstream.seekp(pos,std::ios_base::beg);
  return true;
}
int CL_File::f_tellg()
{
  return m_fstream.tellp();
}
int CL_File::f_tellp()
{
  return m_fstream.tellp();
}


#define CL_NOTGOOD_RETURN() {                   \
    if (!m_fstream.good()) {                    \
      printf("%s:%d ",__FILE__,__LINE__);       \
      CL_PRINT_POS("notgood return");           \
      setErr(CL_ERR,"file handle not good.");   \
      return errNum();                          \
    }                                           \
  }

CL_Err_t CL_File::read_1B(int& val)
{
  //CL_PRINT_POS("read_1B");
  val=0;
  unsigned char c=0;
  m_fstream.read((char*)&c,1);
  CL_NOTGOOD_RETURN();
  val=c;
  CL_PRINT_DAT_D(val);
  //
  return CL_OK;
}

CL_Err_t CL_File::read_4B(int& val)
{
  //CL_PRINT_POS("read_4B");
  val=0;
  unsigned char buf[4];
  m_fstream.read((char*)buf,4);
  CL_NOTGOOD_RETURN();

  val=(buf[0]<<24)|(buf[1]<<16)|(buf[2]<<8)|(buf[3]);
  //
  CL_PRINT_DAT_D(val);
  return CL_OK;
}

CL_Err_t CL_File::read_4B_Nstring(CL_string& val)
{
  CL_PRINT_POS("read_4B_Nstring");
  int char_len;
  read_4B(char_len);
  return read_Nstring_len(val,char_len);
}

CL_Err_t CL_File::read_4B_Wstring(CL_string& str)
{
  CL_PRINT_POS("read_4B_Wstring");
  int char_len;
  read_4B(char_len);
  // convert chars to bytes
  return read_Wstring_len(str,char_len*2);
}

CL_Err_t CL_File::read_4B_Blob(CL_string& str)
{
  str.clear();
  CL_PRINT_POS("read_4B_Blob");
  int byte_len;
  read_4B(byte_len);
  if (byte_len>0) {
    m_fstream.read(str.resizeWstr(byte_len),byte_len);
  }
  CL_PRINT_DAT_S(str);
  //
  if (!m_fstream) {
    return CL_OK;
  }
  return CL_OK;
}

CL_Err_t CL_File::read_Nstring_len(CL_string& str,int byte_len)
{
  assert(byte_len<10000);
  str.clear();
  if (byte_len>0) {
    m_fstream.read(str.resizeNstr(byte_len),byte_len);
  }
  CL_PRINT_DAT_S(str);
  //
  if (!m_fstream) {
    return CL_ERR;
  }
  return CL_OK;
}

CL_Err_t CL_File::read_Wstring_len(CL_string& str,int byte_len)
{
  assert(byte_len<10000);
  str.clear();
  if (byte_len>0) {
    m_fstream.read(str.resizeWstr(byte_len),byte_len);
  }
  CL_PRINT_DAT_S(str);
  //
  if (!m_fstream) {
    return CL_ERR;
  }
  return CL_OK;
}

////////////////////

CL_File* CL_File::openFile(const std::string& filename)
{
  CL_File* dfile=new CL_File();
  if (dfile->open(filename)!=CL_OK) {
    printf("#%%Error='%s'\n",CL_get_err_string(dfile->m_err_num).c_str());
    delete dfile;
    return NULL;
  }
  return dfile;
}

// should share code with the above...
CL_File* CL_File::openFileNoErr(const std::string& filename)
{
  CL_File* dfile=new CL_File();
  if (dfile->open(filename)!=CL_OK) {
    printf("#%%File='%s'\n",filename.c_str());
    printf("#%%Error='%s'\n",CL_get_err_string(dfile->m_err_num).c_str());
    delete dfile;
    assert(0);
    return NULL;
  }
  return dfile;
}

// Methods are in order of containership.
// File
//   Gdh
//      GDH-leaves
//   Group
//     Dataset
//        DataSetCol
//        DataRow

CL_Err_t CL_File::open(const std::string& filename)
{
  return read_File(filename);
}

CL_Err_t CL_File::close()
{
  return CL_OK;
}

//////////

CL_Gdh* CL_File::getGdh()
{
  if (m_hdr_gdh==NULL) {
    m_hdr_gdh=new CL_Gdh();
  }
  return m_hdr_gdh;
}

CL_Gdh* CL_File::getGdhByNumber(const std::string& path) {
  std::vector<std::string> path_vec;
  Util::chopString(path,':',path_vec);

  return getGdhByPath(path_vec);
}

CL_Gdh* CL_File::getGdhByPath(std::vector<std::string>& path_vec) {
  if (path_vec.size()<1) {
    return NULL;
  }
  if (path_vec[0]!="0") {
    return NULL;
  }
  path_vec.erase(path_vec.begin());
  return getGdh()->getGdhByPath(path_vec);
}

//////////

CL_DataGroup* CL_File::newDataGroup(const std::string& dg_name)
{
  CL_DataGroup* dg=new CL_DataGroup();
  dg->setParentFile(this);
  dg->setName(dg_name);
  m_datagroups.push_back(dg);
  return dg;
}

//////////

CL_Err_t CL_File::read_File(const std::string& filename,CL_ReadOpt_t readopt)
{
  int tmp_dg_start;
  clearErr();
  //
  // printf("===== read_File: '%s' =====\n",filename.c_str());
  //
  clear();
  clearErr();
  //
  m_filename=filename;

  if (!f_open_r()) {
    // printf("open failed.\n");
    setErr(CL_ERR_NOTOPEN,"Unable to open");
    goto EXIT;
  }

  CL_PRINT_POS("read_File");

  setFposStart(f_tellg());

  // gotta have the right magic...
  read_1B(m_hdr_magic);

  // here we try and detect other kinds of celfiles so
  // we can give a better error message.

  // The v3 format starts with
  // [CEL]
  // Version=3
  if (m_hdr_magic==0x5b) { // '['
    setErr(CL_ERR_HDRMAGIC,"This is a v3 celfile, not a calvin celfile.");
    goto EXIT;
  }
  // The v4 format start with 64.
  if (m_hdr_magic==0x40) { // 64
    setErr(CL_ERR_HDRMAGIC,"This is a v4 celfile, not a calvin celfile.");
    goto EXIT;
  }

  // if it is something else just bail.
  if (m_hdr_magic!=0x3b) {
    setErr(CL_ERR_HDRMAGIC,"Bad magic number for calvin (1st byte.)");
    goto EXIT;
  }
  read_1B(m_hdr_version);
  if (m_hdr_version!=0x01) {
    setErr(CL_ERR_HDRMAGIC,"Bad magic number (2nd byte.)");
    goto EXIT;
  }
  // the magic check passed now we treat this as a calvin file.

  // @todo should check for resonable values for m_hdr_datagroup_cnt and m_hdr_datagroup_pos

  //
  read_4B(m_hdr_datagroup_cnt);
  read_4B(m_hdr_datagroup_pos);
  // file GDH
  read_Gdh(getGdh());

  setFposEnd(f_tellg());

  if ((readopt&&CL_READOPT_HEADONLY)==CL_READOPT_HEADONLY) {
    clearErr();
    goto EXIT;
  }

  // groups
  tmp_dg_start=m_hdr_datagroup_pos;
  for (int dgi=0;(dgi<m_hdr_datagroup_cnt)&&(tmp_dg_start!=0);dgi++) {
    CL_DataGroup* dg=new CL_DataGroup();
    dg->setParentFile(this);
    m_datagroups.push_back(dg);
    //
    f_seekg(tmp_dg_start);
    if (read_DataGroup(dg)!=CL_OK) {
      // undo the read.
      m_datagroups.pop_back();
      delete dg;
      //m_bad_datagroup=dg;
      break;
    }
    //
    tmp_dg_start=dg->m_next_datagroup_fpos;
  }

  //
 EXIT:
  CL_PRINT_POS("EOF at:");

  // m_file_end_fpos=f_tellg();
  m_fstream.close();

  //
  //dumpSegs();
  //
  return errNum();
}

/////

CL_Err_t CL_File::read_Gdh(CL_Gdh* gdh)
{
  CL_PRINT_POS("read_Gdh");
  gdh->clear();
  //
  gdh->setFposStart(f_tellg());
  setDataDirty(false);

  CL_NOTGOOD_RETURN();

  //
  read_4B_Nstring(gdh->m_datatype_str);
  read_4B_Nstring(gdh->m_uuid_str);
  //
  read_4B_Wstring(gdh->m_datetime_str);
  read_4B_Wstring(gdh->m_locale_str);

  CL_NOTGOOD_RETURN();

  //
  int param_cnt;
  read_4B(param_cnt);
  gdh->m_params.resize(param_cnt);

  for (int i=0;i<param_cnt;i++) {
    //PRINT_POS("param");
    read_Param(&gdh->m_params[i]);
    //gdh->m_params[i].dump("dbg");
  }

  //
  int tmp_leaf_cnt;
  read_4B(tmp_leaf_cnt);
  gdh->m_gdh_leaves.clear();
  for (int i=0;i<tmp_leaf_cnt;i++) {
    CL_Gdh* leaf_gdh=new CL_Gdh();
    if (read_Gdh(leaf_gdh)==CL_OK) {
      gdh->m_gdh_leaves.push_back(leaf_gdh);
    }
    else {
      m_bad_gdh=leaf_gdh;
      return CL_ERR;
    }
  }
  //
  gdh->setFposEnd(f_tellg());
  //
  //gdh->dump("read_GDH: ");
  return CL_OK;
}

CL_Err_t CL_File::read_Param(CL_Param* param)
{
  CL_NOTGOOD_RETURN();
  //
  param->clear();
  //
  param->setFposStart(f_tellg());
  //
  read_4B_Wstring(param->m_name);
  read_4B_Blob(param->m_val_string);
  CL_string tmp_type_mime;
  read_4B_Wstring(tmp_type_mime);
  param->m_type_code=CL_mimeStrToCode(tmp_type_mime.Nstr());
  //
  param->setFposEnd(f_tellg());

  char* ptr=&param->m_val_string.m_wstr[0];
  // should be a function, but only used here
  switch (param->m_type_code) {
  case CL_TC_BYTE:
    if (param->m_val_string.m_wstr.size()>=1) {
      param->m_val_int=CL_mem_read_i1B(ptr+3);
    }
    break;
  case CL_TC_UBYTE:
    if (param->m_val_string.m_wstr.size()>=1) {
      param->m_val_int=CL_mem_read_u1B(ptr+3);
    }
    break;
  case CL_TC_SHORT:
    if (param->m_val_string.m_wstr.size()>=2) {
      param->m_val_int=CL_mem_read_i2B(ptr+2);
    }
    break;
  case CL_TC_USHORT:
    if (param->m_val_string.m_wstr.size()>=2) {
      param->m_val_int=CL_mem_read_i2B(ptr+2);
    }
    break;
  case CL_TC_INT:
    if (param->m_val_string.m_wstr.size()>=4) {
      param->m_val_int=CL_mem_read_i4B(ptr);
    }
    break;
  case CL_TC_FLOAT:
    if (param->m_val_string.m_wstr.size()>=4) {
      param->m_val_double=CL_mem_read_float(ptr);
    }
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
    // We dont need to do any conversion
    break;
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    // We dont need to do any conversion
    break;
  default:
    // Whoa! we shouldnt have this.
    assert(0);
    break;
  }

  //
  return CL_OK;
}

CL_Err_t CL_File::read_DataGroup(CL_DataGroup* dg)
{
  dg->clear();
  dg->setFposStart(f_tellg());
  CL_PRINT_POS("DATA_GROUP");
  //
  int tmp_dataset_cnt=0;
  read_4B(dg->m_next_datagroup_fpos);
  read_4B(dg->m_first_dataset_fpos);
  read_4B(tmp_dataset_cnt);
  read_4B_Wstring(dg->m_name);
  dg->setFposEnd(f_tellg()+1);

  CL_NOTGOOD_RETURN();

  CL_DataSet* ds;
  int next_dataset_fpos=dg->m_first_dataset_fpos;
  for (int dsi=0;dsi<tmp_dataset_cnt;dsi++) {
    //dumpSegs(); // dbg
    f_seekg(next_dataset_fpos);
    CL_PRINT_POS("seekg");
    CL_NOTGOOD_RETURN();
    //
    ds=new CL_DataSet();
    ds->setParentFile(this);
    ds->setParentGroup(dg);
    dg->m_datasets.push_back(ds);
    //
    //dumpSegs(); // dbg
    //
    if (read_DataSet(ds)!=CL_OK) {
      dg->m_datasets.pop_back();
      delete ds;
      break;
    }
    next_dataset_fpos=ds->m_next_dataset_fpos;
    ds->setFposEnd(next_dataset_fpos);
  }

  //
  dg->chainDataSets();

  //
  return CL_OK;
}

CL_Err_t CL_File::read_DataSet(CL_DataSet* ds)
{
  int tmp_fpos;
  int tmp_data_fpos_start;
  //
  ds->clear();
  //
  ds->setParentFile(this);
  ds->setFposStart(f_tellg());
  CL_PRINT_POS("DATA_SET");
  //
  read_4B(tmp_data_fpos_start); // (1)
  read_4B(ds->m_next_dataset_fpos); // (2)
  read_4B_Wstring(ds->m_name); // (3)

  CL_NOTGOOD_RETURN();

  int tmp_param_cnt;
  read_4B(tmp_param_cnt); // (4)
  ds->m_params.resize(tmp_param_cnt);
  for (int i=0;i<tmp_param_cnt;i++) {
    read_Param(&ds->m_params[i]);
  }

  // read each of the column defs.
  int tmp_col_cnt;
  read_4B(tmp_col_cnt); // (6)
  ds->resizeCols(tmp_col_cnt);
  for (int cidx=0;cidx<tmp_col_cnt;cidx++) {
    read_DataSetCol(ds->getColPtr(cidx));
  }

  // now comes the row count.
  int tmp_row_cnt;
  read_4B(tmp_row_cnt); // (8)
  ds->setRowCount(tmp_row_cnt);

  // now we should be at the start of the data.
  // check this fpos against what we read.
  tmp_fpos=f_tellg();
  assert(tmp_data_fpos_start==tmp_fpos);
  ds->setDataFposStart(tmp_data_fpos_start);

  // read all of the data
  char* buf_ptr=(char*)ds->getBufferPtr();
  int   buf_len=ds->getBufferSize();
  m_fstream.read(buf_ptr,buf_len);
  //CL_PRINT_POS("m_data_end:");

  // did we get the all the data we were expecting?
  int gcount=m_fstream.gcount();
  if (gcount!=buf_len) {
    // no, it is an error
    char err_buf[1024];
    sprintf(err_buf,"short read: gcount=%d  expected_len=%d  (approx %d of %d rows)",
            gcount,buf_len,
            int((gcount/ds->rowBytesize())),ds->rowCount());
    // See APT-612 before fiddling with this.
    // Basiclly, we know that classic calvin will create "short" calvin files,
    // where all the data is not there.  We just have to accept that.
    printf("WARNING: %s\n",err_buf);
    //return setErr(CL_ERR,err_buf);
  }

  //
  int tmp_fpos_end=f_tellg();
  ds->m_bytesize_pad=ds->m_next_dataset_fpos-tmp_fpos_end;
  if (ds->m_bytesize_pad!=0) {
    // printf("padding: %d\n",ds->m_bytesize_pad);
    // dumpSegs();
  }
  // preserve the padding for now.
  tmp_fpos_end=ds->m_next_dataset_fpos;
  //
  ds->setDataFposEnd(tmp_fpos_end);
  ds->setFposEnd(tmp_fpos_end);
  //
  return CL_OK;
}

CL_Err_t CL_File::read_DataSetCol(CL_DataSetCol* dcol)
{
  dcol->clear();
  //
  read_4B_Wstring(dcol->m_name);
  read_1B(dcol->m_type_code);
  read_4B(dcol->m_byte_len);
  //
  return CL_OK;
}

//////////////////////////////


int CL_File::bytesize_Nstr(CL_string& str)
{
  str.ensureNstrValid();
  return 4+str.m_nstr.size();
}

int CL_File::bytesize_Wstr(CL_string& str)
{
  str.ensureWstrValid();
  return 4+str.m_wstr.size();
}

int CL_File::bytesize_ParamVec(std::vector<CL_Param>& params)
{
  int bsize=4; // count + the sizes of each of the params.
  for (int i=0;i<params.size();i++) {
    bsize+=bytesize_Param(&params[i]);
  }
  return bsize;
}

int CL_File::bytesize_Param(CL_Param* param)
{
  CL_string tmp_type_mime;
  tmp_type_mime.setNstr(CL_codeToMimeStr(param->m_type_code));

  int tmp_val_size=0;

  switch (param->m_type_code) {
  case CL_TC_BYTE    :
  case CL_TC_UBYTE   :
  case CL_TC_SHORT   :
  case CL_TC_USHORT  :
  case CL_TC_INT     :
  case CL_TC_UINT    :
  case CL_TC_FLOAT   :
    tmp_val_size=CL_SIZEOF_SCALARBUF; // ick - magic number. (see write_Param)
    break;
  case CL_TC_TEXT_ASCII8  :
  case CL_TC_TEXT_PLAIN8  :
    tmp_val_size=param->m_val_string.Nstr().size();
    break;
  case CL_TC_TEXT_ASCII16 :
  case CL_TC_TEXT_PLAIN16 :
    tmp_val_size=param->m_val_string.Wstr().size();
    break;
  default:
    assert(0);
    break;
  }


  int bsize=
    bytesize_Wstr(param->m_name)+
    4+tmp_val_size+ // 4 is the count, then the actual data size.
    bytesize_Wstr(tmp_type_mime);
  //
  //assert(bsize==param->getFposSize());
  //
  return param->setBytesize(bsize);
}

int CL_File::bytesize_File(CL_File* dfile)
{
  int bsize=1+1+ // magic + version
    4+4+ // dg_cnt + dg_pos
    bytesize_Gdh(getGdh());
  //
  dfile->setFposEnd(dfile->getFposStart()+bsize);
  return dfile->setBytesize(bsize);
}

int CL_File::bytesize_Gdh(CL_Gdh* gdh)
{
  int bsize=0;
  //
  bsize+=bytesize_Nstr(gdh->m_datatype_str);
  bsize+=bytesize_Nstr(gdh->m_uuid_str);
  bsize+=bytesize_Wstr(gdh->m_datetime_str);
  bsize+=bytesize_Wstr(gdh->m_locale_str);
  //
  bsize+=bytesize_ParamVec(gdh->m_params);
  //
  bsize+=4;
  for (int i=0;i<gdh->m_gdh_leaves.size();i++) {
    bsize+=bytesize_Gdh(gdh->m_gdh_leaves[i]);
  }
  //
  //assert(bsize==gdh->getFposSize());
  //
  return gdh->setBytesize(bsize);
}

int CL_File::bytesize_DataGroup(CL_DataGroup* dg)
{
  int bsize=
    4+4+4+ // nextdg fpos + first ds fpos + dataset cnt
    bytesize_Wstr(dg->m_name);
  //
  dg->setFposEnd(dg->getFposStart()+bsize);
  return dg->setBytesize(bsize);
}

// break into to
int CL_File::bytesize_DataSet(CL_DataSet* ds)
{
  // the header
  int bsize=4+ // data start fpos (1)
    4+ // next ds fpos (2)
    bytesize_Wstr(ds->m_name); // name (3)

  //
  bsize+=bytesize_ParamVec(ds->m_params); // (4+5)

  // the column defs
  bsize+=4; // col cnt (6)
  for (int cidx=0;cidx<ds->m_col_vec.size();cidx++) {
    bsize+=bytesize_DataSetCol(&ds->m_col_vec[cidx]);
  }
  bsize+=4; // row count (8)

  // remember
  ds->setDataFposStart(ds->getFposStart()+bsize);

  // be sure the size is up to date
  ds->getBufferPtr();
  bsize+=ds->getBufferSize();

  // padding
  // bsize+=ds->m_fpos_pad_bytesize;

  // the padding does not count in the datasize.
  // (Why isnt just the row-count used?)
  ds->setDataFposEnd(ds->getFposStart()+bsize);
  //
  ds->setFposEnd(ds->getFposStart()+bsize);
  return ds->setBytesize(bsize);
}

int CL_File::bytesize_DataSetCol(CL_DataSetCol* dcol)
{
  int bsize=bytesize_Wstr(dcol->m_name)+
    1+ // type code
    4; // col byte len;
  return bsize;
}

//////////////////////////////

int CL_File::roundTo1KB(int val)
{
  return (val+0x3FF)&~0x3FF;
}

int CL_File::roundToPage(int val)
{
  //int bits=10;
  //int mask=(1<<bits)-1;
  //val=(val+mask)&~mask;
  return val;
}

CL_Err_t CL_File::write_Zeros(int cnt) {
  return write_Padding(cnt,0);
}

CL_Err_t CL_File::write_Padding(int cnt,int val)
{
  // get a block of zeros ready...
  char buf[1024];
  memset(buf,val,sizeof(buf));

  // ...write full blocks...
  while (cnt>=sizeof(buf)) {
    m_fstream.write(buf,sizeof(buf));
    cnt-=sizeof(buf);
  }

  // ... and then write the last partial block.
  if (cnt>0) {
    m_fstream.write(buf,cnt);
  }
  return CL_OK;
}

CL_Err_t CL_File::write_1B(int val)
{
  char c;
  c=val;
  m_fstream.put(c);
  //
  return CL_OK;
}

CL_Err_t CL_File::write_4B(int val)
{
  m_fstream.put((val>>24)&0xFF);
  m_fstream.put((val>>16)&0xFF);
  m_fstream.put((val>> 8)&0xFF);
  m_fstream.put((val    )&0xFF);
  //
  return CL_OK;
}

CL_Err_t CL_File::write_4B_Nstring(const std::string& val)
{
  int char_len=val.size();
  write_4B(char_len);
  m_fstream.write(val.c_str(),char_len);
  //
  return CL_OK;
}

CL_Err_t CL_File::write_4B_Nstring(CL_string& val)
{
  val.ensureNstrValid();
  int char_len=val.m_nstr.size();
  write_4B(char_len);
  m_fstream.write(val.m_nstr.c_str(),char_len);
  //
  return CL_OK;
}

CL_Err_t CL_File::write_4B_Wstring(const std::string& val)
{
  // convert and punt
  CL_string tmp_val;
  tmp_val.setNstr(val);
  return write_4B_Wstring(tmp_val);
}

CL_Err_t CL_File::write_4B_Wstring(CL_string& val)
{
  val.ensureWstrValid();
  int byte_len=val.m_wstr.size();
  // Note: convert the byte_len to char_len
  write_4B(byte_len/2);
  m_fstream.write(val.m_wstr.c_str(),byte_len);
  //
  return CL_OK;
}

CL_Err_t CL_File::write_4B_Blob(CL_string& val)
{
  val.ensureWstrValid();
  int byte_len=val.m_wstr.size();
  write_4B(byte_len); // Note!
  m_fstream.write(val.m_wstr.c_str(),byte_len);
  //
  return CL_OK;
}

//////////

CL_Err_t CL_File::write(const std::string& filename)
{
  return write_File(filename);
}

CL_Err_t CL_File::write_File(const std::string& filename)
{
  return write_File(filename,true);
}

CL_Err_t CL_File::write_File(const std::string& filename,bool forceWrite)
{
  computeSizesAndOffsets();
  //dumpSegs();

  //
  m_filename=filename;
  if (!f_open_w()) {
    // printf("open failed!\n");
    return CL_ERR_NOTOPEN;
  }

  CL_DEBUG_BYTESIZE_START("file:1:");

  write_FileHdr(forceWrite);

  CL_DEBUG_BYTESIZE("file:1:",getBytesize());

  for (int dgi=0;dgi<getDataGroupCount();dgi++) {
    CL_DataGroup* dg=getDataGroup(dgi);
    write_DataGroup(dg,forceWrite);
    for (int dsi=0;dsi<dg->getDataSetCount();dsi++) {
      CL_DataSet* ds=dg->getDataSet(dsi);
      write_DataSet(ds,forceWrite);
    }
  }

  //
  m_fstream.close();
  return CL_OK;
}

CL_Err_t CL_File::write_FileHdr(bool forceWrite)
{
  if ((isDataDirty()==false) && (forceWrite==false)) {
    return CL_OK;
  }
  //
  f_seekp(getFposStart());
  //
  write_1B(0x3b);
  write_1B(0x01);

  int tmp_datagroups_size=m_datagroups.size();
  write_4B(tmp_datagroups_size);
  if (tmp_datagroups_size>0) {
    write_4B(m_datagroups[0]->getFposStart());
  }
  else {
    write_4B(0);
  }

  write_Gdh(getGdh(),forceWrite);
  //
  setDataDirty(false);
  return CL_OK;
}

CL_Err_t CL_File::write_Gdh(CL_Gdh* gdh,bool forceWrite)
{
  if ((gdh->isDataDirty()==false) && (forceWrite==false)) {
    return CL_OK;
  }

  //gdh->dump("write:");
  CL_DEBUG_BYTESIZE_START("GDH");
  //
  write_4B_Nstring(gdh->m_datatype_str);
  write_4B_Nstring(gdh->m_uuid_str);
  //
  write_4B_Wstring(gdh->m_datetime_str);
  write_4B_Wstring(gdh->m_locale_str);

  int tmp_param_cnt=gdh->m_params.size();
  write_4B(tmp_param_cnt);
  for (int pidx=0;pidx<tmp_param_cnt;pidx++) {
    write_Param(&gdh->m_params[pidx]);
  }

  int tmp_gdh_leaf_cnt=gdh->m_gdh_leaves.size();
  write_4B(tmp_gdh_leaf_cnt);
  for (int lidx=0;lidx<tmp_gdh_leaf_cnt;lidx++) {
    write_Gdh(gdh->m_gdh_leaves[lidx],forceWrite);
  }

  CL_DEBUG_BYTESIZE("write_GDH:",gdh->m_bytesize);
  //CL_DEBUG_BYTESIZE_FPOSEND("write_GDH:",gdh->getFposEnd());

  return CL_OK;
}

CL_Err_t CL_File::write_Param(CL_Param* param)
{
  CL_DEBUG_BYTESIZE_START("Param");
  //
  write_4B_Wstring(param->m_name);
  // Convert the value and
  switch (param->m_type_code) {
  case CL_TC_BYTE:
  case CL_TC_UBYTE:
    param->m_val_string.resizeBlob(CL_SIZEOF_SCALARBUF);
    CL_mem_write_i1B((void*)param->m_val_string.blob_str(),param->m_val_int);
    break;
  case CL_TC_SHORT:
  case CL_TC_USHORT:
    param->m_val_string.resizeBlob(CL_SIZEOF_SCALARBUF);
    CL_mem_write_i2B((void*)param->m_val_string.blob_str(),param->m_val_int);
    break;
  case CL_TC_INT:
  case CL_TC_UINT:
    param->m_val_string.resizeBlob(CL_SIZEOF_SCALARBUF);
    CL_mem_write_i4B((void*)param->m_val_string.blob_str(),param->m_val_int);
    break;
  case CL_TC_FLOAT:
    param->m_val_string.resizeBlob(CL_SIZEOF_SCALARBUF);
    CL_mem_write_float((void*)param->m_val_string.blob_str(),param->m_val_double);
    break;
  case CL_TC_TEXT_ASCII8  :
  case CL_TC_TEXT_ASCII16 :
  case CL_TC_TEXT_PLAIN8  :
  case CL_TC_TEXT_PLAIN16 :
    // dont need to do anything; The string has the value.
    break;
  default:
    assert(0);
    break;
  }
  // now write the string.
  write_4B_Blob(param->m_val_string);
  //
  std::string tmp_type_mime=CL_codeToMimeStr(param->m_type_code);
  write_4B_Wstring(tmp_type_mime);

  CL_DEBUG_BYTESIZE("write_Param:",param->getBytesize());
  //CL_DEBUG_BYTESIZE_FPOSEND("write_Param:",param->getFposEnd());
  //
  return CL_OK;
}

CL_Err_t CL_File::write_DataGroup(CL_DataGroup* dg,bool forceWrite)
{
  if ((dg->isDataDirty()==false) && (forceWrite==false)) {
    return CL_OK;
  }

  //
  dg->chainDataSets();
  //
  int tmp_datasets_cnt=dg->m_datasets.size();
  //
  f_seekp(dg->getFposStart());
  //
  write_4B(dg->nextDataGroupFposStart());
  write_4B(dg->firstDataSetFposStart());
  write_4B(tmp_datasets_cnt);
  write_4B_Wstring(dg->m_name);

  // this is being done in the file...
  //  for (int dsi=0;dsi<tmp_datasets_cnt;dsi++) {
  //write_DataSet(dg->m_datasets[dsi]);
  //}
  setDataDirty(false);
  return CL_OK;
}

CL_Err_t CL_File::write_DataSet(CL_DataSet* ds,bool forceWrite)
{
  if ((ds->isDataDirty()==false) && (forceWrite==false)) {
    return CL_OK;
  }
  //
  f_seekp(ds->getFposStart());
  //
  int tmp_data_fpos_start=ds->getDataFposStart();
  write_4B(tmp_data_fpos_start); // (1)
  //
  int next_dset_fpos=0;
  if (ds->m_dataset_next==NULL) {
    next_dset_fpos=ds->getFposEnd();
  }
  else {
    next_dset_fpos=ds->m_dataset_next->getFposStart();
  }
  // dumpSegs();
  // printf("next_dset_fpos=0x%08x\n",next_dset_fpos);
  write_4B(next_dset_fpos); // (2)
  //
  write_4B_Wstring(ds->m_name); // (3)

  //
  int tmp_param_cnt=ds->m_params.size();
  write_4B(tmp_param_cnt); // (4)
  for (int i=0;i<tmp_param_cnt;i++) {
    write_Param(&ds->m_params[i]);
  }

  //
  write_4B(ds->m_col_vec.size()); // (6)
  for (int cidx=0;cidx<ds->m_col_vec.size();cidx++) {
    write_DataSetColDef(&ds->m_col_vec[cidx]);
  }
  write_4B(ds->rowCount()); // (8)

  //
  int tmp_fpos=f_tellp();
  assert(tmp_fpos==tmp_data_fpos_start);

  // write the data
  char* buf_ptr=(char*)ds->getBufferPtr();
  int   buf_len=ds->getBufferSize();
  // only do the write if we have data to write.
  // (on unix a zero-len write wont dereference buf_ptr; Windows will)
  if (buf_len>0) {
    assert(buf_ptr!=NULL);
    m_fstream.write(buf_ptr,buf_len);
  }

  //
  write_Zeros(ds->m_bytesize_pad);

  //
  setDataDirty(false);
  return CL_OK;
}

CL_Err_t CL_File::write_DataSetColDef(CL_DataSetCol* dcol)
{
  write_4B_Wstring(dcol->m_name);
  write_1B(dcol->m_type_code);
  write_4B(dcol->m_byte_len);
  return CL_OK;
}

//////////

CL_Err_t CL_File::dump()
{
  std::string prefix;

  CL_DUMP_PREFIX("F   ");
  CL_DUMP_MEMB_INT(m_hdr_magic);
  CL_DUMP_MEMB_INT(m_hdr_version);
  CL_DUMP_MEMB_INT(m_hdr_datagroup_cnt);
  CL_DUMP_MEMB_INT(m_hdr_datagroup_pos);
  CL_DUMP_MEMB_INT(m_err_num);
  CL_DUMP_MEMB_INT(m_file_end_fpos);

  getGdh()->dump(prefix);

  for (int gi=0;gi<m_datagroups.size();gi++) {
    m_datagroups[gi]->dump(prefix);
  }
  return CL_OK;
}

CL_Err_t CL_File::dumpSegs()
{
  printf("F      |      | %p | %08x %08x %08x | '%s'\n",
         this,
         getFposStart(),
         getFposSize(),
         getFposEnd(),
         m_filename.c_str());
  //
  for (int dgi=0;dgi<m_datagroups.size();dgi++) {
    CL_DataGroup* dg=getDataGroup(dgi);
    printf("  DG   | %2d   | %p | %08x %08x %08x |   '%s'\n",
           dgi,dg,
           dg->getFposStart(),
           dg->getFposSize(),
           dg->getFposEnd(),
           dg->m_name.c_str());
    //
    for (int dsi=0;dsi<dg->m_datasets.size();dsi++) {
      CL_DataSet* ds=dg->getDataSet(dsi);
      printf("    DS |   %2d | %p | %08x %08x %08x |     '%s' (%dc x %dr)\n",
             dsi, ds,
             ds->getFposStart(),
             ds->getFposSize(),
             ds->getFposEnd(),
             ds->m_name.c_str(),
             ds->colCount(),
             ds->rowCount());
    }
  }
  return CL_OK;
}

//////////

bool CL_File::requiresRewrite()
{
  computeSizesAndOffsets();

  if (changedFpos()) {
    return true;
  }
  for (int dgi=0;dgi<m_datagroups.size();dgi++) {
    CL_DataGroup* dg=getDataGroup(dgi);
    if (dg->changedFpos()) {
      return true;
    }
    for (int dsi=0;dsi<dg->m_datasets.size();dsi++) {
      CL_DataSet* ds=dg->getDataSet(dsi);
      if (ds->changedFpos()) {
        return true;
      }
    }
  }
  //
  return false;
}

//////////

CL_Err_t CL_File::copyFile(const std::string& from, const std::string& to)
{
  CL_Err_t rv=CL_ERR;
  CL_File* dfile=new CL_File();
  if (dfile->read_File(from)==CL_OK) {
    if (dfile->write_File(to)==CL_OK) {
      rv=CL_OK;
    }
  }
  //
  dfile->close();
  delete dfile;
  return rv;
}

//////////

CL_Err_t CL_File::dumpFile(const std::string& filename)
{
  CL_File* dfile=new CL_File();
  if (dfile->read_File(filename)!=CL_OK) {
    printf("#%%Error='%s'\n",CL_get_err_string(dfile->m_err_num).c_str());
    delete dfile;
    return CL_ERR;
  }

  dfile->dump();

  dfile->close();
  delete dfile;
  return CL_OK;
}

CL_Err_t CL_File::dumpFileSegs(const std::string& filename)
{
  CL_File* dfile=new CL_File();
  if (dfile->open(filename)!=CL_OK) {
    printf("#%%Error='%s'\n",CL_get_err_string(dfile->m_err_num).c_str());
    delete dfile;
    return CL_ERR;
  }

  dfile->dumpSegs();

  dfile->close();
  delete dfile;
  return CL_OK;
}

//////////

CL_Err_t CL_File::changeFile(const std::string& filename,
                             const CL_changevec_t& change_vec,
                             bool substr)
{
  int cnt;
  CL_File* dfile=new CL_File();
  if (dfile->read_File(filename)!=CL_OK) {
    printf("#%%Error='%s'\n",CL_get_err_string(dfile->m_err_num).c_str());
    delete dfile;
    return CL_ERR;
  }

  int total_changes=0;
  for (int i=0;i<change_vec.size();i++) {
    cnt=0;
    std::string from=change_vec[i].first;
    std::string to  =change_vec[i].second;
    printf("# applying %schange '%s' to '%s'...\n",
           ((substr)?"substr ":""),
           from.c_str(),to.c_str());
    changeFileGdh(dfile->getGdh(),from,to,substr,&cnt);
    printf("# ...made %d changes\n",cnt);
    total_changes+=cnt;
  }
  printf("# total changes: %d\n",total_changes);

  dfile->write(filename+".changed");

  delete dfile;
  return CL_OK;
}

CL_Err_t CL_File::changeFileGdh(CL_Gdh* gdh,
                                const std::string& from,
                                const std::string& to,
                                bool substr,
                                int* cnt)
{
  // do the four specials
  changeFileGdhValue("DATATYPE ",&gdh->m_datatype_str ,from,to,substr,cnt);
  changeFileGdhValue("UUID     ",&gdh->m_uuid_str     ,from,to,substr,cnt);
  changeFileGdhValue("DATETIME ",&gdh->m_datetime_str ,from,to,substr,cnt);
  changeFileGdhValue("LOCALE   ",&gdh->m_locale_str   ,from,to,substr,cnt);
  
  // do our params...
  int param_cnt=gdh->m_params.size();
  for (int pidx=0;pidx<param_cnt;pidx++) {
    changeFileGdhValue(gdh->m_params[pidx].m_name.Nstr(),
                       &gdh->m_params[pidx].m_val_string,
                       from,to,
                       substr,cnt);
  }

  // ...then recurse down the tree.
  int leaf_cnt=gdh->m_gdh_leaves.size();
  for (int lidx=0;lidx<leaf_cnt;lidx++) {
    changeFileGdh(gdh->m_gdh_leaves[lidx],from,to,substr,cnt);
  }
  return CL_OK;
}

CL_Err_t CL_File::changeFileGdhValue(const std::string& key,
                                     CL_string* clstr,
                                     const std::string& from,
                                     const std::string& to,
                                     bool substr,
                                     int* cnt)
{
  bool changed=false;
  std::string tmp_str;
  std::string orig_str;

  // this truncs the 'nulls'
  orig_str=clstr->Nstr();
  tmp_str=orig_str;

  // some debugging
  // printf("# %d: '%s' => '%s'/%d\n",idx,param->m_name.c_str(),tmp_str.c_str(),(int)tmp_str.size());

  if (substr) {
    std::string::size_type pos=tmp_str.find(from);
    if (pos!=std::string::npos) {
      tmp_str.replace(pos,from.size(),to);
      clstr->setNstr(tmp_str);
      changed=true;
    }
  }
  else {
    if (tmp_str==from) {
      clstr->setNstr(to);
      changed=true;
    }
  }
  if (changed) {
    printf("### changed: '%s'='%s' to '%s'\n",
           key.c_str(),
           orig_str.c_str(),
           clstr->c_str());
    (*cnt)++;
  }

  return CL_OK;
}

//////////

CL_Err_t CL_File::printFileParamsAll(const std::string& filename)
{
  CL_File* dfile=CL_File::openFile(filename);
  CL_Err_t rv=printFileParamsAll(dfile);
  delete dfile;
  return rv;
}

CL_Err_t CL_File::printFileParamsAll(const CL_File* dfile)
{
  if (dfile==NULL) {
    return CL_ERR;
  }
  dfile->m_hdr_gdh->print("");
  return CL_OK;
}

CL_Err_t CL_File::printFileParams(const std::string& filename,
                                  const std::vector<std::string>& key_vec)
{
  CL_File* dfile=CL_File::openFile(filename);
  CL_Err_t rv=printFileParams(dfile,key_vec);
  delete dfile;
  return rv;
}

CL_Err_t CL_File::printFileParams(const CL_File* dfile,
                                  const std::vector<std::string>& key_vec)
{
  if (dfile==NULL) {
    return CL_ERR;
  }
  dfile->m_hdr_gdh->printMatchingKeys(key_vec,"");
  return CL_OK;
}

//////////

#if CL_WITH_TSVFILE==1

CL_Err_t CL_File::ExportToTsv(const std::string& filename,
                              int dg_idx,int ds_idx,
                              const std::string& out_dir)
{
  CL_Err_t rv;
  CL_File* dfile=new CL_File();

  if ((rv=dfile->open(filename))!=CL_OK) {
    dfile->printErrMsg();
    delete dfile;
    return rv;
  }
  std::string out_base=out_dir+"/"+filename;
  rv=dfile->tsvExport(dg_idx,ds_idx,out_base);
  //
  delete dfile;
  return rv;
}

CL_Err_t CL_File::tsvExport(int dg_idx,int ds_idx,
                            const std::string& out_base)
{
  CL_Err_t rv=CL_OK;
  char buf[2000];

  for (int dgi=0;dgi<getDataGroupCount();dgi++) {
    if ((dg_idx==-1)||(dgi==dg_idx)) {
      for (int dsi=0;dsi<getDataSetCount(dgi);dsi++) {
        if ((ds_idx==-1)||(dsi==ds_idx)) {
          CL_DataSet* ds=getDataSet(dgi,dsi);
          assert(ds!=NULL);
          //
          sprintf(buf,"%s-g%02d-s%02d.tsv",out_base.c_str(),dgi,dsi);
          std::string tsv_name=buf;
          m_tsvexport_filenames.push_back(tsv_name);
          rv=tsvExportDataSet(ds,dgi,dsi,tsv_name);
          if (rv!=CL_OK) {
            break;
          }
        }
      }
    }
  }
  // If all, then do the GDH as well...
  if ((dg_idx==-1) && (ds_idx==-1)) {
    tsvExportGdh(m_hdr_gdh,out_base+"-gdh.tsv");
  }
  
  return rv;
}

//////////

CL_Err_t CL_File::ExportGdhToTsv(const std::string& filename,
                                 const std::string& outdir)
{
  CL_Err_t rv;
  // @todo outdir is unused
  CL_File* dfile=new CL_File();
  if ((rv=dfile->open(filename))!=CL_OK) {
    /// @todo throw an exception
    printf("#%%Error=file '%s' not found.\n",filename.c_str());
    delete dfile;
    return rv;
  }

  rv=dfile->tsvExportGdh(dfile->m_hdr_gdh,filename+"-gdh.tsv");

  delete dfile;
  return rv;
}

CL_Err_t CL_File::tsvExportGdh(CL_Gdh* gdh,
                               const std::string& tsv_name)
{
  printf("Exporting GDH to '%s'...\n",tsv_name.c_str());

  // create the TsvFile
  affx::TsvFile tsv;
  //
  std::string cmd;

  tsv.addHeaderComment("!/usr/bin/env apt-calvinlite-util");
  tsv.addHeader("import-out",m_filename);
  m_tsvexport_filenames.push_back(tsv_name);
  tsv.addHeader("import-file",m_tsvexport_filenames);
  //
  tsv.addHeaderComment("We cant say this due to '#!' limits:");
  cmd="!/usr/bin/env apt-calvinlite-util --tsv-import "+m_filename+"";
  cmd+=" "+tsv_name+" ";
  for (int i=0;i<m_tsvexport_filenames.size();i++) {
    cmd+=" "+m_tsvexport_filenames[i];
  }
  tsv.addHeaderComment(cmd);
  //
  tsv.addHeaderComment("");
  tsv.addHeaderComment(" These are the Generic Data Headers (GDH) for a calvin file.");
  tsv.addHeaderComment(" The 'gdh_number' column is the colon-seperated indexes of the GDH.");
  tsv.addHeaderComment(" '0' is the root GDH; '0:0' is the first child of the root; usw...");
  tsv.addHeaderComment(" The raw_bytes column are the raw bytes. (only used for debugging)");
  tsv.addHeaderComment("");
  //
  tsv.addHeader("calvin-file-name",getFilename());
  tsv.addHeader("calvin-part","gdh");
  //
  tsv.addHeaderComment("");
  //
  tsv.defineColumn(0,0,"gdh_number");
  //
  tsv.defineColumn(1,0,"name");
  tsv.defineColumn(1,1,"type");
  tsv.defineColumn(1,2,"byte_len");
  tsv.defineColumn(1,3,"value");
  //tsv.defineColumn(1,4,"raw_bytes");
  // v2 is required -- multilevel format.
  tsv.writeTsv_v2(tsv_name);

  //
  tsvExportGdh(gdh,"0",&tsv);

  //
  tsv.close();

  return CL_OK;
}

CL_Err_t CL_File::tsvExportGdh(CL_Gdh* gdh,
                               const std::string& gdh_name,
                               affx::TsvFile* tsv)
{
  // What is this GDH called...
  tsv->set(0,0,gdh_name);
  tsv->writeLevel(0);

  // These are special, they dont have types encoded like the params,
  // but we write them out as if they were.
  // the reader can special case them.
  tsvExportGdh_Param(tsv,"datatype" ,gdh->m_datatype_str.Nstr());
  tsvExportGdh_Param(tsv,"uuid"     ,gdh->m_uuid_str.Nstr());
  tsvExportGdh_Param(tsv,"datetime" ,gdh->m_datetime_str.Nstr());
  tsvExportGdh_Param(tsv,"locale"   ,gdh->m_locale_str.Nstr());

  //
  for (int i=0;i<gdh->getParamCount();i++) {
    CL_Param* param=gdh->getParam(i);
    tsvExportGdh_Param(tsv,param);
  }

  // to leaf-gdhs.
  for (int i=0;i<gdh->getLeafCount();i++) {
    CL_Gdh* leaf_gdh=gdh->getLeafGdh(i);
    tsvExportGdh(leaf_gdh,gdh_name+":"+ToStr(i),tsv);
  }
  //
  return CL_OK;
}

CL_Err_t CL_File::tsvExportGdh_Param(affx::TsvFile* tsv,
                                     const std::string& key_str,
                                     const std::string& val_str)
{
  return tsvExportGdh_Param(tsv,
                            key_str,
                            CL_TC_TEXT_PLAIN8,val_str.size(),
                            val_str,"");
}

CL_Err_t CL_File::tsvExportGdh_Param(affx::TsvFile* tsv,
                                     CL_Param* param)
{
  std::string key_str=param->m_name.Nstr();
  int byte_len=param->suggestedByteLen();
  CL_TypeCode_t val_type_code=param->m_type_code;
  std::string val_str=param->valueAsString();
  std::string val_raw=param->m_val_string.Nstr();

  return tsvExportGdh_Param(tsv,
                            key_str,
                            val_type_code,byte_len,
                            val_str,val_raw);
}

CL_Err_t CL_File::tsvExportGdh_Param(affx::TsvFile* tsv,
                                     const std::string& key_str,
                                     CL_TypeCode_t val_type_code,
                                     int byte_len,
                                     const std::string& val_str,
                                     const std::string& val_raw)
{
  tsv->set(1,0,key_str);
  std::string val_type_mime=CL_codeToMimeStr(val_type_code);
  tsv->set(1,1,val_type_mime);
  tsv->set(1,2,byte_len);
  tsv->set(1,3,val_str);
  // for debugging
  //tsv->set(1,4,"'"+val_raw+"'");
  //
  tsv->writeLevel(1);
  //
  return CL_OK;
}

//////////

CL_Err_t CL_File::tsvExportDataSet(CL_DataSet* ds,int dg_idx,int ds_idx,const std::string& tsv_name)
{
  char buf[2000];
  affx::TsvFile tsv;

  printf("Exporting in tsv format to '%s'\n",tsv_name.c_str());
  fflush(NULL);

  tsv.addHeaderComment(" This is a dump of a calvin file.");
  tsv.addHeaderComment(" The headers starting with 'calvin-' describe the size and shape of the dataset.");
  tsv.addHeaderComment(" Changing this data will alter the Calvin DataSet.");
  tsv.addHeaderComment("");

  //
  CL_File* dfile=ds->getParentFile();
  if (dfile!=NULL) {
    tsv.addHeader("calvin-file-name",dfile->getFilename());
  }
  tsv.addHeader("calvin-part","dataset");
  tsv.addHeaderComment("");
  //
  CL_DataGroup* dgroup=ds->getParentGroup();
  if (dgroup!=NULL) {
    tsv.addHeader("calvin-datagroup-name",dgroup->getName());
  }
  tsv.addHeader("calvin-datagroup-idx",dg_idx);
  tsv.addHeaderComment("");
  //
  tsv.addHeader("calvin-dataset-name",ds->getName());
  tsv.addHeader("calvin-dataset-idx",ds_idx);
  // quick ref for the user.
  tsv.addHeader("calvin-dataset-rowcount",ds->rowCount());
  tsv.addHeaderComment("");

  //
  int col_cnt=ds->colCount(); // x
  int row_cnt=ds->rowCount(); // y

  //
  std::string cname;
  for (int cidx=0;cidx<col_cnt;cidx++) {
    ds->getColName(cidx,&cname);
    tsv.defineColumn(0,cidx,cname);
    //
    CL_DataSetCol* dcol=ds->getColPtr(cidx);
    snprintf(buf,sizeof(buf),"calvin-col-%d-datatype",cidx);
    tsv.addHeader(buf,dcol->m_type_code);
    snprintf(buf,sizeof(buf),"calvin-col-%d-bytelen",cidx);
    tsv.addHeader(buf,dcol->m_byte_len);
  }
  tsv.addHeaderComment("");

  //
  if (tsv.writeTsv_v1(tsv_name)!=affx::TSV_OK) {
    return CL_ERR;
  }

  //
  std::string val;
  for (int ridx=0;ridx<row_cnt;ridx++) {
    for (int cidx=0;cidx<col_cnt;cidx++) {
      ds->getAsString(ridx,cidx,&val);
      tsv.set(0,cidx,val);
    }
    tsv.writeLevel(0);
  }

  //
  tsv.close();
  //
  return CL_OK;
}

/////

CL_Err_t CL_File::tsvImport(const std::string& file_in,
                            const std::string& file_out,
                            const std::vector<std::string>& tsvfiles)
{
  CL_File* dfile=new CL_File();

  CL_Err_t rv;

  if (file_in!="") {
    rv=isCalvinFormat(file_in);
  }
  else {
    rv=CL_ERR_NOTOPEN;
  }

  if (rv==CL_OK) {
    rv=dfile->read_File(file_in);
    if (rv!=CL_OK) {
      printf("#%%Error=file '%s' not found.\n",file_in.c_str());
      delete dfile;
      return rv;
    }
  }
  else if (rv==CL_ERR_NOTOPEN) {
    // we will make it.
  }
  else {
    delete dfile;
    return rv;
  }

  // now process all the parts to import.
  for (int i=0;i<tsvfiles.size();i++) {
    dfile->tsvImport(tsvfiles[i]);
  }

  dfile->write(file_out);

  delete dfile;
  return CL_OK;
}

CL_Err_t CL_File::tsvImportScript(const std::string& script_name)
{
  affx::TsvFile script_tsv;
  std::string tmp_import_out;
  std::vector<std::string> tmp_imports;

  printf("#%% importing from script: '%s'\n",script_name.c_str());

  if (script_tsv.open(script_name)!=affx::TSV_OK) {
    printf("#%% bad tsv!\n");
    return CL_ERR;
  }
  if (script_tsv.getHeader("import-out",tmp_import_out)!=affx::TSV_OK) {
    printf("#%% bad header import-out!\n");
    return CL_ERR;
  }

  if (script_tsv.getHeader("import-file",tmp_imports)!=affx::TSV_OK) {
    printf("#%% bad header import-file!\n");
    return CL_ERR;
  }
    
  return CL_File::tsvImport("",tmp_import_out,tmp_imports);
}

//////////

CL_Err_t CL_File::tsvImport(const std::string& tsv_name)
{
  std::string calvin_part;

  affx::TsvFile* tsv=new affx::TsvFile();
  // APT-836: the export files sometimes have binary data in the comment.
  // just try and process them.
  tsv->m_optCheckFormatOnOpen=false;
  tsv->open(tsv_name);
  int rv=tsv->getHeader("calvin-part",calvin_part);
  delete tsv;

  if (rv!=affx::TSV_OK) {
    printf("#%% error: no calvin-part header.\n");
    return setErr(CL_ERR,"no calvin-part header.");
  }

  if (calvin_part=="dataset") {
    return tsvImportData(tsv_name);
  }
  else if (calvin_part=="gdh") {
    return tsvImportGdh(tsv_name);
  }
  // else
  printf("#%% calvin-part should be 'dataset' or 'gdh', not '%s'",calvin_part.c_str());
  return setErr(CL_ERR,"bad calvin-part header");
}

CL_Err_t CL_File::tsvImportData(const std::string& tsv_name)
{
  affx::TsvFile tsv;

  printf("#%% Importing data from '%s'...\n",tsv_name.c_str());

  if (tsv.open(tsv_name)!=affx::TSV_OK) {
    printf("#%% error: tsv open failed.\n");
    return CL_ERR;
  }

  // where does dataset go?
  std::string datagroup_name;
  //int datagroup_idx;
  std::string dataset_name;
  //int dataset_idx;

  if ((tsv.getHeader("calvin-datagroup-name",datagroup_name)!=affx::TSV_OK) ||
      (tsv.getHeader("calvin-dataset-name",dataset_name)!=affx::TSV_OK))
    {
      printf("#%% error: missing required headers: calvin-datagroup-name and calvin-dataset-name.\n");
      return CL_ERR;
    }

  if ((datagroup_name=="")||(dataset_name=="")) {
    printf("#%% '%s' needs both a datagroup name and a dataset name.\n",tsv_name.c_str());
    printf("#%% Got '%s' and '%s'.\n",datagroup_name.c_str(),dataset_name.c_str());
    return CL_ERR;
  }

  /// @todo getDataGroupCreate(name,name);
  CL_DataGroup* dgroup=NULL;
  CL_DataSet* dset=NULL;

  dgroup=getDataGroup(datagroup_name);
  if (dgroup==NULL) {
    dgroup=newDataGroup(datagroup_name);
  }
  dset=dgroup->getDataSet(dataset_name);
  if (dset==NULL) {
    dset=dgroup->newDataSet(dataset_name);
  }

  tsvImportDataCopy(dset,&tsv);

  tsv.close();
  return CL_OK;
}

CL_Err_t CL_File::tsvImportDataCopy(CL_DataSet* dset,affx::TsvFile* tsv)
{
  char buf[100];

  dset->clear();
  const int cidx_max=tsv->getColumnCount(0);

  printf("#%% copying %d columns of data from '%s' to '%s'...\n",
         cidx_max,tsv->m_fileName.c_str(),dset->m_name.c_str());

  std::string col_name;
  int col_dtype_int;
  int col_dlen;
  std::vector<int> col_dtypes;

  for (int cidx=0;cidx<cidx_max;cidx++) {
    col_name=tsv->getColumnName(0,cidx);
    sprintf(buf,"calvin-col-%d-datatype",cidx);
    tsv->getHeader(buf,col_dtype_int);
    sprintf(buf,"calvin-col-%d-bytelen",cidx);
    tsv->getHeader(buf,col_dlen);
    //
    CL_TypeCode_t col_dtype=(CL_TypeCode_t)col_dtype_int;
    dset->newColumn(col_name,col_dtype,col_dlen);
    //
    col_dtypes.push_back(col_dtype);
  }

  //
  int row_cnt=0;
  tsv->rewind();
  while (tsv->nextLevel(0)==affx::TSV_OK) {
    row_cnt++;
  }
  dset->setRowCount(row_cnt);
  tsv->rewind();

  //
  std::string val_string;
  int         val_int;
  double      val_double;

  //
  int ridx=0;
  while (tsv->nextLevel(0)==affx::TSV_OK) {
    for (int cidx=0;cidx<cidx_max;cidx++) {
      switch (col_dtypes[cidx]) {
      case CL_TC_BYTE    :
      case CL_TC_UBYTE   :
      case CL_TC_SHORT   :
      case CL_TC_USHORT  :
      case CL_TC_INT     :
      case CL_TC_UINT    :
        tsv->get(0,cidx,val_int);
        dset->set(ridx,cidx,val_int);
        break;
      case CL_TC_FLOAT    :
        tsv->get(0,cidx,val_double);
        dset->setFloat(ridx,cidx,val_double);
        break;
      case CL_TC_TEXT_ASCII8  :
      case CL_TC_TEXT_ASCII16 :
      case CL_TC_TEXT_PLAIN8  :
      case CL_TC_TEXT_PLAIN16 :
        tsv->get(0,cidx,val_string);
        dset->set(ridx,cidx,val_string);
        break;
      default:
        assert(0);
        return CL_ERR;
      }
    }
    ridx++;
  }

  //
  printf("#%% imported %d rows.\n",ridx);
  return CL_OK;
}

//////////

CL_Err_t CL_File::tsvImportGdh(const std::string& tsv_gdh_name) {
  affx::TsvFile tsv_gdh;
  // dont abort when seeing binary data in "raw_bytes" col
  tsv_gdh.m_optCheckFormatOnOpen=false;
  tsv_gdh.open(tsv_gdh_name);

  //
  std::string tmp_calvin_part;
  tsv_gdh.getHeader("calvin-part",tmp_calvin_part);
  if (tmp_calvin_part!="gdh") {
    printf("#%% not a gdh file.\n");
    return CL_ERR;
  }

  //
  CL_Gdh* gdh;
  std::string tmp_gdh_number;
  std::string tmp_name;
  std::string tmp_type;
  int tmp_byte_len;
  std::string tmp_value;
  int param_import_cnt;

  while (tsv_gdh.nextLevel(0)==affx::TSV_OK) {
    //
    param_import_cnt=0;
    tsv_gdh.get(0,"gdh_number",tmp_gdh_number);
    //
    gdh=getGdhByNumber(tmp_gdh_number);

    while (tsv_gdh.nextLevel(1)==affx::TSV_OK) {
      param_import_cnt++;
      tsv_gdh.get(1,"name",tmp_name);
      tsv_gdh.get(1,"type",tmp_type);
      tsv_gdh.get(1,"byte_len",tmp_byte_len);
      tsv_gdh.get(1,"value",tmp_value);
      //
//      printf("%s:%s = (%s) '%s'/%d\n",tmp_gdh_number.c_str(),
//             tmp_name.c_str(),tmp_type.c_str(),
//             tmp_value.c_str(),tmp_byte_len);

      // @todo fix the double "setWstr" "setNstr".
      //       just want to get it out.
      // the "special" four headers.
      if (tmp_name=="datatype") {
        gdh->setDataType(tmp_value);
        continue;
      }
      if (tmp_name=="uuid") {
        gdh->setUuid(tmp_value);
        continue;
      }
      if (tmp_name=="datetime") {
        gdh->setDateTime(tmp_value);
        continue;
      }
      if (tmp_name=="locale") {
        gdh->setLocale(tmp_value);
        continue;
      }
      // normal params
      // ...
      gdh->addParam(tmp_name,tmp_type,tmp_value,tmp_byte_len);
    }
    printf("#%% imported GDH '%s' with %d params.\n",
           tmp_gdh_number.c_str(),param_import_cnt);
  }  
  //
  return CL_OK;
}

#endif
