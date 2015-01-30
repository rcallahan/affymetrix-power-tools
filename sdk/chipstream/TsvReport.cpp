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
// affy/sdk/chipstream/TsvReport.cpp ---
//
// $Id: TsvReport.cpp,v 1.27 2009-09-22 00:25:02 mspald Exp $
//

//
#include "chipstream/TsvReport.h"
//
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/Util.h"
//
#include <cassert>


#define TSVREPORT_COMMENT_MARKER ""

//
affx::TsvReport::TsvReport() {
  init();
}

affx::TsvReport::~TsvReport() {
  close();
  //
  delete m_tsv;
  m_tsv=NULL;
  //
  // dont free the shared file, someone else owns it.
  if (m_a5_group!=m_a5_shared_group) {
    delete m_a5_group;
  }
  m_a5_group=NULL;
  m_a5_shared_group=NULL;
  //
  delete m_a5_file;
  m_a5_file=NULL;

  // if we were made just to buffer, erase m_header_buffer -- we dont care.
  if (m_is_header_buffer==1) {
    m_header_buffer.resize(0);
  }
  //
  if (m_header_buffer.size()!=0) {
    //printf("### TsvReport::~TsvReport: '%s' has %d unflushed headers.\n",
    //       getFilePath().c_str(), int(m_header_buffer.size()));
    //
    // Err::errAbort("TsvReport::~TsvReport: unflushed headers.");
  }
}

void affx::TsvReport::init() {
  m_format=FMT_UNSET;
  m_precision=5;
  //
  m_tsv=NULL;
  //
  m_a5_tsv=NULL;
  m_a5_shared_group=NULL;
  m_a5_group=NULL;
  m_a5_file=NULL;
  //
  m_is_open=false;
  m_has_guid_header=0;
  m_is_header_buffer=0;
  m_use_default_suffix=1;
  //
  m_header_buffer.resize(0);
}

void affx::TsvReport::close() {
  m_is_open=false;
  //
  if (m_format==affx::TsvReport::FMT_TSV) {
    if (m_tsv!=NULL) {
      m_tsv->close();
      delete m_tsv;
      m_tsv=NULL;
    }
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    // do close our tsv
    if (m_a5_tsv!=NULL) {
      m_a5_tsv->close();
      delete m_a5_tsv;
      m_a5_tsv=NULL;
    }
    // but dont close the shared file.
    if (m_a5_group!=m_a5_shared_group) {
      if (m_a5_group!=NULL) {
        m_a5_group->close();
        delete m_a5_group;
        m_a5_group=NULL;
      }
    }
    //
    if (m_a5_file!=NULL) {
      m_a5_file->close();
      delete m_a5_file;
      m_a5_file=NULL;
    }
  }
  else {
    // do nothing! We might be destroying a report which was never used.
    // these two lines were good at flushing out errors for debugging.
    // printf("### TsvReport::close(): no format set: '%s'\n",getFilePath().c_str());
    // Err::errAbort("TsvReport::close(): no format set for report.");
  }
}

void affx::TsvReport::setA5SharedGroup(affx::File5_Group* shared_group)
{
  if (m_is_open==1) {
    Err::errAbort("TsvReport::setA5SharedGroup(): Shouldnt be setting the shared file while open."); 
  }  
  // dont free the shared file, someone else owns it.
  // remember it.
  m_a5_shared_group=shared_group;
}

bool affx::TsvReport::is_open()
{
  return m_is_open;
}

void affx::TsvReport::setDirPath(const std::string& dirpath)
{
  m_dir_path=dirpath;
}
std::string affx::TsvReport::getFileprefix()
{
  return m_file_prefix;
}
void affx::TsvReport::setFileprefix(const std::string& fileprefix)
{
  m_file_prefix=fileprefix;
}
void affx::TsvReport::setFilename(const std::string& filename)
{
  m_file_name=filename;
}
void affx::TsvReport::setGroupname(const std::string& groupname)
{
  m_group_name=groupname;
}

void affx::TsvReport::setIsHeaderBuffer(int val) {
  m_is_header_buffer=val;
}
void affx::TsvReport::setUseDefaultSuffix(int val) {
  m_use_default_suffix=val;
}

/////

void affx::TsvReport::setPrecision(int precision)
{
  // this isnt a big error, but might as well catch it while we are here.
  if (m_precision<0) {
    APT_ERR_ABORT("precision must be >=0.");
  }
  m_precision=precision;
}

void affx::TsvReport::setPrecision(int clvl,int cidx,int precision) {
  if (m_tsv!=NULL) {
    m_tsv->setPrecision(clvl,cidx,precision);
  }
}

/////

affx::File5_File* affx::TsvReport::openA5File(const std::string& filename, bool replace)
{
  // printf("### TsvReport::openA5File('%s')\n",filename.c_str());
  File5_File* a5_file=new affx::File5_File();
  assert(a5_file!=NULL);
  if(replace)
    a5_file->open(filename,affx::FILE5_REPLACE);
  else
    a5_file->open(filename,affx::FILE5_OPEN|affx::FILE5_CREATE);
  return a5_file;
}

void affx::TsvReport::closeA5File(affx::File5_File*& a5_file)
{
  if (a5_file!=NULL) {
    a5_file->close();
    delete a5_file;
  }
  a5_file=NULL;
}

//////////

void affx::TsvReport::setFormat(affx::TsvReport::TsvReportFmt_t format)
{
  // remember
  m_format=format;
  // if we are just a buffer, we are done.
  if (m_is_header_buffer) {
    return;
  }
  // TsvFile is just in memory at the moment. (No IO)
  if (m_format==affx::TsvReport::FMT_TSV) {
    assert(m_tsv==NULL);
    m_tsv=new affx::TsvFile();
    assert(m_tsv!=NULL);
  }
  // A5-TsvFile goes to disk right away.
  else if (m_format==affx::TsvReport::FMT_A5) {
    // have a shared group?
    if (m_a5_shared_group==NULL) {
      // no, make our own A5 file.
      assert(m_a5_file==NULL);
      m_a5_file=new affx::File5_File();
      assert(m_a5_file!=NULL);
      // printf("### setFormat(%d): '%s'\n",format,getFilePath().c_str());
      m_a5_file->open(getFilePath(),affx::FILE5_REPLACE);
      // and use the root group.
      if (m_group_name=="") {
        m_group_name="/";
      }
      m_a5_group=m_a5_file->openGroup(m_group_name,affx::FILE5_CREATE|affx::FILE5_OPEN);
    }
    else {
      // yes, use it.
      m_a5_group=m_a5_shared_group;
    }

    // we need Fs::basename to make the group name the last part of the filename,
    // if someone puts a path into the filename.
    m_a5_tsv=m_a5_group->openTsv(Fs::basename(m_file_name),affx::FILE5_CREATE);
    assert(m_a5_tsv!=NULL);
  }
  else {
    Err::errAbort("TsvReport::setFormat: bad format");
  }
}

//////////

#define DISPATCH_GET()                            \
  {                                               \
    if (m_format==affx::TsvReport::FMT_TSV) {     \
      assert(m_tsv!=NULL);                        \
      return m_tsv->get(clvl,cidx,*val);          \
    }                                             \
    else if (m_format==affx::TsvReport::FMT_A5) { \
      assert(m_a5_tsv!=NULL);                     \
      return m_a5_tsv->get(clvl,cidx,val);        \
    }                                             \
    else {                                        \
      Err::errAbort("File format not set.");      \
    }                                             \
    return 0;                                     \
  }

int affx::TsvReport::get(int clvl,int cidx,int* val)
{
  DISPATCH_GET();
}

int affx::TsvReport::get(int clvl,int cidx,float* val)
{
  DISPATCH_GET();
}

int affx::TsvReport::get(int clvl,int cidx,double* val)
{
  DISPATCH_GET();
}
int affx::TsvReport::get(int clvl,int cidx,std::string* val)
{
  DISPATCH_GET();
}

//////////

#define DISPATCH_SET(X)                                          \
  {                                                              \
    if (m_format==affx::TsvReport::FMT_TSV) {                    \
      assert(m_tsv!=NULL);                                       \
      return m_tsv->set(clvl,cidx,val);                          \
    }                                                            \
    else if (m_format==affx::TsvReport::FMT_A5) {                \
      assert(m_a5_tsv!=NULL);                                    \
      return m_a5_tsv->set_##X(clvl,cidx,val);                   \
    }                                                            \
    else {                                                       \
      Err::errAbort("Format not set on TsvReport");              \
    }                                                            \
    return 0;                                                    \
  }

//
int affx::TsvReport::set_c(int clvl,int cidx,char val)
{
  DISPATCH_SET(c);
}
int affx::TsvReport::set_i(int clvl,int cidx,int val)
{
  DISPATCH_SET(i);
}
int affx::TsvReport::set_f(int clvl,int cidx,float val)
{
  DISPATCH_SET(f);
}
int affx::TsvReport::set_d(int clvl,int cidx,double val)
{
  DISPATCH_SET(d);
}
int affx::TsvReport::set_string(int clvl,int cidx,const std::string& val)
{
  DISPATCH_SET(string);
}

int affx::TsvReport::defineColumn(const int clvl,
                                  const int cidx,
                                  const std::string& cname,
                                  const File5_dtype_t ctype)
{
  return defineColumn(clvl,cidx,cname,ctype,0);
}

int affx::TsvReport::defineColumn(const int clvl,
                                  const int cidx,
                                  const std::string& cname,
                                  const File5_dtype_t ctype,
                                  int str_size)
{
  // we require a name.
  if (cname=="") {
    Err::errAbort("TsvReport::defineColumn: Column name required.");
  }
  //
  if (m_format==affx::TsvReport::FMT_TSV) {
    assert(m_tsv!=NULL);
    m_tsv->defineColumn(clvl,cidx,cname,affx::TSV_TYPE_UNKNOWN);
    m_tsv->setPrecision(clvl,cidx,m_precision);
    return 0;
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    return m_a5_tsv->defineColumn(clvl,cidx,cname,ctype,str_size);
  }
  else {
    Err::errAbort("File format not set.");
  }
  return 0;
}
// no need for the dtype
int affx::TsvReport::defineStringColumn(const int clvl,
                                        const int cidx,
                                        const std::string& cname,
                                        int str_size)
{
  if (m_format==affx::TsvReport::FMT_TSV) {
    assert(m_tsv!=NULL);
    return m_tsv->defineColumn(clvl,cidx,cname,affx::TSV_TYPE_STRING);
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    assert(m_a5_tsv!=NULL);
    return m_a5_tsv->defineStringColumn(clvl,cidx,cname,str_size);
  }
  else {
    Err::errAbort("");
  }
  return 0;
}

int affx::TsvReport::getColumnCount(int clvl)
{
  if (m_format==affx::TsvReport::FMT_TSV) {
    assert(m_tsv!=NULL);
    return m_tsv->getColumnCount(clvl);
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    assert(m_a5_tsv!=NULL);
    return m_a5_tsv->getColumnCount(clvl);
  }
  else {
    Err::errAbort("TsvReport::getColumnCount(): unset format.");
  }
  return 0;
}

void affx::TsvReport::defineColumns(const std::vector<std::string> &colNames,File5_dtype_t ctype,int str_size)
{
  int cidx_start=getColumnCount(0);
  for (int cidx=0;cidx<colNames.size();cidx++) {
    defineColumn(0,cidx_start+cidx,colNames[cidx],ctype,str_size);
  }
}

//////////

int affx::TsvReport::addHeader(const std::string& key, const std::string& val)
{
  // remember if 
  if (key==TSVREPORT_GUID_HDR) {
    m_has_guid_header=1;
  }

  if ((m_is_header_buffer==1)||(m_format==affx::TsvReport::FMT_UNSET)) {
    m_header_buffer.push_back(std::pair<std::string,std::string>(key,val));
    //printf("### addHeader: '%s': buffering: '%s'='%s'\n",getFilePath().c_str(),key.c_str(),val.c_str());
  }
  else if (m_format==affx::TsvReport::FMT_TSV) {
    return m_tsv->addHeader(key,val);
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    return m_a5_tsv->addHeader(key,val);
  }
  else {
    Err::errAbort("TsvReport::addHeader: Format or buffer not set.");
  }
  return 0;
}

int affx::TsvReport::addHeaderComment(const std::string& comment)
{
  if (m_is_header_buffer==1) {
    m_header_buffer.push_back(std::pair<std::string,std::string>(TSVREPORT_COMMENT_MARKER,comment));
    // printf("### addHeaderComment: '%s': buffering: '%s'\n",getFilePath().c_str(),comment.c_str());
  }
  else if (m_format==affx::TsvReport::FMT_TSV) {
    return m_tsv->addHeaderComment(comment);
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    // return m_a5_tsv->addHeaderComment(comment);
  }
  else {
    Err::errAbort("TsvReport::addHeaderComment: Format or buffer not set.");
  }
  return 0;
}

int affx::TsvReport::addHeaderComments(const std::vector<std::string>& comments) {
  for (int i=0;i<comments.size();i++) {
    addHeaderComment(comments[i]);
  }
  return 0;
}
//////////

int affx::TsvReport::writeLevel(int lvl)
{
  if (!m_is_open) {
    Err::errAbort("TsvReport::writeLevel: attempt to write to an unopened TsvReport. " + m_file_name);
  }
  //
  if (m_format==affx::TsvReport::FMT_TSV) {
    return m_tsv->writeLevel(lvl);
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    return m_a5_tsv->writeLevel(lvl);
  }
  else {
    Err::errAbort("TsvReport::writeLevel: bad or unknown format.");
  }
  return 0;
}

int affx::TsvReport::writeTsv_v1(const std::string& filename,bool abort_on_err)
{
  ensureGuidHeader();
  //
  if (m_format==affx::TsvReport::FMT_TSV) {
    // TsvFile wants the headers before the file is opened for writing.
    flushHeaderBuffer();
    int rv=m_tsv->writeTsv_v1(filename);
    if (rv==affx::TSV_OK) {
      m_is_open=true;
      return affx::TSV_OK;
    }
    else {
      if (abort_on_err) {
        Err::errAbort("TsvReport::writeTsv_v1: unable to open '"+filename+"' for writing'");
      }
      return affx::TSV_ERR_FILEIO;
    }
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    // A5_Tsv wants the headers after the file is opened (which we did in prepare)
    flushHeaderBuffer();
    m_is_open=true;
  }
  else {
    Err::errAbort("TsvReport::writeTsv_v1: bad or unset format.");
  }
  //
  return affx::TSV_OK;
}
int affx::TsvReport::writeTsv_v2(const std::string& filename,bool abort_on_error)
{
  ensureGuidHeader();
  //
  if (m_format==affx::TsvReport::FMT_TSV) {
    flushHeaderBuffer();
    int rv=m_tsv->writeTsv(filename);
    if (rv!=affx::TSV_OK) {
      Err::errAbort("TsvReport::writeTsv_v2: Cant open '"+filename+"' to write.");
    }
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    flushHeaderBuffer();
  }
  else {
    Err::errAbort("TsvReport::writeTsv_v2: bad or unset format.");
  }
  //
  m_is_open=true;
  //
  return affx::TSV_OK;
}

void affx::TsvReport::copyOptionsTo(affx::TsvReport& tsv)
{
  // force the call to "setFormat()"
  // tsv.m_format=m_format;
  if (tsv.m_dir_path == "")
    tsv.m_dir_path=m_dir_path;
  if (tsv.m_file_prefix == "")
    tsv.m_file_prefix=m_file_prefix;
  if (tsv.m_file_name == "") 
    tsv.m_file_name=m_file_name;
  if (tsv.m_group_name == "")
    tsv.m_group_name=m_group_name;
  //
  tsv.m_precision=m_precision;
  // if there is a shared group, remember to use it.
  if (tsv.m_a5_shared_group == NULL)
    tsv.m_a5_shared_group=m_a5_shared_group;
}

std::string affx::TsvReport::getFilePath() {
  std::string file_path="";

  if (m_dir_path!="") {
    file_path = m_dir_path;
  }
  //
  //if (m_file_name=="") {
  //  Err::errAbort("TsvReport::getFilePath: m_file_name unset");
  //}
  file_path=Fs::join(file_path,m_file_name);
  //
  if ((m_use_default_suffix==1)&&(m_format!=affx::TsvReport::FMT_UNSET)) {
    file_path=file_path+getFmtSuffix();
  }
  //
  return file_path;
}

std::string affx::TsvReport::getFmtSuffix()
{
  if (m_use_default_suffix==0) {
    return "";
  }
  if (m_format==affx::TsvReport::FMT_TSV) {
    return ".txt";
  }
  else if (m_format==affx::TsvReport::FMT_A5) {
    return ".a5";
  }
  else if (m_format==affx::TsvReport::FMT_UNSET) {
    return ".DEBUG-UNSET";
  }
  else {
    Err::errAbort("TsvReport::getFileSuffix: no format set.");
    return "";
  }
}

void affx::TsvReport::flushHeaderBuffer()
{
  if (m_format==FMT_UNSET) {
    // @todo: take this out sometime later.
    // you shouldnt see this in correct usage.
    // printf("### TsvReport::flushHeaderBuffer(): '%s': flushing w/o format.\n",getFilePath().c_str());
    Err::errAbort("flushHeaderBuffer() called without a report format set");
    return;
  }
  if (m_is_header_buffer==1) {
    // @todo: take this out sometime later.
    // you shouldnt see this in correct usage.
    //printf("### TsvReport::flushHeaderBuffer(): '%s': attempting to flush header buffer.\n",getFilePath().c_str());
    Err::errAbort("flushHeaderBuffer() called withm_is_header_buffer set to 1");
    return;
  }
  // copy the buffered headers back out.
  for (int i=0;i<m_header_buffer.size();i++) {
    if (m_header_buffer[i].first==TSVREPORT_COMMENT_MARKER) {
      addHeaderComment(m_header_buffer[i].second);
    }
    else {
      addHeader(m_header_buffer[i].first,m_header_buffer[i].second);
    }
  }
  // empty the buffer
  m_header_buffer.resize(0);
}

void affx::TsvReport::addHeadersFrom(const std::vector<keyval_t>& header_vec)
{
  for (int i=0;i<header_vec.size();i++) {
    addHeader(header_vec[i].first,header_vec[i].second);
  }
};

void affx::TsvReport::ensureGuidHeader()
{
  if (m_has_guid_header==0) {
    std::string guid=affxutil::Guid::GenerateNewGuid();
    addHeader(TSVREPORT_GUID_HDR,guid);
  }
}
