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
// affy/sdk/file5/File5_Tsv.cpp ---
//
// $Id: File5_Tsv.cpp,v 1.41 2009-10-28 19:06:41 csugne Exp $
//

//
#include "file5/File5_Tsv.h"
//
#include "file5/File5.h"
#include "file5/File5_Vector.h"
//
#include "file/TsvFile/TsvFile.h"
//
#include <cassert>
#include <stdio.h>
//

#define TSV_LINE_CLVL     "tsv-line-clvl"
#define TSV_LINE_CLVL_IDX "tsv-line-clvl-idx"

//////////

affx::File5_TsvColumn::File5_TsvColumn()
{
  init();
}

affx::File5_TsvColumn::~File5_TsvColumn()
{
  close();
}

int affx::File5_TsvColumn::init()
{
  affx::File5_Vector::init();
  //
  m_kind=FILE5_KIND_TSV_COLUMN;
  m_kind_char='C';
  //
  m_col_name="";
  m_clvl=-1;
  m_cidx=-1;
  //
  m_val_i=0;
  m_val_f=0;
  m_val_d=0;
  m_val_string="";
  //
  return 0;
}

char affx::File5_TsvColumn::file5_kind_char() {
  return 'C';
}

//////////

int affx::File5_TsvColumn::flush()
{
  return affx::File5_Vector::flush();
}

//////////

#define TSV_FIELD_SET(SUFFIX1,SUFFIX2,TYPE)          \
  int affx::File5_TsvColumn::set_##SUFFIX1(TYPE val) \
  {                                                  \
    m_val_##SUFFIX2=val;                             \
    return affx::FILE5_OK;                           \
  }

//
TSV_FIELD_SET(c,i,char);
TSV_FIELD_SET(s,i,short);
TSV_FIELD_SET(i,i,int);
TSV_FIELD_SET(f,f,float);
TSV_FIELD_SET(d,d,double);
TSV_FIELD_SET(string,string,const std::string&);

#undef TSV_FIELD_SET

//////////

int affx::File5_TsvColumn::push_back_val()
{
  //printf("tsv_col(%d,%d)->write_col_val(%d)\n",m_clvl,m_cidx,idx);

  switch (m_dtype) {
  case FILE5_DTYPE_CHAR:
    return push_back_c(m_val_i);
    break;
  case FILE5_DTYPE_SHORT:
    return push_back_s(m_val_i);
    break;
  case FILE5_DTYPE_INT:
    return push_back_i(m_val_i);
    break;
  case FILE5_DTYPE_FLOAT:
    return push_back_f(m_val_f);
    break;
  case FILE5_DTYPE_DOUBLE:
    return push_back_d(m_val_d);
    break;
  case FILE5_DTYPE_STRING:
    return push_back_string(m_val_string);
    break;
  default:
    FILE5_ABORT("Unhandled case.");
    return affx::FILE5_ERR;
    break;
  }
}

//////////

int
affx::File5_TsvColumn::open(const std::string& tsv_name)
{
  int rv;
  // open the base vector
  rv=affx::File5_Vector::open(tsv_name,affx::FILE5_DTYPE_ANY,affx::FILE5_OPEN);
  // should test that this is set.
  // set_file5_kind_string("file5-tsv-column");
  // init the meta-data
  attrib_get("tsv-col-name",&m_col_name);
  attrib_get("tsv-clvl",&m_clvl);
  attrib_get("tsv-cidx",&m_cidx);
  //
  return 0;
}

int
affx::File5_TsvColumn::open(const std::string& tsv_name,
                            const std::string& col_name,
                            int clvl,
                            int cidx,
                            affx::File5_dtype_t dtype,
                            int flags)
{
  int rv;
  // open the base vector
  rv=affx::File5_Vector::open(tsv_name,dtype,flags);
  //
  set_file5_kind_string("file5-tsv-column");
  // init the meta-data
  m_col_name=col_name;
  m_clvl=clvl;
  m_cidx=cidx;
  // save the meta-data for later.
  attrib_set("tsv-col-name",m_col_name);
  attrib_set("tsv-clvl",clvl);
  attrib_set("tsv-cidx",cidx);
  //
  return 0;
}

int affx::File5_TsvColumn::close()
{
  //
  affx::File5_Vector::close();
  //
  return 0;
}

std::string affx::File5_TsvColumn::getColumnName()
{
  return m_col_name;
}

int
affx::File5_TsvColumn::dump()
{
  printf("#%3d:%3d: %p: type=%2d name='%s'\n",m_clvl,m_cidx,this,m_dtype,m_col_name.c_str());
  //
  return 0;
}

//////////

affx::File5_Tsv::File5_Tsv()
{
  init();
}

int affx::File5_Tsv::init()
{
  File5_Group::init();
  //
  m_kind=FILE5_KIND_TSV;
  m_kind_char='T';
  //
  m_line_clvl=NULL;
  m_line_clvl_idx=NULL;
  //
  m_lineNum=-1;
  //
  return 0;
}

char affx::File5_Tsv::file5_kind_char() {
  return 'T';
}

affx::File5_Tsv::~File5_Tsv()
{
  close();
}

// This works as long as it is opened before data is opened.
void
affx::File5_Tsv::openLineIndexes()
{
  int flags=affx::FILE5_CREATE|affx::FILE5_OPEN;
  //
  if (m_line_clvl==NULL) {
    m_line_clvl=this->openVector(TSV_LINE_CLVL,affx::FILE5_DTYPE_INT,flags);
    if (!m_line_clvl->is_open()) {
      FILE5_ABORT("File5_Tsv: cant open: " TSV_LINE_CLVL);
    }
  }
  //
  if (m_line_clvl_idx==NULL) {
    m_line_clvl_idx=this->openVector(TSV_LINE_CLVL_IDX,affx::FILE5_DTYPE_INT,flags);
    if (!m_line_clvl_idx->is_open()) {
      FILE5_ABORT("File5_Tsv: cant open: " TSV_LINE_CLVL_IDX);
    }
  }
}

int affx::File5_Tsv::open(const std::string& name,int flags)
{
  close();
  rewind();

  int rv;

  if ((flags&affx::FILE5_REPLACE)==affx::FILE5_REPLACE) {
    flags|=affx::FILE5_CREATE;
  }

  //
  rv=affx::File5_Group::open(name,flags);
  if (rv!=0) {
    return rv;
  }
  set_file5_kind_string("file5-tsv");

  // scan for columns.
  int num_obj=get_num_objs();

  for (int idx=0;idx<num_obj;idx++) {
    affx::File5_Object* obj=get_object_by_idx(idx);
    //
#ifdef FILE5_DEBUG_PRINT
    printf("### %2d: %-30s: %-20s\n",idx,obj->m_name.c_str(),obj->file5_kind_string().c_str());
#endif
    // skip items which cant be cast.
    if (obj==NULL) {
      continue;
    }
    //
    if (obj->file5_kind_string()=="file5-tsv-column") {
      std::string col_name;
      int col_clvl;
      int col_cidx;
      obj->attrib_get("tsv-col-name",&col_name);
      obj->attrib_get("tsv-clvl",&col_clvl);
      obj->attrib_get("tsv-cidx",&col_cidx);
      // printf("tsv-open: found: %s %d %d\n",col_name.c_str(),col_clvl,col_cidx);
      //
      File5_TsvColumn* col=new affx::File5_TsvColumn();
      col->setParent(this);
      col->open(obj->m_name);
      registerColumn(col,
                     col->m_clvl,
                     col->m_cidx,
                     col->m_col_name);
    }
    // release the tmp object.
    obj->close();
    delete obj;
  }

  // @todo: set up the last indexes
  // m_clvl_last_idx[clvl]++;

  //
  m_state=affx::FILE5_STATE_OPEN;
  return 0;
}

int affx::File5_Tsv::close()
{
  flush();

  // close the indexes if they are open
  if (m_line_clvl!=NULL) {
    m_line_clvl->close();
    delete m_line_clvl;
    m_line_clvl=NULL;
  }
  if (m_line_clvl_idx!=NULL) {
    m_line_clvl_idx->close();
    delete m_line_clvl_idx;
    m_line_clvl_idx=NULL;
  }

  // walk the columns and free them...
  for (int clvl=0;clvl<m_columns.size();clvl++) {
    for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
      affx::File5_TsvColumn* col=m_columns[clvl][cidx];
      if (col!=NULL) {
        col->close();
        delete col;
        m_columns[clvl][cidx]=NULL;
      }
    }
    m_columns[clvl].resize(0);
  }
  m_columns.resize(0);
  //
  affx::File5_Group::close();
  //
  return 0;
}

/////

int affx::File5_Tsv::getLevelCount() {
  return m_columns.size();
}

int affx::File5_Tsv::getLineCount() {
  // if we have an index, that tells us the total number of lines.
  if (m_line_clvl_idx!=NULL) {
    return m_line_clvl_idx->size();
  }
  // if not, we use the number of items
  affx::File5_TsvColumn* colptr=getColumnPtr(0,0);
  if (colptr!=NULL) {
    return colptr->size();
  }
  // dont have any data...
  return -1;
}

int affx::File5_Tsv::getColumnCount(int clvl) {
  if (m_columns.size() == 0) {
    return 0;
  }
  else if (clvl>=m_columns.size()) {
    return -1;
  }
  return m_columns[clvl].size();
}

int affx::File5_Tsv::getColumnIdx(int clvl,const std::string& cname) {
  if ((clvl<0) || (clvl>=m_columns.size())) {
    return -1;
  }
  std::map<std::string, int>::iterator it=m_column_name_map[clvl].find(cname);
  if (it==m_column_name_map[clvl].end()) {
    // not found
    return -1;
  }
  return it->second;
}

/////

void affx::File5_Tsv::resize(size_t new_size) {
  // the problem is what to do about the 2nd level columns.
  // really this should translate the size as a line number
  // and adjust the columns correctly.
  // for the moment we just punt and say we cant have them.

  if (m_columns.size()>1) {
    FILE5_ABORT("File5_Tsv::resize: Cant resize a multi-level File5_Tsv.");
  }

  if (m_line_clvl!=NULL) {
    m_line_clvl->resize(new_size);
  }
  
  //
  for (int clvl=0;clvl<m_columns.size();clvl++) {
    for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
      affx::File5_TsvColumn* col=m_columns[clvl][cidx];
      if (col!=NULL) {
        col->resize(new_size);
      }
    }
  }
}

/////

affx::File5_TsvColumn* affx::File5_Tsv::getColumnPtr(int clvl,int cidx)
{
  if ((clvl<0) || (clvl>=m_columns.size()) ||
      (cidx<0) || (cidx>=m_columns[clvl].size())) {
    return NULL;
  }
  return m_columns[clvl][cidx];
}

affx::File5_TsvColumn* affx::File5_Tsv::getColumnPtr(int clvl,const std::string& cname) {
  return getColumnPtr(clvl,getColumnIdx(clvl,cname));
}

/////

affx::File5_dtype_t affx::File5_Tsv::getColumnDtype(int clvl,int cidx) {
  affx::File5_TsvColumn* col=getColumnPtr(clvl,cidx);
  if (col==NULL) {
    return affx::FILE5_DTYPE_UNKNOWN;
  }
  return col->file5_dtype();
}

affx::File5_dtype_t affx::File5_Tsv::getColumnDtype(int clvl,const std::string& cname) {
  return getColumnDtype(clvl,getColumnIdx(clvl,cname));
}

/////

int affx::File5_Tsv::getColumnName(int clvl,int cidx,std::string* cname) {
  affx::File5_TsvColumn* col=getColumnPtr(clvl,cidx);
  if (col==NULL) {
    return -1;
  }
  //
  *cname=col->getColumnName();
  return 0;
}

////

int
affx::File5_Tsv::defineColumn(const int clvl,
                              const int cidx,
                              const std::string& cname,
                              const File5_dtype_t ctype,
                              int str_size)
{
  // what should it be stored in the file as?
  char name_buf[100];
  sprintf(name_buf,"tsv-col-%03d-%03d",clvl,cidx);

#ifdef FILE5_DEBUG_PRINT
  if (cidx==0) {
    printf("### \n");
  }
  printf("### defineColumn: '%s': cname='%s' type=%d size=%d\n",
         name_buf,cname.c_str(),ctype,str_size);
#endif

  // create it
  File5_TsvColumn* col=new File5_TsvColumn();
  col->setParent(this);
  col->setOptStringSize(str_size);
  col->m_col_name=cname;
  col->open(name_buf,cname,clvl,cidx,ctype,affx::FILE5_REPLACE);
  registerColumn(col,clvl,cidx,cname);

  //
  return 0;
}

int
affx::File5_Tsv::defineColumn(const int clvl,
                              const int cidx,
                              const std::string& cname,
                              const File5_dtype_t ctype)
{
  return defineColumn(clvl,cidx,cname,ctype,-1);
}

int
affx::File5_Tsv::defineStringColumn(const int clvl,
                              const int cidx,
                              const std::string& cname,
                              int size)
{
  return defineColumn(clvl,cidx,cname,affx::FILE5_DTYPE_STRING,size);
}

int
affx::File5_Tsv::registerColumn(affx::File5_TsvColumn* col,
                                const int clvl,
                                int cidx,
                                const std::string& cname)
{
  // resize m_columns needed
  if (m_columns.size()<=clvl) {
    m_columns.resize(clvl+1);
    m_column_name_map.resize(clvl+1);
    m_clvl_last_idx.resize(clvl+1);
  }
  if (m_columns[clvl].size()<=cidx) {
    m_columns[clvl].resize(cidx+1,NULL);
  }

  //
  if (m_columns.size()>1) {
    openLineIndexes();
  }

  // @todo: method to remove old column.
  // discard old column if there.
  if (m_columns[clvl][cidx]!=NULL) {
#ifdef FILE5_DEBUG_PRINT
    printf("### File5_Tsv::registerColumn(%d,%d,'%s'): replacing!\n",
           clvl,cidx,m_columns[clvl][cidx]->m_col_name.c_str());
#endif
    m_columns[clvl][cidx]->close();
    delete m_columns[clvl][cidx];
  }

  // remember the mapping
  m_columns[clvl][cidx]=col;
  m_column_name_map[clvl][cname]=cidx;

  //
  return 0;
}

int
affx::File5_Tsv::dump()
{
  for (int clvl=0;clvl<m_columns.size();clvl++) {
    for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
      File5_TsvColumn* colptr=m_columns[clvl][cidx];
      printf("#%3d:%3d: ptr=%p\n",clvl,cidx,colptr);
      if (colptr!=NULL) {
        colptr->dump();
      }
    }
  }

  //
  return 0;
}

//////////

int
affx::File5_Tsv::writeLevel(int clvl)
{
  if (m_line_clvl!=NULL) {
    // set the map of line_no -> level.
    m_line_clvl->push_back_i(clvl);
    // set the map of line_no -> vector_idx
    int vec_idx;
    vec_idx=m_clvl_last_idx[clvl];
    m_line_clvl_idx->push_back_i(vec_idx);
  }

  //
  for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
    File5_TsvColumn* col=m_columns[clvl][cidx];
    if (col!=NULL) {
      col->push_back_val();
    }
    else {
      // @todo: dont print a warning.
      // printf("writeLevel(%d): cidx %d is NULL.\n",clvl,cidx);
    }
  }
  // last index used for this level
  m_clvl_last_idx[clvl]++;
  // bump to next line.
  m_lineNum++;
  //
  return 0;
}

int
affx::File5_Tsv::writeLevel(int clvl, int cidx)
{
  if (m_line_clvl!=NULL) {
    // set the map of line_no -> level.
    m_line_clvl->push_back_i(clvl);
    // set the map of line_no -> vector_idx
    int vec_idx;
    vec_idx=m_clvl_last_idx[clvl];
    m_line_clvl_idx->push_back_i(vec_idx);
  }

  //
  //for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
    File5_TsvColumn* col=m_columns[clvl][cidx];
    if (col!=NULL) {
      col->push_back_val();
    }
    else {
      // @todo: dont print a warning.
      // printf("writeLevel(%d): cidx %d is NULL.\n",clvl,cidx);
    }
  //}
  // last index used for this level
  m_clvl_last_idx[clvl]++;
  // bump to next line.
  m_lineNum++;
  //
  return 0;
}

/// @brief     Copy the "format" of the file
/// @param     f_tsv     the tsv to copy from
void
affx::File5_Tsv::copyFormat(affx::File5_Tsv& tsv5)
{
  //
  for (int clvl=0;clvl<tsv5.getLevelCount();clvl++) {
    for (int cidx=0;cidx<tsv5.getColumnCount(clvl);cidx++) {
      std::string cname;
      tsv5.getColumnName(clvl,cidx,&cname);

      affx::File5_TsvColumn* column = tsv5.getColumnPtr(clvl,cidx);
      affx::File5_dtype_t type = column->getDType();
      int size = column->getOptStringSize();

      defineColumn(clvl,cidx,cname,type,size);
    }
  }
}

/////

int affx::File5_Tsv::getAsStr(int clvl,int cidx,std::string* str)
{
  affx::File5_TsvColumn* col=getColumnPtr(clvl,cidx);
  if (col==NULL) {
    return -1;
  }
  col->getAsStr(m_lineNum,str);
  int lidx;
  if (m_line_clvl_idx==NULL) {
    lidx=m_lineNum;
  }
  else {
    m_line_clvl_idx->get(m_lineNum,&lidx);
  }
  return col->getAsStr(lidx,str);
}

/////

#define TSV_GET_BODY() {                              \
    affx::File5_TsvColumn* col=getColumnPtr(clvl,cidx);  \
    if (col==NULL) {                                  \
      return -1;                                      \
    }                                                 \
    int lidx;                                         \
    if (m_line_clvl_idx==NULL) {                      \
      lidx=m_lineNum;                                 \
    }                                                 \
    else {                                            \
      m_line_clvl_idx->get(m_lineNum,&lidx);          \
    }                                                 \
    return col->get(lidx,val);                        \
  }

#define TSV_GET_FUNC(TYPE)                                              \
  int affx::File5_Tsv::get(int clvl,int cidx,TYPE* val)                 \
  {                                                                     \
    TSV_GET_BODY();                                                     \
  }                                                                     \
  int affx::File5_Tsv::get(int clvl,const std::string& cidx,TYPE* val)  \
  {                                                                     \
    TSV_GET_BODY();                                                     \
  }

#define TSV_GET_LINE_BODY() {                         \
    affx::File5_TsvColumn* col=m_columns[clvl][cidx]; \
    return col->get(lineIx,val);                      \
    }

#define TSV_GET_LINE_FUNC(TYPE)                                            \
    int affx::File5_Tsv::getLine(int clvl,int cidx, int lineIx, TYPE* val) \
    {                                                                      \
        TSV_GET_LINE_BODY();                                               \
    }                                                                      \


TSV_GET_FUNC(std::string);
TSV_GET_FUNC(char);
TSV_GET_FUNC(int);
TSV_GET_FUNC(float);
TSV_GET_FUNC(double);

TSV_GET_LINE_FUNC(float);

#undef TSV_GET_FUNC
#undef TSV_GET_BODY

/////

#define TSV_SET_BODY(SUFFIX) {                        \
    affx::File5_TsvColumn* col=getColumnPtr(clvl,cidx);  \
    if (col==NULL) {                                  \
      return -1;                                      \
    }                                                 \
    return col->set_##SUFFIX(val);                    \
  }

#define TSV_SET_FUNC(SUFFIX,TYPE)                                       \
  int affx::File5_Tsv::set_##SUFFIX(int clvl,int cidx,TYPE val)         \
  {                                                                     \
    TSV_SET_BODY(SUFFIX);                                               \
  }                                                                     \
  int affx::File5_Tsv::set_##SUFFIX(int clvl,const std::string& cidx,TYPE val) \
  {                                                                     \
    TSV_SET_BODY(SUFFIX);                                               \
  }

TSV_SET_FUNC(c,char);
TSV_SET_FUNC(i,int);
TSV_SET_FUNC(f,float);
TSV_SET_FUNC(d,double);
TSV_SET_FUNC(string,const std::string&);

#undef TSV_SET_FUNC
#undef TSV_SET_BODY

/////

int affx::File5_Tsv::flush()
{
  if (m_state!=affx::FILE5_STATE_OPEN) {
    return 0;
  }

  for (int clvl=0;clvl<m_columns.size();clvl++) {
    for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
      affx::File5_TsvColumn* col=m_columns[clvl][cidx];
      if (col!=NULL) {
        col->flush();
      }
    }
  }

  return 0;
}

//
int affx::File5_Tsv::rewind()
{
  m_lineNum=-1; // point to
  return affx::FILE5_OK;
}

int affx::File5_Tsv::lineNum()
{
  return m_lineNum;
}

int affx::File5_Tsv::lineLevel()
{
  int clvl;
  // if not set, everything is "0" level
  if (m_line_clvl==NULL) {
    return 0;
  }
  m_line_clvl->get(m_lineNum,&clvl);
  return clvl;
}

int affx::File5_Tsv::nextLevel(int clvl)
{
  // @todo: this is only for a single level.
  return nextLine();
}

// No binding so just inc
int affx::File5_Tsv::nextLine()
{
  // advance the current line.
  m_lineNum++;

  // and check if we looking at a line.

  if (m_line_clvl==NULL) {
    // dont have a m_line_clvl index, so use the first column
    // as a proxy for the size.
    // do we have a [0][0] column?
    if ((m_columns.size()==0) || (m_columns[0].size()==0)) {
      return affx::FILE5_EOF;
    }
    // not off the end of this column?
    if (m_lineNum<m_columns[0][0]->m_vec_fill_idx) {
      return affx::FILE5_OK;
    }
  }
  else {
    if (m_lineNum<m_line_clvl->m_vec_fill_idx) {
      return affx::FILE5_OK;
    }
  }

  // backup
  m_lineNum--;
  return affx::FILE5_EOF;
}


//////////

int affx::File5_Tsv::gotoLine(int line)
{
  m_lineNum=line;
  return affx::FILE5_OK;
}

int affx::File5_Tsv::findFirst(int clvl,const std::string& cname,const std::string& val)
{
  return findFirst(clvl,getColumnIdx(clvl,cname),val);
}

int affx::File5_Tsv::findFirst(int clvl,int cidx,const std::string& val)
{
  std::string tmpval;

  // @todo: BAD! linear scan.
  rewind();
  while (nextLine()==affx::FILE5_OK) {
    getAsStr(0,cidx,&tmpval);
    // debugging
    // printf("%3d: ",lineNum()); printfRow();
    //
    if (tmpval==val) {
      return affx::FILE5_OK;
    }
  }
  return affx::FILE5_ERR;
}

void affx::File5_Tsv::printfTabs(int cnt)
{
  for (int i=0;i<cnt;i++) {
    printf("\t");
  }
}

int affx::File5_Tsv::printfHeader()
{
  for (int clvl=0;clvl<m_columns.size();clvl++) {
    printfTabs(clvl);
    //
    for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
      if (cidx!=0) {
        printf("\t");
      }
      File5_TsvColumn* colptr=m_columns[clvl][cidx];
      if (colptr==NULL) {
        printf("NULL");
      }
      else {
        printf("%s",colptr->m_col_name.c_str());
      }
    }
    printf("\n");
  }
  return affx::FILE5_OK;
}

int affx::File5_Tsv::printfRow()
{
  int clvl=0;
  std::string val;

  for (int cidx=0;cidx<m_columns[clvl].size();cidx++) {
    if (cidx!=0) {
      printf("\t");
    }
    getAsStr(0,cidx,&val);
    printf("%s",val.c_str());
  }
  printf("\n");
  return affx::FILE5_OK;
}


/// @brief     Look inside the file and get the count of rows of a tsv
/// @param     file_name 
/// @param     tsv_name  
/// @return    Number of rows, -1 if there is any error.

int affx::File5_Tsv::getFileTsvLineCount(const std::string& file_name,
                                         const std::string& tsv_name)
{
  affx::File5_File* file5;
  affx::File5_Tsv* tsv5;
  int rowCount=-1;

  //
  file5=new affx::File5_File();
  if (file5->open(file_name,affx::FILE5_OPEN_RO)==0) {
    // "/tsv" is the same as "/tsv" because the tsv is the root group.
    if (file5->name_exists(tsv_name)) {
      tsv5=file5->openTsv(tsv_name,affx::FILE5_OPEN);
      if (tsv5->is_open()) {
        rowCount=tsv5->getLineCount();
      }
      tsv5->close();
      delete tsv5;
    }
  }
  file5->close();
  delete file5;

  return rowCount;
}

/// @brief     Look inside the file and get the count of rows of a tsv
/// @param     file_name 
/// @param	   group_name	
/// @param     tsv_name  
/// @return    Number of rows, -1 if there is any error.

int affx::File5_Tsv::getFileTsvLineCount(const std::string& file_name,
										 const std::string& group_name,
                                         const std::string& tsv_name)
{
  affx::File5_File* file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5;
  int rowCount=-1;

  //
  file5=new affx::File5_File();
  if (file5->open(file_name,affx::FILE5_OPEN_RO)==0) {
	if (file5->name_exists(group_name))
	{
	  group5 = file5->openGroup(group_name, affx::FILE5_OPEN_RO);
	  if (group5 != NULL)	{	
	    if (group5->name_exists(tsv_name)) {
	      tsv5=group5->openTsv(tsv_name,affx::FILE5_OPEN_RO);
	      if (tsv5->is_open()) {
	        rowCount=tsv5->getLineCount();
	      }
	      tsv5->close();
	      delete tsv5;
	    }
	    group5->close();
	    delete group5;
	  }
	}
  }
  file5->close();
  delete file5;

  return rowCount;
}
