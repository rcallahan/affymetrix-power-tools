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
// CalvinLite/CL_DataSet.cpp ---
//
// $Id: CL_DataSet.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

//
#include "calvinlite/CL_DataSet.h"
//
#include "calvinlite/CL_File.h"
#include "calvinlite/CL_util.h"
//
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

CL_DataSet::CL_DataSet()
{
  init();
  clear();
}

CL_DataSet::~CL_DataSet()
{
  clear();
}

void CL_DataSet::init()
{
  CL_Object::init();
  // in the same order as the .h
  m_parent_file=NULL;
  m_parent_group=NULL;
  //
  m_fpos_data_start=0;
  m_fpos_data_end=0;
  //m_fpos_pad_bytesize=0;
  //
  m_next_dataset_fpos=0;
  //
  m_row_cnt=0;
  m_row_bytesize=0;
  //
  m_buffer_ptr=NULL;
  m_buffer_size=0;
  //
  m_pushed_dataset=NULL;
  //
  m_dataset_next=NULL;
}

void CL_DataSet::clear()
{
  CL_Object::clear();
  //
  //m_data_fpos_start=0;
  //m_data_fpos_end=0;
  //
  m_next_dataset_fpos=0;
  /// @todo break into clear_data clear_all
  //m_name.clear();
  m_col_vec.clear();
  //
  m_row_cnt=0;
  m_row_bytesize=-1;
  //
  freeBuffer();
}

int CL_DataSet::close()
{
  return CL_OK;
}

//////////

std::string CL_DataSet::getName()
{
  return m_name.Nstr();
}
void CL_DataSet::setName(const std::string& name)
{
  m_name.setNstr(name);
}

//////////

void CL_DataSet::setParentGroup(CL_DataGroup* cl_group)
{
  m_parent_group=cl_group;
}

CL_DataGroup* CL_DataSet::getParentGroup()
{
  return m_parent_group;
}

//////////

int  CL_DataSet::getDataFposStart() const
{
  return m_fpos_data_start;
}
int CL_DataSet::setDataFposStart(int fpos)
{
  m_fpos_data_start=fpos;
  return fpos;
}
int  CL_DataSet::getDataFposEnd() const
{
  return m_fpos_data_end;
}
int CL_DataSet::setDataFposEnd(int fpos)
{
  m_fpos_data_end=fpos;
  return fpos;
}

//////////

int CL_DataSet::getBufferSize()
{
  if (m_buffer_ptr==NULL) {
    allocBuffer();
  }
  return m_buffer_size;
}

void* CL_DataSet::getBufferPtr()
{
  if (m_buffer_ptr==NULL) {
    allocBuffer();
  }
  return m_buffer_ptr;
}

void* CL_DataSet::allocBuffer()
{
  int old_buffer_size=m_buffer_size;
  m_buffer_size=rowCount()*rowBytesize();
  // same size!
  if (m_buffer_size==old_buffer_size) {
    return m_buffer_ptr;
  }
  // zero
  if (m_buffer_size==0) {
    if (m_buffer_ptr!=NULL) {
      free(m_buffer_ptr);
      m_buffer_ptr=NULL;
    }
    return m_buffer_ptr;
  }
  //
  m_buffer_ptr=realloc(m_buffer_ptr,m_buffer_size);
  //printf("allocBuffer():  size=%d ptr=%p\n",m_buffer_size,m_buffer_ptr);
  //
  if (m_buffer_size>old_buffer_size) {
    memset((char*)m_buffer_ptr+old_buffer_size,0,m_buffer_size-old_buffer_size);
  }
  //
  return m_buffer_ptr;
}

void CL_DataSet::freeBuffer()
{
  if (m_buffer_ptr!=NULL) {
    free(m_buffer_ptr);
  }
  m_buffer_ptr=NULL;
  m_buffer_size=0;
}

//////////

int CL_DataSet::newColumn(const std::string& cname,CL_TypeCode_t ctype)
{
  // default the bytelen.
  return newColumn(cname,ctype,CL_codeToBytelen(ctype));
}
int CL_DataSet::newColumn(const std::string& cname,CL_TypeCode_t ctype,int cbytelen)
{
  //printf("newColumn('%s',%d,%d)\n",cname.c_str(),ctype,cbytelen);

  // we gotta have a non-neg bytelen.
  assert(cbytelen>=0);

  // add 4 bytes of overhead for the length - maybe this needs a helper function?
  // but we dont want the user to have to think about the extra 4B.
  if ((ctype==CL_TC_TEXT_ASCII8) ||
      (ctype==CL_TC_TEXT_PLAIN16)||
      (ctype==CL_TC_TEXT_PLAIN8) ||
      (ctype==CL_TC_TEXT_PLAIN16)) {
    cbytelen+=4;
  }

  //
  maybePushData();
  setDataDirty();

  //
  CL_DataSetCol dcol;
  dcol.m_name.setNstr(cname);
  dcol.m_type_code=ctype;
  dcol.m_byte_len=cbytelen;
  dcol.m_src_cidx=-2; // sentinel value for debugging.
  m_col_vec.push_back(dcol);
  //
  computeRowBytesize();
  //
  maybePullData();
  allocBuffer();
  // return the the cidx of this column
  return m_col_vec.size()-1;
}

void CL_DataSet::resizeCols(int cnt)
{
  m_row_bytesize=-1;
  m_col_vec.resize(cnt);
}

CL_DataSetCol* CL_DataSet::getColPtr(int cidx)
{
  if ((cidx<0)||(cidx>=m_col_vec.size())) {
    return NULL;
  }
  return &m_col_vec[cidx];
}

int CL_DataSet::colCount() const
{
  return m_col_vec.size();
}

int CL_DataSet::rowCount() const
{
  return m_row_cnt;
}

int CL_DataSet::setRowCount(int val)
{
  if (m_row_cnt==val) {
    return m_row_cnt;
  }
  setDataDirty();
  m_row_cnt=val;
  allocBuffer();
  return m_row_cnt;
}

int CL_DataSet::rowBytesize()
{
  if (m_row_bytesize==-1) {
    computeRowBytesize();
  }
  return m_row_bytesize;
}

int CL_DataSet::computeRowBytesize()
{
  m_row_bytesize=0;
  for (int cidx=0;cidx<m_col_vec.size();cidx++) {
    m_col_vec[cidx].m_byte_offset=m_row_bytesize;
    m_row_bytesize+=m_col_vec[cidx].m_byte_len;
  }
  return m_row_bytesize;
}

//void CL_DataSet::dirtyMeta() {
//  m_row_bytesize=-1;
//}

int CL_DataSet::getColName(int cidx,std::string* cname)
{
  *cname="";
  if ((cidx>=0)&&(cidx<m_col_vec.size())) {
    *cname=m_col_vec[cidx].m_name.Nstr();
  }
  return CL_OK;
}

unsigned char* CL_DataSet::getValuePtr(int ridx,int cidx)
{
  if ((ridx<0)||(ridx>=m_row_cnt)) {
    return NULL;
  }
  if ((cidx<0)||(cidx>=m_col_vec.size())) {
    return NULL;
  }
  // be sure our buffer data is up to date.
  if (m_row_bytesize==-1) {
    rowBytesize();
  }
  if (m_buffer_ptr==NULL) {
    allocBuffer();
  }

  unsigned char* ptr=(unsigned char*)m_buffer_ptr+(ridx*m_row_bytesize)+m_col_vec[cidx].m_byte_offset;
  return ptr;
}

//////////

int CL_DataSet::get(int ridx,int cidx,int* val)
{
  *val=0;
  unsigned char* ptr=getValuePtr(ridx,cidx);
  if (ptr==NULL) {
    return CL_ERR;
  }

  switch (m_col_vec[cidx].m_type_code) {
  case CL_TC_BYTE    :
    *val=*(char*)ptr;
    break;
  case CL_TC_UBYTE   :
    *val=*(unsigned char*)ptr;
    break;
  case CL_TC_SHORT   :
    *val=((char)ptr[0])<<8|ptr[1];
    break;
  case CL_TC_USHORT  :
    *val=(ptr[0]<<8)|ptr[1];
    break;
  case CL_TC_INT     :
  case CL_TC_UINT    :
    *val=(ptr[0]<<24)|(ptr[1]<<16)|(ptr[2]<<8)|(ptr[3]);
    break;
  default:
    assert(0);
    return CL_ERR;
  }

  //
  return CL_OK;
}

//
int CL_DataSet::get(int ridx,int cidx,float* val)
{
  *val=0.0;
  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }
  unsigned char* ptr=getValuePtr(ridx,cidx);
  if (ptr==NULL) {
    return CL_ERR;
  }
  if (dcol->m_type_code!=CL_TC_FLOAT) {
    return CL_ERR;
  }
  //
  type_punned pun;
  pun.v_uint32=(ptr[0]<<24)|(ptr[1]<<16)|(ptr[2]<<8)|(ptr[3]);
  *val=pun.v_float;
  //
  return CL_OK;
}

int CL_DataSet::get(int ridx,int cidx,std::string* val)
{
  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }

  if (!((dcol->m_type_code==CL_TC_TEXT_ASCII8 )||
        (dcol->m_type_code==CL_TC_TEXT_ASCII16)||
        (dcol->m_type_code==CL_TC_TEXT_PLAIN8 )||
        (dcol->m_type_code==CL_TC_TEXT_PLAIN16))) {
    return CL_ERR;
  }

  unsigned char* valptr=getValuePtr(ridx,cidx);
  if (valptr==NULL) {
    return CL_ERR;
  }

  // how long could it be for this column?
  int str_len_max=dcol->m_byte_len-4; // deduct the first 4B of length count.
  // how long is this string?
  int str_len=CL_mem_read_u4B(valptr);
  // do these values make sense?
  assert(0<=str_len_max);
  assert(0<=str_len);
  assert(str_len<=str_len_max);

  unsigned char* str_ptr;
  str_ptr=valptr+4;
  val->resize(str_len);
  // we should handle the wide to narrow here too.
  memcpy(&((*val)[0]),str_ptr,str_len);

  // maybe this should be optional?
  // strnlen?
  // ust stl to do the find.
  for (int i=0;i<str_len;i++) {
    if ((*val)[i]==0) {
      (*val).resize(i);
      break;
    }
  }

  return CL_OK;
}

int CL_DataSet::getAsString(int ridx,int cidx,std::string* val)
{
  *val="";

  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }
  unsigned char* valptr=getValuePtr(ridx,cidx);
  if (valptr==NULL) {
    return CL_ERR;
  }

  // used when converting non strings to strings.
  char buf[64];
  memset(buf,0,sizeof(buf));

  int tmp_int;
  float tmp_float;

  switch (dcol->m_type_code) {
  case CL_TC_BYTE    :
  case CL_TC_UBYTE   :
  case CL_TC_SHORT   :
  case CL_TC_USHORT  :
  case CL_TC_INT     :
  case CL_TC_UINT    :
    get(ridx,cidx,&tmp_int);
    sprintf(buf,"%d",tmp_int);
    *val=buf;
    break;
  case CL_TC_FLOAT   :
    get(ridx,cidx,&tmp_float);
    sprintf(buf,"%.10f",tmp_float);
    *val=buf;
    break;
  case CL_TC_TEXT_ASCII8  :
  case CL_TC_TEXT_ASCII16 :
  case CL_TC_TEXT_PLAIN8  :
  case CL_TC_TEXT_PLAIN16 :
    return get(ridx,cidx,val);
    break;
  case CL_TC_UNSET   :
    strcpy(buf,"UNSET");
    *val=buf;
    break;
  default:
    strcpy(buf,"DEFAULT");
    *val=buf;
    break;
  }

  return CL_OK;
}

/////

CL_Err_t CL_DataSet::set(int ridx,int cidx,int val)
{
  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }
  unsigned char* ptr=getValuePtr(ridx,cidx);
  if (ptr==NULL) {
    return CL_ERR;
  }
  //
  switch (dcol->m_type_code) {
  case CL_TC_BYTE:
  case CL_TC_UBYTE:
    CL_mem_write_i1B(ptr,val);
    break;
  case CL_TC_SHORT:
  case CL_TC_USHORT:
    CL_mem_write_i2B(ptr,val);
    break;
  case CL_TC_INT:
  case CL_TC_UINT:
    CL_mem_write_i4B(ptr,val);
    break;
  case CL_TC_FLOAT:
    CL_mem_write_float(ptr,val);
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    // this is an invalid usage.
    // should we just abort?
    // assert(0);
    return CL_ERR;
    break;
  default:
    assert(0);
  }
  //
  setDataDirty();
  return CL_OK;
}

CL_Err_t CL_DataSet::setFloat(int ridx,int cidx,float val)
{
  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }
  unsigned char* ptr=getValuePtr(ridx,cidx);
  if (ptr==NULL) {
    return CL_ERR;
  }
  //
  if (dcol->m_type_code!=CL_TC_FLOAT) {
    return CL_ERR;
  }
  //
  CL_mem_write_float(ptr,val);
  //
  setDataDirty();
  return CL_OK;
}

CL_Err_t CL_DataSet::set(int ridx,int cidx,const std::string& val)
{
  CL_DataSetCol* dcol=getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }
  
  if (!((dcol->m_type_code==CL_TC_TEXT_ASCII8 )||
        (dcol->m_type_code==CL_TC_TEXT_ASCII16)||
        (dcol->m_type_code==CL_TC_TEXT_PLAIN8 )||
        (dcol->m_type_code==CL_TC_TEXT_PLAIN16))) {
    // again, should we just abort?
    // assert(0);
    return CL_ERR;
  }
  
  unsigned char* ptr=getValuePtr(ridx,cidx);
  if (ptr==NULL) {
    return CL_ERR;
  }

  // zero out the buffer.
  memset(ptr,0,dcol->m_byte_len);

  int str_len=val.size();
  int str_len_max=dcol->m_byte_len-4;
  assert(str_len_max>=0);
  
  // clip length
  if (str_len_max<str_len) {
    str_len=str_len_max;
  }
  // 
  CL_mem_write_i4B(ptr,str_len);
  memcpy(ptr+4,&val[0],str_len);

  //
  setDataDirty();
  return CL_OK;
}

/////

#define SETFROMARRAY_COPY(_type,_write_func)    \
  for (int i=0;i<arrayLen;i++) {                \
    _type tmp_val=arrayData[i];                 \
    _write_func(row_ptr,tmp_val);               \
    row_ptr+=row_bytesize;                      \
  }

template<typename T1>  
CL_Err_t CL_DataSet__setFromArray(CL_DataSet* ds,
                                  int ridx_start,
                                  int cidx,
                                  T1* arrayData,
                                  int arrayLen)
{
  CL_DataSetCol* dcol=ds->getColPtr(cidx);
  if (dcol==NULL) {
    return CL_ERR;
  }

  if ((ridx_start<0)||
      (ridx_start+arrayLen>ds->rowCount())) {
    return CL_ERR; // better code needed.
  }

  unsigned char* row_ptr=ds->getValuePtr(ridx_start,cidx);
  if (row_ptr==NULL) {
    return CL_ERR;
  }
  int row_bytesize=ds->rowBytesize();

  switch (dcol->m_type_code) {
  case CL_TC_BYTE:
    SETFROMARRAY_COPY(unsigned char,CL_mem_write_i1B);
    break;
  case CL_TC_UBYTE:
    SETFROMARRAY_COPY(unsigned char,CL_mem_write_i1B);
    break;
  case CL_TC_SHORT:
    SETFROMARRAY_COPY(unsigned short,CL_mem_write_i2B);
    break;
  case CL_TC_USHORT:
    SETFROMARRAY_COPY(unsigned short,CL_mem_write_i2B);
    break;
  case CL_TC_INT:
    SETFROMARRAY_COPY(unsigned int,CL_mem_write_i4B);
    break;
  case CL_TC_UINT:
    SETFROMARRAY_COPY(unsigned int,CL_mem_write_i4B);
    break;
  case CL_TC_FLOAT:
    SETFROMARRAY_COPY(float,CL_mem_write_float);
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    return CL_ERR;
    break;
  default:
    assert(0);
  }

  ds->setDataDirty();
  return CL_OK;
}

CL_Err_t CL_DataSet::setFromArray(int ridx_start,
                                  int cidx,
                                  char* arrayData,
                                  int arrayLen)
{
  return CL_DataSet__setFromArray(this,ridx_start,cidx,arrayData,arrayLen);
}
CL_Err_t CL_DataSet::setFromArray(int ridx_start,
                                  int cidx,
                                  int* arrayData,
                                  int arrayLen)
{
  return CL_DataSet__setFromArray(this,ridx_start,cidx,arrayData,arrayLen);
}
CL_Err_t CL_DataSet::setFromArray(int ridx_start,
                                  int cidx,
                                  float* arrayData,
                                  int arrayLen)
{
  return CL_DataSet__setFromArray(this,ridx_start,cidx,arrayData,arrayLen);
}

/////

void CL_DataSet::dump(const std::string& prefix_in)
{
  std::string prefix=prefix_in;

  rowBytesize();
  //
  CL_DUMP_PREFIX("DS  ");
  CL_DUMP_MEMB_STR(m_name);
  CL_DUMP_MEMB_INT(m_row_cnt);
  CL_DUMP_MEMB_INT(m_row_bytesize);
  CL_DUMP_MEMB_INT(m_fpos_data_start);
  CL_DUMP_MEMB_INT(m_fpos_data_end);
  //CL_DUMP_MEMB_INT(m_fpos_pad_bytesize);
  CL_DUMP_MEMB_INT(m_next_dataset_fpos);
  CL_DUMP_MEMB_INT(m_buffer_size);

  for (int i=0;i<m_col_vec.size();i++) {
    m_col_vec[i].dump(prefix);
  }

  CL_DUMP_BREAK();
}

void CL_DataSet::maybePushData()
{
  //printf("maybePushData()\n");
  if (m_buffer_ptr!=NULL) {
    assert(m_pushed_dataset==NULL);
    m_pushed_dataset=new CL_DataSet();
    pushDataInto(m_pushed_dataset);
  }
}

void CL_DataSet::pushDataInto(CL_DataSet* dst)
{
  // printf("pushDataInto(%p)\n",dst);
  //
  dst->clear();

  // copy
  dst->m_row_cnt=m_row_cnt;
  dst->m_row_bytesize=m_row_bytesize;
  dst->m_buffer_size=m_buffer_size;
  dst->m_buffer_ptr=m_buffer_ptr;
  dst->m_col_vec=m_col_vec;

  //dst->dump("push ");

  // clear the memory ptr ourselves, dont call clear.
  m_row_cnt=0;
  m_row_bytesize=-1;
  m_buffer_size=0;
  m_buffer_ptr=NULL;

  // remember order of the columns, in case we shuffle the order around
  for (int cidx=0;cidx<m_col_vec.size();cidx++) {
    m_col_vec[cidx].m_src_cidx=cidx;
  }
}


void CL_DataSet::copyDataFrom(CL_DataSet* src)
{
  CL_DataSet::copyData(this,src);
}

void CL_DataSet::copyData(CL_DataSet* dst,CL_DataSet* src)
{
  //dst->dump("dst");
  //src->dump("src");
  //
  int row_cnt=src->rowCount();
  dst->setRowCount(row_cnt);
  //
  dst->allocBuffer();
  //
  char* src_ptr=(char*)src->m_buffer_ptr;
  char* dst_ptr=(char*)dst->m_buffer_ptr;
  int dst_cidx_size=dst->m_col_vec.size();
  //
  for (int ri=0;ri<row_cnt;ri++) {
    for (int dst_cidx=0;dst_cidx<dst_cidx_size;dst_cidx++) {
      CL_DataSetCol* dst_col=&dst->m_col_vec[dst_cidx];
      int dst_src_cidx=dst_col->m_src_cidx;
      if (dst_src_cidx>=0) {
        CL_DataSetCol* src_col=&src->m_col_vec[dst_src_cidx];
        //
        if (CL_File_debug_flags==1) {
          printf("%3d | %p %3d  | %p %3d | %3d\n",
                 dst_src_cidx,
                 dst_ptr,dst_col->m_byte_offset,
                 src_ptr,src_col->m_byte_offset,
                 dst_col->m_byte_len);
        }
        //
        memcpy((void*)&dst_ptr[dst_col->m_byte_offset],
               (void*)&src_ptr[src_col->m_byte_offset],
               dst_col->m_byte_len);
      }
    }
    // advance to the next row
    src_ptr+=src->m_row_bytesize;
    dst_ptr+=dst->m_row_bytesize;
  }
}

void CL_DataSet::maybePullData()
{
  if (m_pushed_dataset==NULL) {
    return;
  }
  //
  //printf("pullData()\n");
  //
  copyDataFrom(m_pushed_dataset);
  m_pushed_dataset->clear();
  delete m_pushed_dataset;
  m_pushed_dataset=NULL;
}
