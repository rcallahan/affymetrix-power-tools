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
// CalvinLite/CL_Object.cpp ---
// 
// $Id: CL_Object.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

//
#include "calvinlite/CL_Object.h"
//
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_string.h"
#include "calvinlite/CL_util.h"
//
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

CL_Object::CL_Object()
{
  init();
  clear();
}

void CL_Object::init()
{
  m_parent_file=NULL;
  m_err_should_throw=true;
}

void CL_Object::clear()
{
  clearErr();
  clearFpos();
  clearOrigFpos();
  //
  m_bytesize_dirty=true;
  m_bytesize=0;
  m_bytesize_pad=0;
  m_data_dirty=false;
}

//////////

CL_Err_t CL_Object::setErr(CL_Err_t err_num,
                           const std::string& err_msg,
                           const std::string& err_path)
{
  m_err_num=err_num;
  m_err_msg=err_msg;
  m_err_path=err_path;

  // Waah! we have an error and we should complain...
  if ((m_err_num!=CL_OK) &&  (m_err_should_throw==true)) {
    printErrMsg();
    throw CalvinLiteException(m_err_num,m_err_msg,m_err_path);
  }
  return m_err_num;
}

bool CL_Object::throwOnErr(bool val) {
  m_err_should_throw=val;
  return m_err_should_throw;
}

void CL_Object::clearErr() {
  setErr(CL_OK,"","");
}

CL_Err_t  CL_Object::errNum() {
  return m_err_num;
}
std::string CL_Object::errMsg() {
  return m_err_msg;
}
std::string CL_Object::errPath() {
  return m_err_path;
}
void CL_Object::printErrMsg() {
  if (m_err_num!=CL_OK) {
    printf("CalvinLite file '%s' had an error: %d: '%s'\n",
           m_err_path.c_str(),
           m_err_num,
           m_err_msg.c_str());
  }
}

//////////

void CL_Object::setParentFile(CL_File* cl_file)
{
  m_parent_file=cl_file;
}

CL_File* CL_Object::getParentFile()
{
  return m_parent_file;
}

//////////

void CL_Object::setDataDirty()
{
  m_data_dirty=true;
}
void CL_Object::setDataDirty(bool val)
{
  m_data_dirty=val;
}
bool CL_Object::isDataDirty()
{
  return m_data_dirty;
}

//////////

int CL_Object::getBytesize() {
  return m_bytesize;
}
bool CL_Object::isBytesizeDirty() {
  return m_bytesize_dirty;
}
int CL_Object::setBytesize(int val) {
  m_bytesize=val;
  m_bytesize_dirty=false;
  return val;
}

//////////

void CL_Object::clearFpos()
{
  m_fpos_start=0;
  m_fpos_size=0;
  m_fpos_end=0;
}

// start is only changed by calling "getFposStart"
// when "Size" or "End" is updated, the other is changed to match.
int CL_Object::getFposStart() const
{
  return m_fpos_start;
}
int CL_Object::setFposStart(int val)
{
  // display updates to m_fpos_start when debugging.
  // if (m_fpos_start!=0) {
  //   printf("setFposStart(%08x)\n",val);
  // }

  //assert(m_fpos_start==0); // only set it once.
  //assert((m_fpos_start==0)||(val==m_fpos_start));
  m_fpos_start=val;
  // scoot the end down
  m_fpos_end=m_fpos_start+m_fpos_size;
  //
  return m_fpos_size;
}
int CL_Object::getFposSize() const
{
  return m_fpos_size;
}
int CL_Object::setFposSize(int val)
{
  assert(0<=val);
  m_fpos_size=val;
  // adjust the end to match
  m_fpos_end=m_fpos_start+m_fpos_size;
  //
  return m_fpos_size;
}
int CL_Object::getFposEnd() const
{
  return m_fpos_end;
}
int CL_Object::setFposEnd(int val)
{
  assert(m_fpos_start<=val);
  m_fpos_end=val;
  m_fpos_size=m_fpos_end-m_fpos_start;
  return val;
}

/////

void CL_Object::clearOrigFpos()
{
  m_orig_fpos_start=0;
  m_orig_fpos_size=0;
  m_orig_fpos_end=0;
}

int CL_Object::getOrigFposStart() const
{
  return m_orig_fpos_start;
}
int CL_Object::setOrigFposStart(int val)
{
  m_orig_fpos_start=val;
  m_orig_fpos_end=m_orig_fpos_start+m_orig_fpos_size;
  return m_orig_fpos_size;
}
int CL_Object::getOrigFposSize() const
{
  return m_orig_fpos_size;
}
int CL_Object::setOrigFposSize(int val)
{
  assert(0<=val);
  m_orig_fpos_size=val;
  m_orig_fpos_end=m_orig_fpos_start+m_orig_fpos_size;
  return m_orig_fpos_size;
}
int CL_Object::getOrigFposEnd() const
{
  return m_orig_fpos_end;
}
int CL_Object::setOrigFposEnd(int val)
{
  assert(m_orig_fpos_start<=val);
  m_orig_fpos_end=val;
  m_orig_fpos_size=m_orig_fpos_end-m_orig_fpos_start;
  return val;
}

//
bool CL_Object::changedFpos() const
{
  if ((m_fpos_start==m_orig_fpos_start) &&
      (m_fpos_size==m_orig_fpos_size)) {
    return false;
  }
  return true;
}
    
