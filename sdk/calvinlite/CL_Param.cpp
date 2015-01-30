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
// ~/CL/trunk/CL_Param.cpp ---
//
// $Id: CL_Param.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//


#ifdef _MSC_VER
// suppress sprintf warning.
#define _CRT_SECURE_NO_WARNINGS
#endif

//
#include "calvinlite/CL_Param.h"
//
#include "calvinlite/CL_util.h"
#include "util/Convert.h"
//
#include <stdio.h>

CL_Param::CL_Param()
{
  clear();
}
CL_Param::CL_Param(const std::string& key,const std::string& val)
{
  clear();
  m_name.setNstr(key);
  m_type_code=CL_TC_TEXT_PLAIN16;
  m_val_string.setNstr(val);
}
CL_Param::CL_Param(const std::string& key,int val)
{
  clear();
  m_name.setNstr(key);
  m_type_code=CL_TC_INT;
  m_val_int=val;
}

CL_Param::CL_Param(const std::string& key,double val)
{
  clear();
  m_name.setNstr(key);
  m_type_code=CL_TC_FLOAT;
  m_val_double=val;
}

CL_Param::CL_Param(const std::string& key,
                   const std::string& type,
                   const std::string& value,
                   int byte_len)
{
  clear();
  m_name.setNstr(key);
  m_type_code=CL_mimeStrToCode(type);
  setValueFromString(value);
  m_suggested_byte_len=byte_len;
}

//////////

void CL_Param::clear()
{
  m_start_fpos=0;
  m_name.clear();
  //
  m_type_code=CL_TC_UNSET;
  //
  m_val_int=0;
  m_val_double=0.0;
  m_val_string.clear();
  //
  m_suggested_byte_len=0;
}

std::string CL_Param::valueAsString()
{
  char buf[100];
  std::string val;

  switch (m_type_code) {
   case CL_TC_BYTE:
   case CL_TC_UBYTE:
   case CL_TC_SHORT:
   case CL_TC_USHORT:
   case CL_TC_INT:
   case CL_TC_UINT:
    sprintf(buf,"%d",m_val_int);
    val=buf;
    break;
  case CL_TC_FLOAT:
    sprintf(buf,"%f",m_val_double);
    val=buf;
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    // just turn it into a narrow (8b) string.
    val=m_val_string.Nstr();
    break;
  default:
    val="";
    break;
  }
  //
  return val;
}

CL_Err_t CL_Param::setValueFromString(const std::string& str)
{
  switch (m_type_code) {
  case CL_TC_BYTE:
  case CL_TC_UBYTE:
  case CL_TC_SHORT:
  case CL_TC_USHORT:
  case CL_TC_INT:
  case CL_TC_UINT:
    m_val_int=Convert::toInt(str);
    break;
  case CL_TC_FLOAT:
    m_val_double=Convert::toDouble(str);
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    // just turn it into a narrow (8b) string.
    m_val_string.setNstr(str);
    break;
  default:
    return CL_ERR;
    break;
  }
  //
  return CL_OK;
}

// @todo should this be "computeByteLen"? "bytelenfortype"?
//       and have another method to return the actual bytelen?
int CL_Param::valueByteLen()
{
  int byte_len=-222;

  switch (m_type_code) {
    // all these are 8 Bytes
  case CL_TC_BYTE   :
  case CL_TC_UBYTE  :
  case CL_TC_SHORT  :
  case CL_TC_USHORT :
  case CL_TC_INT    :
  case CL_TC_UINT   :
  case CL_TC_FLOAT  :
    byte_len=8;
    break;
    //
  case CL_TC_TEXT_ASCII16:
    byte_len=m_val_string.byteLen();
    break;
  case CL_TC_TEXT_PLAIN16:
    byte_len=m_val_string.byteLen();
    break;
    // whoa! bad length.
  default:
    byte_len=-1;
    break;
  }
  
  return byte_len;
}
int CL_Param::suggestedByteLen() {
  int byte_len=valueByteLen();
  // round up to the byte size.
  if (byte_len<m_suggested_byte_len) {
    byte_len=m_suggested_byte_len;
  }
  return byte_len;
}

void CL_Param::dump(const std::string& prefix_in)
{
  std::string prefix=prefix_in;

  CL_DUMP_PREFIX("P   ");
  CL_DUMP_MEMB_STR(m_name);
  CL_DUMP_MEMB_INT(m_type_code);
  //CL_DUMP_MEMB_STR(m_type_mime);
  //printf("m_val_string.m_nstr.size()=%d\n",m_val_string.m_nstr.size());
  CL_DUMP_MEMB_STR(m_val_string);
  CL_DUMP_MEMB_INT(m_val_int);
  CL_DUMP_MEMB_DBL(m_val_double);
}
