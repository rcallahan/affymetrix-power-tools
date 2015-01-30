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
// CalvinLite/CL_util.cpp ---
//
// $Id: CL_util.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#include "calvinlite/CL_util.h"
#include "portability/affy-system-api.h"
//
#include <assert.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <time.h>

//////////

int CL_gen_uuid_seed=0;

std::string CL_gen_uuid()
{
  // we should use "random()" but that causes problems for Matlab.
  // so here we fiddle with global state -- baaaad.
  // you may want to set CL_gen_uuid_seed to "1" if this has already been done.
  if (CL_gen_uuid_seed==0) {
    CL_gen_uuid_seed=getpid()^(int)time(NULL);
    srand(CL_gen_uuid_seed);
  }

  // generate the random bits of the guid...
  int data1=rand();
  int data2=(rand()%0xFFFF);
  int data3=rand();
  // mark this guid as version 4 (a pseudo-random guid)
  // the version is in the high byte.
  data3=data3&0x0FFF;
  data3=data3|0x4000;
  // data4 is 8 bytes. We generate it in three groups.
  int data41=(rand()&0xFFFF);   // 2B
  int data42=(rand()&0xFFFFFF); // 3B
  int data43=(rand()&0xFFFFFF); // 3B

  // f81d4fae-7dec-11d0-a765-00a0c91e6bf6
  char buf[100];
  snprintf(buf,sizeof(buf),"%08x-%04x-%04x-%04x-%06x%06x",
           data1,data2,data3,data41,data42,data43);
  
  std::string uuid=buf;
  return uuid;
}


std::string CL_gen_timestr() {
	time_t now;
	time(&now);
  return CL_gen_timestr(now);
}

std::string CL_gen_timestr(time_t when) {
	char buf[256];
	struct tm* timeinfo;

	timeinfo = gmtime (&when);
  // approved format seems to be 2009-10-05T16:24:16Z
	int rv = strftime(buf, sizeof(buf), "%Y-%m-%dT%XZ", (timeinfo));
	if (rv<=0){
    assert(0);
	}
  std::string tmp_str=buf;
  return tmp_str;
}

//////////

// Keep in sync with
const char* CL_error_strings[]=
  {
    "No error. (OK).",   // CL_OK
    //
    "Generic Calvin Error.", // CL_ERR
    //
    "File not found.",       // CL_ERR_NOTFOUND
    "File not opened.",      // CL_ERR_NOTOPEN
    "Bad magic number.",     // CL_ERR_HDRMAGIC
    //
    NULL, // __LAST
  };

std::string CL_get_err_string(CL_Err_t err)
{
  if ((err<CL_OK)&&(err>=CL_ERR__LAST)) {
    assert(0);
    return "";
  }
  
  return CL_error_strings[err];
}

//////////

CL_TypeCode_t
CL_mimeStrToCode(const std::string& mimeStr)
{
  if (mimeStr=="text/x-calvin-integer-8") {
    return CL_TC_BYTE;
  }
  if (mimeStr=="text/x-calvin-unsigned-integer-8") {
    return CL_TC_UBYTE;
  }
  if (mimeStr=="text/x-calvin-integer-16") {
    return CL_TC_SHORT;
  }
  if (mimeStr=="text/x-calvin-unsigned-integer-16") {
    return CL_TC_USHORT;
  }
  if (mimeStr=="text/x-calvin-integer-32") {
    return CL_TC_INT;
  }
  if (mimeStr=="text/x-calvin-unsigned-integer-32") {
    return CL_TC_UINT;
  }
  if (mimeStr=="text/x-calvin-float") {
    return CL_TC_FLOAT;
  }
  if (mimeStr=="text/ascii") {
    return CL_TC_TEXT_ASCII16;
  }
  if (mimeStr=="text/plain") {
    return CL_TC_TEXT_PLAIN16;
  }
  //
  printf("CL_mimeStrToCode(): cant convert '%s'.\n",mimeStr.c_str());
  assert(0);
  return CL_TC_BAD;
}

std::string
CL_codeToMimeStr(CL_TypeCode_t tcode)
{
  switch (tcode) {
  case CL_TC_BYTE:
    return std::string("text/x-calvin-integer-8");
    break;
  case CL_TC_UBYTE:
    return std::string("text/x-calvin-unsigned-integer-8");
    break;
  case CL_TC_SHORT:
    return std::string("text/x-calvin-integer-16");
    break;
  case CL_TC_USHORT:
    return std::string("text/x-calvin-unsigned-integer-16");
    break;
  case CL_TC_INT:
    return std::string("text/x-calvin-integer-32");
    break;
  case CL_TC_UINT:
    return std::string("text/x-calvin-unsigned-integer-32");
    break;
  case CL_TC_FLOAT:
    return std::string("text/x-calvin-float");
    break;
  case CL_TC_TEXT_ASCII8:
  case CL_TC_TEXT_ASCII16:
    return std::string("text/ascii");
    break;
  case CL_TC_TEXT_PLAIN8:
  case CL_TC_TEXT_PLAIN16:
    return std::string("text/plain");
    break;
  default:
    printf("CL_codeToMimeStr(): cant convert tcode '%d' to string.\n",tcode);
    assert(0);
    break;
  }
  assert(0);
  return std::string("");
}

int
CL_codeToBytelen(CL_TypeCode_t tcode)
{
  switch (tcode) {
  case CL_TC_BYTE:
  case CL_TC_UBYTE:
    return 1;
    break;
  case CL_TC_SHORT:
  case CL_TC_USHORT:
    return 2;
    break;
  case CL_TC_INT:
  case CL_TC_UINT:
  case CL_TC_FLOAT:
    return 4;
    break;
    // This really shouldnt be hit as the length should be computed elsewhere.
  case CL_TC_TEXT_ASCII16:
  case CL_TC_TEXT_PLAIN16:
    return -1;
    break;
  default:
    assert(0);
    break;
  }
  assert(0);
  return -1;
}

//
int CL_mem_read_i1B(void* ptr) {
  char* p=(char*)ptr;
  char val=*p;
  return val;
}
int CL_mem_read_u1B(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  unsigned char val=*p;
  return val;
}
int CL_mem_read_i2B(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  short val=(p[0]<<8)|(p[1]);
  return val;
}
int CL_mem_read_u2B(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  unsigned short val=(p[0]<<8)|(p[1]);
  return val;
}
int CL_mem_read_i4B(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  int val=(p[0]<<24)|(p[1]<<16)|(p[2]<<8)|(p[3]);
  return val;
}
int CL_mem_read_u4B(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  int val=(p[0]<<24)|(p[1]<<16)|(p[2]<<8)|(p[3]);
  return val;
}

float CL_mem_read_float(void* ptr) {
  unsigned char* p=(unsigned char*)ptr;
  type_punned pun;
  pun.v_int32=(p[0]<<24)|(p[1]<<16)|(p[2]<<8)|(p[3]);
  return pun.v_float;
}

//
void CL_mem_write_i1B(void* ptr,int val) {
  unsigned char* p=(unsigned char*)ptr;
  *p=(val)&0xFF;
}
void CL_mem_write_i2B(void* ptr,int val) {
  unsigned char* p=(unsigned char*)ptr;
  p[0]=(val>> 8)&0xFF;
  p[1]=(val    )&0xFF;
}
void CL_mem_write_i4B(void* ptr,int val) {
  unsigned char* p=(unsigned char*)ptr;
  p[0]=(val>>24)&0xFF;
  p[1]=(val>>16)&0xFF;
  p[2]=(val>> 8)&0xFF;
  p[3]=(val    )&0xFF;
}
void CL_mem_write_float(void* ptr,float val) {
  unsigned char* p=(unsigned char*)ptr;
  type_punned pun;
  pun.v_float=val;
  p[0]=(pun.v_int32>>24)&0xFF;
  p[1]=(pun.v_int32>>16)&0xFF;
  p[2]=(pun.v_int32>> 8)&0xFF;
  p[3]=(pun.v_int32    )&0xFF;
}


//////////

std::string CL_basename(const std::string& path)
{
  std::string::size_type pos=path.rfind("/");
  if (pos==std::string::npos) {
    return path;
  }
  return path.substr(pos+1);
}
