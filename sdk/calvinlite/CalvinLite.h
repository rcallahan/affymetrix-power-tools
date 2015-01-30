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
// CalvinLite/CalvinLite.h ---
//
// $Id: CalvinLite.h,v 1.1 2009-10-27 16:53:53 harley Exp $
//

#ifndef _CALVINLITE_H_
#define _CALVINLITE_H_

// @todo:
//  move the byteSizeCompute methods to CL_File
//  raise an error if the option is set and when CL_File::setErr() is called.
//  write_File should be flush_buffer()

// Options for CalvinLite.
#ifndef CL_WITH_TSVFILE
#define CL_WITH_TSVFILE 1
#endif

// fwd decls
class CL_DataColCol;
class CL_DataGroup;
class CL_DataRow;
class CL_DataSet;
class CL_File;
class CL_Gdh;
class CL_Object;
class CL_ObjectWithParams;
class CL_Param;
class CL_string;

// CL_TC_DOUBLE is not defined in the calvin spec; Nor are 64b ints.
enum CL_TypeCode_t {
  CL_TC_BYTE       = 0,
  CL_TC_UBYTE,
  CL_TC_SHORT,
  CL_TC_USHORT,
  CL_TC_INT,
  CL_TC_UINT,
  CL_TC_FLOAT,
  CL_TC_TEXT_ASCII8,
  CL_TC_TEXT_ASCII16,
  CL_TC_TEXT_PLAIN8,
  CL_TC_TEXT_PLAIN16,
  //
  CL_TC_UNSET      = 100,
  CL_TC_BAD        = 101
};

enum CL_ReadOpt_t {
  CL_READOPT_ALL      =0x00,
  CL_READOPT_HEADONLY =0x01,
};

// keep the strings in sync with the CL_ERR codes!
enum CL_Err_t {
  CL_OK  = 0,
  //
  CL_ERR,
  //
  CL_ERR_NOTFOUND,
  CL_ERR_NOTOPEN,
  CL_ERR_HDRMAGIC,
  //
  CL_ERR__LAST // place holder
};

#ifndef _AFFY_TYPE_PUNNED_
#define _AFFY_TYPE_PUNNED_
union type_punned {
  float v_float;
  int v_int32;
  unsigned int v_uint32;
};
#endif

////

#define CL_DUMP_PREFIX(_lbl) { \
    printf("%s%s: %p ----------\n",prefix.c_str(),_lbl,this); \
    prefix=prefix+_lbl; \
  }
#define CL_DUMP_BREAK() { \
    printf("%s\n%s\n",prefix.c_str(),prefix.c_str()); \
  }
#define CL_DUMP_MEMB_INT(_memb) { \
    printf("%s%-30s : %10d   %08x\n",prefix.c_str(),#_memb,_memb,_memb); \
  }
#define CL_DUMP_MEMB_DBL(_memb) { \
    printf("%s%-30s : %f\n",prefix.c_str(),#_memb,_memb); \
  }
#define CL_DUMP_MEMB_STR(_memb) { \
    printf("%s%-30s : '%s'\n",prefix.c_str(),#_memb,_memb.c_str()); \
  }

#endif
