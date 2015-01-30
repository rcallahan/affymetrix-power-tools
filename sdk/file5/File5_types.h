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
// sdk/file5/File5_types.h ---
// 
// $Id: File5_types.h,v 1.7 2009-09-11 21:31:12 csugne Exp $
// 

#ifndef _FILE5_TYPES_H_
#define _FILE5_TYPES_H_

// the types defined by HDF5
#include "../external/hdf5/src/hdf5.h"

// the types file5 (aka A5) will use.
namespace affx {
  class File5_File;
  //
  class File5_Group;
  //
  class File5_Object;
  //
  class File5_Matrix;
  //
  class File5_Tsv;
  class File5_TsvColumn;
  //
  class File5_Vector;
  //class File5_Vector_Port;

  //
  enum File5_dtype_t {
    FILE5_DTYPE_UNKNOWN = 0,
    FILE5_DTYPE_ERR     ,
    FILE5_DTYPE_ANY     ,
    FILE5_DTYPE_STRING  ,
    // this should be "BYTE" or I8,I16,I32
    FILE5_DTYPE_CHAR    ,
    FILE5_DTYPE_SHORT   ,
    FILE5_DTYPE_INT     ,
    FILE5_DTYPE_FLOAT   ,
    FILE5_DTYPE_DOUBLE  ,
  };

  enum File5_flag_t {
    FILE5_REPLACE    = 0x01,
    FILE5_CREATE     = 0x02,
    FILE5_OPEN       = 0x04,
    //FILE5_TRUNCATE   = 0x08,
    FILE5_RO         = 0x10,
    FILE5_RW         = 0x11,
    FILE5_NOERR      = 0x80,
    // OPEN||RO
    FILE5_OPEN_RO    = 0x14,
    FILE5_OPEN_CREATE = FILE5_OPEN | FILE5_CREATE, // Open file, creating if necessary
  };

  enum File5_kind_t {
    FILE5_KIND_UNKNOWN = 0,
    FILE5_KIND_FILE,
    FILE5_KIND_GROUP,
    FILE5_KIND_MATRIX,
    FILE5_KIND_OBJECT,
    FILE5_KIND_TSV,
    FILE5_KIND_TSV_COLUMN,
    FILE5_KIND_VECTOR,
  };

  enum File5_state_t {
    FILE5_STATE_UNKNOWN = 0,
    FILE5_STATE_OPEN,
    FILE5_STATE_CLOSED,
  };

  enum File5_return_t {
    FILE5_OK = 0,
    // end of the file/table
    FILE5_EOF,
    //
    FILE5_ERR,         // a generic error
    FILE5_ERR_BADID,
    FILE5_ERR_NOTFOUND,
    FILE5_ERR_REFCNT,
  };
};

#endif
