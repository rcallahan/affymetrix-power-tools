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
// affy/sdk/file5/File5.h ---
// 
// $Id: File5.h,v 1.16 2009-09-18 03:37:27 mspald Exp $
// 

// The "sdk/file5" library, aka "file5" is intented to be a
// simple set of objects for APT to use in reading and
// writing HDF5 format files.

// TODO:
// goals today:
//   * convert to m_parent.
//   * reading of tsv
//   * opening for Matrix/Tsv/Vector

// * inserting / update of rows
// * adding removing columns
// * merging columns
// * File5_Group
// * apt-file5-util finished
// * usage examples

// use "H5Eset_auto(H5E_auto_tfunc, void *client_data )" to supress error printing.

// Class layout
//
// * File5_Object
// |
// +--* File5_File
// |
// +--* File5_Group
// |
// +--* File5_Matrix
// |
// +--* File5_Tsv
// |
// +--* File5_Vector  == resizeable vector
//    |
//    +--* File5_TsvColumn

#ifndef _FILE5_H_
#define _FILE5_H_

// #define FILE5_DEBUG_PRINT 1

/// used
#define TSV5_PREFIX "tsv-"
/// when exporting and importing meta info, use this header prefix.
#define FILE5_TSVMETA_PREFIX "file5-tsv-meta"

// @todo:
//  export all the Tsv parts from a file5 to Tsv at once.

//
#include "file5/File5_File.h"
#include "file5/File5_Group.h"
#include "file5/File5_Matrix.h"
#include "file5/File5_Object.h"
#include "file5/File5_Tsv.h"
#include "file5/File5_Vector.h"
#include "file5/File5_types.h"
#include "file5/File5_util.h"
//

#endif // _FILE5_H_
