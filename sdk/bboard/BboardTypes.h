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
//  ~/Affy/apt2/trunk/BboardTypes.h ---
// 
//  $Id: BboardTypes.h,v 1.5 2009-11-05 20:41:39 harley Exp $
// 

//  Author:  harley <harley@mahalito.net>
//  Keywords: 

//  Commentary:
//  * 
//  Code:

#ifndef _BBOARDTYPES_H_
#define _BBOARDTYPES_H_

class Bboard;
class BboardBox;
class BboardBoxRef;
class BboardObj;

#include <string>

/// @todo replace doCreate, doParents with these.
/// unused
enum BboardGetFlags_t {
  BB_LOCAL,
  BB_CREATE,
  BB_PARENT
};

/// The codes of what can be stored into a BboardBox.
enum BboardType_t {
  // @todo BBT_
  BBT_UNSET=0,
  BBT_ERR=1,
  BBT_NULLREF,
  //
  BBT__START_OF_TYPES=10,
  // base types
  BBT_CHAR,
  BBT_INT,
  BBT_FLOAT,
  BBT_DOUBLE,
  BBT_STRING,
  //
  BBT_CHAR_PTR,
  // vectors of X
  BBT_VEC_CHAR,
  BBT_VEC_INT,
  BBT_VEC_FLOAT,
  BBT_VEC_DOUBLE,
  BBT_VEC_STRING,
  // vectors of vectors of X
  BBT_VECVEC_CHAR,
  BBT_VECVEC_INT,
  BBT_VECVEC_FLOAT,
  BBT_VECVEC_DOUBLE,
  BBT_VECVEC_STRING,
  // APT2 class types
  BBT_BBOARD,
  BBT_BBOARDOBJ,
  // an example species object
  BBT_BIOSPECIES,
  // chipstream stuff
  BBT_CHIPLAYOUT, 
  BBT_PGOPTIONS,
  //
  BBT_FSPATH,
  BBT_DAO_FILE,
  BBT_DAO_GROUP,
  BBT_DAO_TABLE,
  //
  BBT_EXAMPLETESTOBJ,
  // When adding a new type, search for "BBT__NEWTYPE"
  // and add your new type there.
  // BBT__NEWTYPE
  //
  BBT__END_OF_TYPES
};

/// these should be in "libaffyutil", but they are here for now.

/// Translate codes to strings...
std::string BboardTypeToString(BboardType_t bbt);
/// ...and strings to codes.
BboardType_t BboardTypeFromString(const std::string& bbt_str);

#endif // _BBOARDTYPES_H_
