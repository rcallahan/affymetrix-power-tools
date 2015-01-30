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
// affy/sdk/bboard/BboardTypes.cpp ---
//
// $Id: BboardTypes.cpp,v 1.2 2009-11-05 20:41:39 harley Exp $
//

#include "BboardTypes.h"

std::string BboardTypeToString(BboardType_t bbt)
{
  switch (bbt) {
  case BBT_CHAR:
    return "char";
    break;
  case BBT_INT:
    return "int";
    break;
  case BBT_FLOAT:
    return "float";
    break;
  case BBT_DOUBLE:
    return "double";
    break;
  case BBT_STRING:
    return "string";
    break;

//  // vectors of
//  BBT_VEC_CHAR,
//  BBT_VEC_INT,
//  BBT_VEC_FLOAT,
//  BBT_VEC_DOUBLE,
//  BBT_VEC_STRING,
  case BBT_VEC_INT:
    return "vec-int";
    break;
//  // vectors of vectors.
//  BBT_VECVEC_CHAR,
//  BBT_VECVEC_INT,
//  BBT_VECVEC_FLOAT,
//  BBT_VECVEC_DOUBLE,
//  BBT_VECVEC_STRING,
  //
  case BBT_BBOARD:
    return "bboard";
    break;
  case BBT_BBOARDOBJ:
    return "bboardobj";
    break;
  case BBT_DAO_FILE:
    return "dao-file";
    break;
  // an example species object
  case BBT_BIOSPECIES:
    return "biospecies";
    break;
  case BBT_CHIPLAYOUT:
    return "chiplayout";
    break;
  case BBT_PGOPTIONS:
    return "pgoptions";
    break;
    //
    // BBT__NEWTYPE
    //
  default:
    return "*unknown*";
  };
  return "";
}

BboardType_t BboardTypeFromString(const std::string& bbt_str)
{
  if (bbt_str=="char"  ) { return BBT_CHAR  ; }
  if (bbt_str=="int"   ) { return BBT_INT   ; }
  if (bbt_str=="float" ) { return BBT_FLOAT ; }
  if (bbt_str=="double") { return BBT_DOUBLE; }
  if (bbt_str=="string") { return BBT_STRING; }
  //
  if (bbt_str=="vec-double") { return BBT_VEC_DOUBLE; }
  if (bbt_str=="vec-float") { return BBT_VEC_FLOAT; }
  if (bbt_str=="vec-int") { return BBT_VEC_INT; }
  //
  if (bbt_str=="bboard") { return BBT_BBOARD; }
  if (bbt_str=="bboardobj") { return BBT_BBOARDOBJ; }
  //
  if (bbt_str=="chiplayout") { return BBT_CHIPLAYOUT; }
  if (bbt_str=="dao-file") { return BBT_DAO_FILE; }
  if (bbt_str=="pgoptions") { return BBT_PGOPTIONS; }
  //   if (bbt_str=="") { return BBT_; }
  // BBT__NEWTYPE

  return BBT_ERR;
}
