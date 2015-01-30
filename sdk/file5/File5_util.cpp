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
// affy/sdk/file5/File5_util.cpp ---
//
// $Id: File5_util.cpp,v 1.17 2009-10-27 21:07:57 harley Exp $
//

//
#include "file5/File5_util.h"
//
#include "file5/File5.h"
//
#include <iostream>
#include <stdio.h>

//
std::vector<int> affx::file5_dims(int d0) {
  std::vector<int> d;
  d.push_back(d0);
  return d;
}
std::vector<int> affx::file5_dims(int d1,int d0) {
  std::vector<int> d;
  d.push_back(d0);
  d.push_back(d1);
  return d;
}
std::vector<int> affx::file5_dims(int d2,int d1,int d0) {
  std::vector<int> d;
  d.push_back(d0);
  d.push_back(d1);
  d.push_back(d2);
  return d;
}
std::vector<int> affx::file5_dims(int d3,int d2,int d1,int d0) {
  std::vector<int> d;
  d.push_back(d0);
  d.push_back(d1);
  d.push_back(d2);
  d.push_back(d3);
  return d;
}

void affx::file5_dims_hsize_copy(std::vector<hsize_t>& dst,const std::vector<int>& src)
{
  dst.resize(0);
  //
  for (int i=0;i<src.size();i++) {
    dst.push_back(src[i]);
  }
}

//////////

affx::File5_dtype_t affx::as_file5_dtype(affx::tsv_type_t tsv_dtype)
{
  switch (tsv_dtype) {
  case affx::TSV_TYPE_UNKNOWN:
    return affx::FILE5_DTYPE_UNKNOWN;
    break;
  case affx::TSV_TYPE_ERR:
    return affx::FILE5_DTYPE_ERR;
    break;
  case affx::TSV_TYPE_STRING:
    return affx::FILE5_DTYPE_STRING;
    break;
  case affx::TSV_TYPE_INT:
    return affx::FILE5_DTYPE_INT;
    break;
  case affx::TSV_TYPE_FLOAT:
    return affx::FILE5_DTYPE_FLOAT;
    break;
  case affx::TSV_TYPE_DOUBLE:
    return affx::FILE5_DTYPE_DOUBLE;
    break;
  default:
    FILE5_ABORT("Unhandled case.");
    return affx::FILE5_DTYPE_ERR;
    break;
  }
  //
  FILE5_ABORT("Unhandled case. (default didnt catch it.)");
  return affx::FILE5_DTYPE_ERR;
}

// Cant say this as the type_ids arent the  same.
//  if ((hdf5_dtype==H5T_STD_I32LE) ||
//      (hdf5_dtype==H5T_STD_I32BE)) {
//    return affx::FILE5_DTYPE_INT;
//  }
affx::File5_dtype_t affx::as_file5_dtype(hid_t h5_dtype)
{
  //dump_h5_dtype(h5_dtype);

  H5T_class_t dtype_class=H5Tget_class(h5_dtype);
  int dtype_size=H5Tget_size(h5_dtype);

  //
  if (dtype_class==H5T_INTEGER) {
    if (dtype_size==4) {
      return affx::FILE5_DTYPE_INT;
    }
    if (dtype_size==2) {
      return affx::FILE5_DTYPE_SHORT;
    }
    if (dtype_size==1) {
      return affx::FILE5_DTYPE_CHAR;
    }
  }
  else if (dtype_class==H5T_FLOAT) {
    if (dtype_size==4) {
      return affx::FILE5_DTYPE_FLOAT;
    }
    if (dtype_size==8) {
      return affx::FILE5_DTYPE_DOUBLE;
    }
  }
  else if (dtype_class==H5T_STRING) {
    return affx::FILE5_DTYPE_STRING;
  }

  //
  FILE5_ABORT("affx::as_file5_dtype(): cant map h5_dtype.");
  return affx::FILE5_DTYPE_ERR;
}

/////

hid_t affx::as_hdf5_dtype(affx::tsv_type_t tsv_dtype)
{
  return as_hdf5_dtype(as_file5_dtype(tsv_dtype));
}

hid_t affx::as_hdf5_dtype(affx::File5_dtype_t file5_dtype)
{
  switch (file5_dtype) {
  case affx::FILE5_DTYPE_STRING:
    //assert(0);
    break;
  case affx::FILE5_DTYPE_CHAR:
    return H5T_NATIVE_CHAR;
    break;
  case affx::FILE5_DTYPE_SHORT:
    return H5T_NATIVE_SHORT;
    break;
  case affx::FILE5_DTYPE_INT:
    return H5T_NATIVE_INT;
    break;
  case affx::FILE5_DTYPE_FLOAT:
    return H5T_NATIVE_FLOAT;
    break;
  case affx::FILE5_DTYPE_DOUBLE:
    return H5T_NATIVE_DOUBLE;
    break;
  default:
    FILE5_ABORT("Unhandled type.");
    break;
  }
  //
  FILE5_ABORT("internal error! shoudnt be reached.");
  return H5T_NATIVE_DOUBLE;
}

hid_t affx::as_h5t_native_type(hid_t dtype)
{
  int d_class=H5Tget_class(dtype);
  int d_size=H5Tget_size(dtype);

  // translate to native types.
  if (d_class==H5T_INTEGER) {
    if (d_size==4) {
      return H5T_NATIVE_INT;
    }
    if (d_size==2) {
      return H5T_NATIVE_SHORT;
    }
    if (d_size==1) {
      return H5T_NATIVE_CHAR;
    }
  }
  else if (d_class==H5T_FLOAT) {
    if (d_size==4) {
      return H5T_NATIVE_FLOAT;
    }
    if (d_size==8) {
      return H5T_NATIVE_DOUBLE;
    }
  }
  // cant convert this to H5T_NATIVE_*
  return -1;
}

/////

void affx::dump_h5_dset(hid_t dset_id,std::ostream& ostm)
{
  ostm << "dset_id=" << dset_id << " --";
  if (dset_id>0) {
    //
    hid_t dspace_id=H5Dget_space(dset_id);
    FILE5_CHECKID(dspace_id,"dump_h5_dset");
    affx::dump_h5_dspace(dspace_id,ostm);
    H5Sclose(dspace_id);
    //
    if (0) {
      hid_t dtype_id=H5Dget_type(dset_id);
      affx::dump_h5_dtype(dtype_id,ostm);
      H5Tclose(dtype_id);
    }
  }
  else {
    ostm << "  invalid";
  }
  ostm << "\n";
}

void dump_h5_dims(int ndims,hsize_t* dims,std::ostream& ostm)
{
  ostm <<"(";
  for (int d=ndims-1;d>=0;d--) {
    if (dims[d]==H5S_UNLIMITED) {
      ostm << "*";
    }
    else {
      int x=dims[d];
      ostm << x;
    }
    if (d!=0) {
      ostm << ",";
    }
  }
  ostm << ")";
}

void affx::dump_h5_dspace(hid_t dspace_id,std::ostream& ostm)
{
  ostm << "dspace_id=" << dspace_id << " -- ";
  if (dspace_id>0) {
    if (H5Sis_simple(dspace_id)) {
      int ndims=H5Sget_simple_extent_ndims(dspace_id);
      hsize_t dims[20];
      hsize_t dims_max[20];
      //
      H5Sget_simple_extent_dims(dspace_id,dims,dims_max);
      //
      dump_h5_dims(ndims,dims,ostm);
      ostm << "--";
      dump_h5_dims(ndims,dims_max,ostm);
      //
      int npoints=H5Sget_select_npoints(dspace_id);
      if (npoints>=0) {
        ostm << "points=" << npoints;
      }
      //
      int nblocks;
      switch (H5Sget_select_type(dspace_id)) {
      case H5S_SEL_ALL:
        ostm << "{all}";
        break;
      case H5S_SEL_HYPERSLABS:
        nblocks=H5Sget_select_hyper_nblocks(dspace_id);
        if (nblocks>0) {
          hsize_t blocksize[2];
          char str_buf[100];
          ostm << "{";
          for (int b=0;b<nblocks;b++) {
            H5Sget_select_hyper_blocklist(dspace_id,b,1,blocksize);
            sprintf(str_buf,"%d:%d",int(blocksize[0]),int(blocksize[1]));
            ostm << str_buf;
            if (b<nblocks-1) {
              ostm << ",";
            }
          }
          ostm << "}";
          break;
        default:
          break;
        }
      } // switch
    }
    else {
      ostm << "not simple";
    }
  }
  //printf("\n");
}

void affx::dump_h5_dtype(hid_t dtype_id,std::ostream& ostm)
{
  ostm << "dtype_id=" << dtype_id << " -- ";
  if (dtype_id>0) {
    switch (H5Tget_class(dtype_id)) {
    case H5T_INTEGER:
      ostm << "int";
      break;
    case H5T_FLOAT:
      ostm << "float";
      break;
    case H5T_STRING:
      ostm << "string";
      break;
    case H5T_VLEN:
      ostm << "vlen";
      break;
    default:
      ostm << "???";
      break;
    }
    ostm << ":";
    //
    switch (H5Tget_order(dtype_id)) {
    case H5T_ORDER_LE:
      ostm<<"LE";
      break;
    case H5T_ORDER_BE:
      ostm<<"BE";
      break;
    default:
      ostm<<"???";
      break;
    }
    //
    char buf[100];
    sprintf(buf,"%4dB",int(H5Tget_size(dtype_id)));
    ostm << buf;
  }
  else {
    ostm << "invalid";
  }
}

std::string affx::file5_dtype_to_string(affx::File5_dtype_t dtype)
{
  if (dtype==affx::FILE5_DTYPE_CHAR) {
    return "char";
  }
  if (dtype==affx::FILE5_DTYPE_SHORT) {
    return "short";
  }
  if (dtype==affx::FILE5_DTYPE_INT) {
    return "int";
  }
  if (dtype==affx::FILE5_DTYPE_FLOAT) {
    return "float";
  }
  if (dtype==affx::FILE5_DTYPE_DOUBLE) {
    return "double";
  }
  if (dtype==affx::FILE5_DTYPE_STRING) {
    return "string";
  }
  return "unknown";
}

affx::File5_dtype_t affx::file5_string_to_dtype(const std::string& str)
{
  if (str=="char") {
    return affx::FILE5_DTYPE_CHAR;
  }
  if (str=="short") {
    return affx::FILE5_DTYPE_SHORT;
  }
  if (str=="int") {
    return affx::FILE5_DTYPE_INT;
  }
  if (str=="float") {
    return affx::FILE5_DTYPE_FLOAT;
  }
  if (str=="double") {
    return affx::FILE5_DTYPE_DOUBLE;
  }
  if (str=="string") {
    return affx::FILE5_DTYPE_STRING;
  }
  //
  return affx::FILE5_DTYPE_ERR;
}
