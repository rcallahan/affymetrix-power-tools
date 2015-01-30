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
// affy/sdk/file5/File5_Vector.h ---
//
// $Id: File5_Vector.h,v 1.26 2009-10-28 19:06:41 csugne Exp $
//


#ifndef _FILE5_VECTOR_H_
#define _FILE5_VECTOR_H_

//
#include "file5/File5_File.h"
#include "file5/File5_types.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

// The default chunk sizes.
#define FILE5_CHUNKSIZE        5000
#define FILE5_CHUNKSIZE_STRING  500
#define FILE5_DEFAULT_COMPRESS    9

//////////

// template <typename T1>
// class affx::File5_Vector_Port : public {
// public:
//   //T1 pop_back();
//   virtual T1& at(size_t idx);
//   virtual void assign(size_t idx,const T1& val);
//   virtual void push_back(const T1& val);
//   //
//   virtual int read_vector(size_t idx,std::vector<T1>* vec);
//   virtual int write_vector(size_t idx,const std::vector<T1>& vec);
// };

class affx::File5_Vector : public affx::File5_Object {
public:
  //
  size_t m_file_end_idx;
  size_t m_vec_end_idx; // 
  size_t m_vec_fill_idx;  // idx of where the NEXT value will be written. (==0 -> empty)

  //
  int m_h5_dtype_size;
  int m_dtype_size;
  
  // buffer
  char* m_buf_ptr;
  //int m_buf_cnt;
  int m_buf_max_cnt;
  size_t m_buf_start_idx;
  size_t m_buf_end_idx;

  //
  int m_opt_autoresize;
  int m_opt_string_size;
  int m_opt_buffer_size;
  int m_opt_resize_on_flush;

  //
  File5_Vector();
  virtual ~File5_Vector();

  //
  static affx::File5_Vector* allocate(affx::File5_dtype_t type);

  //
  static affx::File5_Vector* new_vector(int dtype);

  //
  void setOptAutoResize(int flag);
  void setOptChunkSize(int size);
  void setOptStringSize(int size);
  /// Maximum string length allowed.
  int getOptStringSize();
  //
  void setBufferSize(int size);

  //
  void buffer_free();
  void buffer_alloc(int cnt);
  void buffer_seek(size_t idx);
  char* buffer_idx2ptr(size_t idx);
  char* buffer_pushback2ptr();
  void buffer_zero();
  void buffer_clear();
  void buffer_print();
  //
  void check_sanity();

  //
  int open(const std::string& name,affx::File5_dtype_t dtype,int flags);
  int open(const std::string& name);
  //
  int init();
  int flush();
  int close();

  //
  size_t fillIndex();
  size_t setFillIndex(size_t idx);

  //
  int read_scalar( size_t idx,void* ptr);
  int write_scalar(size_t idx,void* ptr);

  // These are internal functions; Please dont use.
  // they bypass the File5 caches.
  int read_array_io( size_t idx,size_t cnt,void* buf);
  int write_array_io(size_t idx,size_t cnt,const void* buf);

  // these methods can be used if you know what you are doing.
  // (IE: you have checked the types and done the allocations.)
  // they flush the File5 cache before doing the IO.
  // you should know the size of the elements.
  int read_array( size_t idx,size_t cnt,void* buf);
  int write_array(size_t idx,size_t cnt,const void* buf);

  //
  int dump();
  int dump1();

  //
  size_t size();
  bool empty();

  //
  virtual char file5_kind_char();

  //
  size_t resize(size_t size);
  size_t reserve(size_t size);
  size_t reserved();
  //
  size_t resizeFile(size_t size);
  size_t reserveFile(size_t size);
  size_t setFileSize(size_t size);
  //
  size_t dtype_size();
  
  inline affx::File5_dtype_t getDType() { return m_dtype; }
  //
  affx::File5_dtype_t file5_dtype();


  
  //// Accessors
  //int pop_back();
  int push_back(int val);
  int get(size_t idx,int* val);
  int get(size_t idx,char* val);
  int get(size_t idx,float* val);
  int get(size_t idx,double* val);

  //
  //int pop_back_c();
  int get_c(size_t idx,char* val);
  //int& at_i(size_t idx);
  int assign_c(size_t idx,char val);
  int push_back_c(char val);
  //
  int read_vector(size_t idx,std::vector<char>* vec);
  int write_vector(size_t idx,const std::vector<char>* vec);

  //
  //int pop_back_s();
  int get_s(size_t idx,short* val);
  //int& at_i(size_t idx);
  int assign_s(size_t idx,short val);
  int push_back_s(short val);
  //
  int read_vector(size_t idx,std::vector<short>* vec);
  int write_vector(size_t idx,const std::vector<short>* vec);

  // ints
  //int pop_back_i();
  int get_i(size_t idx,int* val);
  //int& at_i(size_t idx);
  int assign_i(size_t idx,int val);
  int push_back_i(int val);
  //
  int read_vector(size_t idx,std::vector<int>* vec);
  int write_vector(size_t idx,const std::vector<int>* vec);

  // floats
  //int pop_back_f();
  int get_f(size_t idx,float* val);
  //float& at_f(size_t idx);
  int assign_f(size_t idx,float val);
  int push_back_f(float val);
  int read_vector(size_t idx,std::vector<float>* vec);
  int write_vector(size_t idx,const std::vector<float>* vec);

  // doubles
  //int pop_back_d();
  int get_d(size_t idx,double* val);
  //double& at_d(size_t idx);
  int assign_d(size_t idx,double val);
  int push_back_d(double val);
  //
  int read_vector(size_t idx,std::vector<double>* vec);
  int write_vector(size_t idx,const std::vector<double>* vec);

  // strings
  // int pop_back_string();
  int get(size_t idx,std::string* val);
  //std::string& at_string(size_t idx);
  int assign_string(size_t idx,const std::string& val);
  int push_back_string(const std::string& val);
  //
  int read_vector(size_t idx,std::vector<std::string>* vec);
  int write_vector(size_t idx,const std::vector<std::string>* vec);

  //
  int getAsStr(size_t idx,std::string* str);
  //
  int getAsStrVector(std::vector<std::string>* strvec);
  int getAsStrVector(std::vector<std::string>* strvec,size_t read_start,size_t read_end);

  //
  void assert_dtype(affx::File5_dtype_t dtype,const std::string& msg);

};

#endif // _FILE5_VECTOR_H_
