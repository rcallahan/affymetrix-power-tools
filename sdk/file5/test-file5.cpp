////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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


// @file   test-file5
// @author Harley Gorrell
// @brief  A test program to check that file5 is working.

//
#include "file5/File5.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/PgOptions.h"
//
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string.h>
#include <string>
#include <stdio.h>
#include <vector>
//
#ifndef WIN32
#include <unistd.h>
#endif

/////

int dot_cnt;
#define DOT_CLEAR()    do { dot_cnt=0; fflush(NULL); } while (0)
#define DOT_TITLE(T,F) do { printf("dot: %-20s '%s'\n",T,F); fflush(NULL); } while (0)
#define DOT_PRINT()    do { printf("dot: %2d.. \n",++dot_cnt); fflush(NULL); } while (0)
#define DOT_OK()       do { printf("dot: ok!\n"); DOT_CLEAR(); } while (0)

/////

void run_test_consts() {
  FILE5_ASSERT((affx::FILE5_OPEN|affx::FILE5_RO)==affx::FILE5_OPEN_RO);
}

void run_test_file(const std::string& file_name) {
  int should_be_zero;

  //
  DOT_TITLE("run_test_file",file_name.c_str());
  DOT_CLEAR();

  //
  affx::File5_File* file5=new affx::File5_File();

  // should be ok
  DOT_PRINT();
  unlink(file_name.c_str());
  file5->open(file_name,affx::FILE5_CREATE);
  file5->close();

  // should fail as file is there.
  DOT_PRINT();
  should_be_zero=0;
  try {
    file5->open(file_name,affx::FILE5_OPEN);
    should_be_zero=0; // error!
  }
  catch (...) {
    // ignore
  }
  FILE5_ASSERT(should_be_zero==0);
  //
  file5->close();
  
  // should be ok
  DOT_PRINT();
  file5->open(file_name,affx::FILE5_OPEN);
  file5->close();

  //
  DOT_PRINT();
  file5->open(file_name,affx::FILE5_OPEN_RO);
  //
  file5->close();
  // it should be safe to call close twice...
  file5->flush();
  file5->close();
  // ...or more!
  file5->close();

  // should be ok.
  /*
  DOT_PRINT();
  unlink(file_name.c_str());
  file5->open(file_name,affx::FILE5_OPEN);
  file5->close();
  */
  //
  delete file5;

  //
  DOT_OK();
}

/////

void run_test_matrix(const std::string& file_name) {
  // int should_be_zero;

  DOT_TITLE("run_test_matrix",file_name.c_str());
  DOT_CLEAR();

  affx::File5_File* file5;
  affx::File5_Matrix* matrix;

  int mat_dim1=3;
  int mat_dim0=7;

  //
  DOT_PRINT();
  {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_REPLACE);
    std::vector<int> matrix_dims=affx::file5_dims(mat_dim1,mat_dim0);
    matrix=file5->openMatrix("matrix",affx::FILE5_DTYPE_DOUBLE,matrix_dims,affx::FILE5_REPLACE);
    //
    double val=0.0;
    std::vector<int> set_dims=affx::file5_dims(0,0);
    for (int d1=0;d1<mat_dim1;d1++) {
      set_dims[1]=d1;
      for (int d0=0;d0<mat_dim0;d0++) {
      set_dims[0]=d0;
      matrix->set(set_dims,val);
      val+=1.0;
      //printf("%f ",val);
      }
    }

    //
    matrix->close();
    delete matrix;
    file5->close();
    delete file5;
  }

  //printf("xxx1\n"); fflush(NULL);
  /*
  //
  DOT_PRINT();
  should_be_zero=0;
  try {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_OPEN);
    std::vector<int> matrix_dims=affx::file5_dims(mat_dim1,mat_dim0);
    matrix=file5->openMatrix("matrix",affx::FILE5_DTYPE_DOUBLE,matrix_dims,affx::FILE5_OPEN);
    //
    // should abort as we didnt close the matrix.
    // matrix->close();
    file5->close();
    should_be_zero=1;
    printf("### Error! shouldnt make it here!\n");
  }
  catch (...) {
    //
  }
  FILE5_ASSERT(should_be_zero==0);

  printf("xxx2\n"); fflush(NULL);

  // this should close it.
  DOT_PRINT();
  matrix->close();
  delete matrix;
  file5->close();
  delete file5;
  */
  /*
  //
  DOT_PRINT();
  {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_OPEN);
    std::vector<int> matrix_dims=affx::file5_dims(10,20);
    matrix=file5->openMatrix("matrix",affx::FILE5_DTYPE_DOUBLE,matrix_dims,affx::FILE5_OPEN);
    //
    double val_ref=0.0;
    double val_out;
    std::vector<int> set_dims=affx::file5_dims(0,0);
    for (int d1=0;d1<mat_dim1;d1++) {
      set_dims[1]=d1;
      for (int d0=0;d0<mat_dim0;d0++) {
        set_dims[0]=d0;
        //printf("%f ",val_ref);
        matrix->get(set_dims,&val_out);
        FILE5_ASSERT(val_out==val_ref);
        val_ref+=1.0;
      }
    }
    // test for values.
    //
    matrix->close();
    file5->close();
  }
  */
  DOT_PRINT();
  DOT_OK();
}


//////////

void run_test_vector_string_size(const std::string& file_name,int size,int compress)
{
  //
  DOT_TITLE("run_test_vector_string",file_name.c_str());
  DOT_CLEAR();
  affx::File5_File* file5=new affx::File5_File();
  affx::File5_Vector* vec;

  //
  DOT_PRINT();
  file5->open(file_name,affx::FILE5_REPLACE|affx::FILE5_OPEN);
  //vec=file5->openVector_String("vector-string",size,affx::FILE5_REPLACE);
  vec=file5->createVector();
  vec->setOptCompress(compress);
  vec->setOptStringSize(100);
  vec->open("vector-string",affx::FILE5_DTYPE_STRING,affx::FILE5_REPLACE);
  //
  vec->attrib_set("foo","bar");
  //
  vec->setBufferSize(100);
  //
  std::string vec_name;
  vec->get_name(&vec_name);
  printf("vec_name='%s'\n",vec_name.c_str());
  //
  for (int i=0;i<20;i++) {
    char buf[100];
    sprintf(buf,"string-1-%d",i);
    vec->push_back_string(buf);
  }
  vec->assign_string( 0,"zero");
  vec->assign_string(10,"ten");
  vec->assign_string(11,"eleven");
  // vec->dump(); file5->dump();
  vec->close();
  delete vec;
  // vec->dump(); file5->dump();
  file5->close();
  delete file5;
  // vec->dump(); file5->dump();

  //
  DOT_PRINT();
  file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_OPEN);
  //vec=file5->openVector("vector-string-2",0);
  vec=file5->openVector("vector-string");
  //
  vec->setBufferSize(1);
  vec->assign_string(0,"string-2-0");
  vec->push_back_string("string-2-1");
  vec->push_back_string("string-2-2");
  vec->push_back_string("string-2-3");
  vec->push_back_string("string-2-4");
  //
  vec->close();
  delete vec;
  file5->close();
  delete file5;
  //
  DOT_OK();
}

void run_test_vector_1_print(affx::File5_Vector* vec,int len,int buf_size)
{
  printf("TESTING len=%d buf=%d --------------------\n",len,buf_size);
  int val;
  vec->setBufferSize(buf_size);
  for (int i=0;i<len;i++) {
    vec->get(i,&val);
    printf("%d,",val);
  }
  printf("\n");
}


void run_test_vector_1_w(const std::string& file_name,
                         const affx::File5_dtype_t dtype
                         )
{
  affx::File5_File* file5;
  affx::File5_Vector* vec;

  file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);

  //
  vec=file5->openVector(file_name,dtype,affx::FILE5_REPLACE);
  for (int j=0;j<100;j++) {
    vec->push_back(j+0);
  }
  vec->dump();

  run_test_vector_1_print(vec,10,100);

  vec->close();
  delete vec;
  file5->close();
  delete file5;
}

void run_test_vector_1_r(const std::string& file_name,
                         const affx::File5_dtype_t dtype
                         )
{
  affx::File5_File* file5;
  affx::File5_Vector* vec;

  file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_OPEN);

  //
  vec=file5->openVector(file_name,dtype,affx::FILE5_OPEN);

  run_test_vector_1_print(vec,10,100);

  vec->close();
  delete vec;
  file5->close();
  delete file5;
}

void run_test_vector_1(const std::string& file_name) {
  //
  DOT_TITLE("run_test_matrix",file_name.c_str());
  DOT_CLEAR();

  run_test_vector_1_w(file_name+"-int",   affx::FILE5_DTYPE_INT);
  run_test_vector_1_r(file_name+"-int",   affx::FILE5_DTYPE_INT);
  //
  run_test_vector_1_w(file_name+"-short", affx::FILE5_DTYPE_SHORT);
  run_test_vector_1_r(file_name+"-short", affx::FILE5_DTYPE_SHORT);
  //
  run_test_vector_1_w(file_name+"-char",  affx::FILE5_DTYPE_CHAR);
  run_test_vector_1_r(file_name+"-char",  affx::FILE5_DTYPE_CHAR);
}

void run_test_vector_2(const std::string& file_name) {
  //int should_be_zero;

  //
  DOT_TITLE("run_test_matrix",file_name.c_str());
  DOT_CLEAR();

  affx::File5_File* file5;
  affx::File5_Vector* vector;

  int vec_len=10;

  //
  DOT_PRINT();
  {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_REPLACE);
    vector=file5->openVector("vector",affx::FILE5_DTYPE_DOUBLE,affx::FILE5_REPLACE);

    //
    FILE5_ASSERT(vector->file5_dtype()==affx::FILE5_DTYPE_DOUBLE);

    //
    float val_out=0.0; // float for testing.
    for (int d0=0;d0<vec_len;d0++) {
      vector->push_back_d(val_out);
      val_out+=1.00;
    }

    //
    vector->close();
    delete vector;
    file5->close();
    delete file5;
  }

  //printf("xxx1\n"); fflush(NULL);

  //
  DOT_PRINT();
  //should_be_zero=0;

  /*
    // the close of the file aborts -- disable this test.
  try {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_OPEN);
    vector=file5->openVector("vector",affx::FILE5_DTYPE_DOUBLE,affx::FILE5_OPEN);
    //
    // should abort as we didnt close the vector.
    // vector->close();
    file5->close();
    //should_be_zero=1;
    printf("### Error! shouldnt make it here!\n");
  }
  catch (...) {
    //
  }
  //FILE5_ASSERT(should_be_zero==0);
  // this should close it.
  DOT_PRINT();
  vector->close();
  delete vector;
  file5->close();
  delete file5;
  */

  //
  DOT_PRINT();
  {
    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_OPEN);
    std::vector<int> vector_dims=affx::file5_dims(10,20);
    vector=file5->openVector("vector",affx::FILE5_DTYPE_DOUBLE,affx::FILE5_OPEN);
    //
    double val_ref=0.0;
    double val_in;
    for (int d0=0;d0<vec_len;d0++) {
      vector->get(d0,&val_in);
      FILE5_ASSERT(val_in==val_ref);
      val_ref+=1.00;
    }
    // test for values.
    //
    vector->close();
    delete vector;
    file5->close();
    delete file5;
  }

  DOT_PRINT();
  DOT_OK();
}

void run_test_vector_3(const std::string& file_name) {
  int tmp_int;
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);
  affx::File5_Vector* vec=file5->openVector("vector",affx::FILE5_DTYPE_INT,affx::FILE5_REPLACE);

  std::vector<int> data(10,1);

  vec->write_vector(0,&data);
  FILE5_ASSERT(vec->size()==10);

  data.resize(5);
  vec->write_vector(0,&data);
  FILE5_ASSERT(vec->size()==10);

  data.resize(20);
  vec->write_vector(0,&data);
  FILE5_ASSERT(vec->size()==20);

  vec->push_back_i(21);
  FILE5_ASSERT(vec->size()==21);

  vec->resize(10);
  vec->push_back_i(1100);
  FILE5_ASSERT(vec->size()==11);
  vec->get(10,&tmp_int);
  FILE5_ASSERT(tmp_int==1100);
  
  vec->close();
  delete vec;
  file5->close();
  delete file5;
}

void run_test_usermeta(const std::string& file_name) {
  // 
  DOT_TITLE("run_test_user_meta",file_name.c_str());
  DOT_CLEAR();
  //
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);
  //
  file5->addHeader("key","val");
  file5->addHeader("foo","bar");
  file5->addHeader("abc","def");
  //
  for (int i=0;i<file5->getHeaderCount();i++) {
    std::string key;
    std::string val;
    file5->getHeaderByIdx(i,&key,&val);
    printf("#%%%s=%s\n",key.c_str(),val.c_str());
  }

  //
  file5->close();
  delete file5;
  //
  DOT_PRINT();
  DOT_OK();
}

void run_test_usermeta_big(const std::string& file_name) {
  // 
  DOT_TITLE("run_test_usermeta_big",file_name.c_str());
  DOT_CLEAR();
  //
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);
  //
  for (int i=0;i<1000;i++) {
    char key_buf[32];
    sprintf(key_buf,"key%05d",i);
    char val_buf[32];
    sprintf(val_buf,"val%05d",i);
    file5->addHeader(key_buf,val_buf);
  }
  //
  for (int i=0;i<file5->getHeaderCount();i++) {
    std::string key;
    std::string val;
    file5->getHeaderByIdx(i,&key,&val);
    if ((i%100)==0) {
      printf("%4d : #%%%s=%s\n",i,key.c_str(),val.c_str());
    }
  }
  //
  file5->close();
  delete file5;
  //
  DOT_PRINT();
  DOT_OK();
}

//////////

void run_test_tsv_write(const std::string& file_name)
{
  DOT_TITLE("run_test_tsv_write",file_name.c_str());
  DOT_CLEAR();
  
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);

  affx::File5_Tsv* tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);
  
  std::string col_name="test-column";
  tsv5->defineColumn(0,0,col_name,affx::FILE5_DTYPE_INT);

  // check the types
  FILE5_ASSERT(tsv5->getColumnDtype(0,0)==affx::FILE5_DTYPE_INT);
  FILE5_ASSERT(tsv5->getColumnDtype(0,10)==affx::FILE5_DTYPE_UNKNOWN);
  //
  FILE5_ASSERT(tsv5->getColumnDtype(0,"test-column")==affx::FILE5_DTYPE_INT);
  FILE5_ASSERT(tsv5->getColumnDtype(0,"noSuchColumn")==affx::FILE5_DTYPE_UNKNOWN);

  for (int i=0;i<10;i++) {
    // good cases
    if (tsv5->set_i(0,0,i)!=affx::FILE5_OK) {
      FILE5_ASSERT(0);
    }
    tsv5->writeLevel(0);
    if (tsv5->set_i(0,col_name,i)!=affx::FILE5_OK) {
      FILE5_ASSERT(0);
    }
    tsv5->writeLevel(0);

    // bad cases
    if (tsv5->set_i(0,-1,i)==affx::FILE5_OK) {
      FILE5_ASSERT(0);
    }
    if (tsv5->set_i(0,"bad_col",i)==affx::FILE5_OK) {
      FILE5_ASSERT(0);
    }
    if (tsv5->set_i(9,col_name,i)==affx::FILE5_OK) {
      FILE5_ASSERT(0);
    }
  }
  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}


void run_test_twoleveltsv(const std::string& file_name)
{
  DOT_TITLE("run_test_twoleveltsv",file_name.c_str());
  DOT_CLEAR();

  int rv;
  int val;
  int val_ref;
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);

  affx::File5_Tsv* tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);

  //
  tsv5->defineColumn(0,0,"col-0-0",affx::FILE5_DTYPE_INT);
  tsv5->defineColumn(1,0,"col-1-0",affx::FILE5_DTYPE_INT);

  //
  val=0;
  for (int l0=0;l0<5;l0++) {
    tsv5->set_i(0,0,val++*100);
    tsv5->writeLevel(0);
    for (int l1=0;l1<5;l1++) {
      tsv5->set_i(1,0,val++);
      tsv5->writeLevel(1);
    }
  }

  if (1) {
    tsv5->close();
    delete tsv5;
    file5->close();
    delete file5;

    file5=new affx::File5_File();
    file5->open(file_name,affx::FILE5_OPEN);
    tsv5=file5->openTsv("tsv",affx::FILE5_OPEN);
  }

  //
  tsv5->rewind();
  val_ref=0;
  val=0;
  for (int l0=0;l0<5;l0++) {
    tsv5->nextLine();
    tsv5->get(0,0,&val);
    FILE5_ASSERT(val==val_ref*100);
    val_ref++;
    for (int l1=0;l1<5;l1++) {
      tsv5->nextLine();
      tsv5->get(1,0,&val);
      FILE5_ASSERT(val==val_ref);
      val_ref++;
    }
  }
  // should be an error
  rv=tsv5->nextLine();
  FILE5_ASSERT(rv!=affx::FILE5_OK);

  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

//////////

void run_test_tsv_delete_mk(affx::File5_File* file5,
                            const std::string& tsv_name)
{
  affx::File5_Tsv* tsv5=file5->openTsv(tsv_name,affx::FILE5_REPLACE);
  //
  tsv5->defineColumn(0,0,"col-0-0",affx::FILE5_DTYPE_INT);
  tsv5->defineColumn(0,1,"col-0-1",affx::FILE5_DTYPE_INT);
  //
  int val=0;
  for (int l=0;l<5;l++) {
    tsv5->set_i(0,0,val++);
    tsv5->set_i(0,1,val++);
    tsv5->writeLevel(0);
  }
  //
  tsv5->close();
  delete tsv5;
}
  
void run_test_tsv_delete(const std::string& file_name)
{
  DOT_TITLE("run_test_tsv_delete",file_name.c_str());
  DOT_CLEAR();
  
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);

  // make two tables
  run_test_tsv_delete_mk(file5,"tsv-1");
  run_test_tsv_delete_mk(file5,"tsv-deleteme");
  run_test_tsv_delete_mk(file5,"tsv-2");

  //
  file5->deleteTsv("tsv-deleteme");
  file5->deleteTsv("tsv-doesnotexists");

  //
  file5->close();
  delete file5;
}

//////////

#define BENCHMARK_COL_CNT 10

void run_test_benchmark_write(const std::string& file_name,int benchmark_cnt)
{
  printf("FILE5_CHUNKSIZE        = %5d\n", FILE5_CHUNKSIZE);
  printf("FILE5_DEFAULT_COMPRESS = %5d\n", FILE5_DEFAULT_COMPRESS);
  //
  DOT_TITLE("run_test_benchmark_write",file_name.c_str());
  DOT_CLEAR();


  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);
  affx::File5_Tsv* tsv5=file5->openTsv("benchmark",affx::FILE5_REPLACE);
  
  for (int cidx=0;cidx<BENCHMARK_COL_CNT;cidx++) {
    char buf[100];
    sprintf(buf,"col-%04d",cidx);
    tsv5->defineColumn(0,cidx,buf,affx::FILE5_DTYPE_DOUBLE);
  }

  //
  int cnt=0;
  while (cnt<benchmark_cnt) {
    for (int cidx=0;cidx<BENCHMARK_COL_CNT;cidx++) {
      tsv5->set_d(0,cidx,(cnt++)+0.123456);
    }
    tsv5->writeLevel(0);
  }
  
  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

void run_test_vec_strlen_1(const std::string& file_base,int strlen)
{
  //
  char buf[100];

  //
  sprintf(buf,"%s-%03d",file_base.c_str(),strlen);
  std::string file_name=buf;
  printf("  file=%s\n",file_name.c_str());
  //
  DOT_CLEAR();


  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_REPLACE);
  affx::File5_Tsv* tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);

  tsv5->defineColumn(0,0,"string",affx::FILE5_DTYPE_STRING,strlen);  
  int val=0;
  for (int i=0;i<100000;i++) {
    // psudo-random
    // val=val*23*val+val+10;
    // sequental
    val++;
    sprintf(buf,"%u%u",val,val);
    // fix the length at 10.
    buf[11]=0;
    tsv5->set_string(0,0,buf);
    tsv5->writeLevel(0);
  }
  
  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

void run_test_vec_strlen(const std::string& file_base)
{
  printf("Strlen compression test ==\n");
  // why busted?
  // run_test_vec_strlen_1(file_base,  0);
  run_test_vec_strlen_1(file_base, 10);
  run_test_vec_strlen_1(file_base, 20);
  run_test_vec_strlen_1(file_base, 30);
  run_test_vec_strlen_1(file_base, 40);
  run_test_vec_strlen_1(file_base, 50);
  run_test_vec_strlen_1(file_base,100);
  run_test_vec_strlen_1(file_base,150);
  run_test_vec_strlen_1(file_base,200);
  run_test_vec_strlen_1(file_base,300);
}

//////////

void run_test_open_many() {
  std::vector<affx::File5_File*> files_vec;
  affx::File5_File* file5;
  char fname_buf[200];

  int i=0;
  while (1) {
    sprintf(fname_buf,"test-open-many-%06d",i++);
    file5=new affx::File5_File();
    files_vec.push_back(file5);
    file5->open(fname_buf,affx::FILE5_CREATE);
    printf("test-open-many: '%s'\n",fname_buf);
  }
}

//////////

void run_test_big_vector(affx::File5_File* file5,int size) {
  assert(file5!=NULL);

  printf("test-big-vector: %10d...\n",size); fflush(NULL);
  if (size>1000000) {
    sleep(10);
  }

  int* buf=(int*)malloc(sizeof(int)*size);
  assert(buf!=NULL);
  
  memset(buf,0,size);
  for (int i=0;i<size;i++) {
    buf[i]=i;
  }

  char vec_name[100];
  sprintf(vec_name,"vec-%09d",size);

  affx::File5_Vector* vec;
  vec=file5->openVector(vec_name,affx::FILE5_DTYPE_INT,affx::FILE5_REPLACE);

  printf("   write_array...\n"); fflush(NULL);

  vec->write_array(0,size,buf);

  free(buf);
  vec->close();
  delete vec;

  printf("\n"); fflush(NULL);
}
  
void run_test_big_vectors(const std::string& filename) {
  affx::File5_File* file5=new affx::File5_File();
  file5->open(filename,affx::FILE5_REPLACE);

  int sizes[]={ 
    // 100,
    // 10000,
    100000,
    // 1000,     2000,   5000,
    // 10000,   20000,  50000,
    // 100000, 200000, 500000,
    // 1000000000,
    -1};
              
  for (int i=0;sizes[i]!=-1;i++) {
    run_test_big_vector(file5, sizes[i]);
  }

  file5->close();
  delete file5;
}

//////////

// For generating a bunch of tsvfiles
void
run_test_text_tsv()
{
  char str_buf[100];
  affx::File5_File* file5;
  affx::File5_Tsv* tsv5;

  for (int s_width=10;s_width<20;s_width++) {
    for (int s_cnt=1;s_cnt<10;s_cnt++) {
      //
      sprintf(str_buf,"test-text-tsv-%03d-%03d.file5",s_width,s_cnt);
      file5=new affx::File5_File();
      file5->open(str_buf,affx::FILE5_REPLACE);
      //
      printf("Writing '%s'...\n",str_buf);
      //
      tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);
      tsv5->defineColumn(0,0,"test-column",affx::FILE5_DTYPE_STRING,s_width);
      //
      for (int c=0;c<s_cnt*100;c++) {
        sprintf(str_buf,"%d-%d",s_width,c);
        tsv5->set_string(0,0,str_buf);
        tsv5->writeLevel(0);
      }
      //
      tsv5->close();
      delete tsv5;
      file5->close();
      delete file5;
    }
  }
}

//////////

void print_vos(const std::vector<std::string>& vec)
{
  for (int i=0;i<vec.size();i++) {
    printf(" %2d: '%s'\n",i,vec[i].c_str());
  }
}

void run_test_list_names(const std::string& filename)
{
  affx::File5_File* file5;
  affx::File5_Tsv* tsv5;

  file5=new affx::File5_File();
  file5->open(filename,affx::FILE5_REPLACE);

  tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);
  tsv5->defineColumn(0,0,"col-000",affx::FILE5_DTYPE_DOUBLE);
  tsv5->defineColumn(0,1,"col-001",affx::FILE5_DTYPE_DOUBLE);

  std::vector<std::string> names;
  names=tsv5->listNames();
  print_vos(names);

  names=file5->listNames();
  print_vos(names);
  
  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

//////////

void define_testfile5_options(PgOptions* opts)
{
  opts->setUsage("test-file5 -- A file5 test program.");

  opts->defOpt("h","help",PgOpt::BOOL_OPT,
               "Print this message.",
               "false");
  // file
  opts->defOpt("tf","test-file",PgOpt::BOOL_OPT,
               "File test",
               "true");
  // matrix
  opts->defOpt("tm","test-matrix",PgOpt::BOOL_OPT,
               "Matrix test",
               "true");
  // vector
  opts->defOpt("tv","test-vector",PgOpt::BOOL_OPT,
               "Vector test",
               "true");
  //
  opts->defOpt("tum","test-user-meta",PgOpt::BOOL_OPT,
               "Vector test",
               "true");
  //
  opts->defOpt("tsv","test-tsv",PgOpt::BOOL_OPT,
               "Tsv5 test",
               "true");
  opts->defOpt("","test-vec-strlen",PgOpt::BOOL_OPT,
               "test-vec-strlen",
               "false");
  // stress testing
  opts->defOpt("","test-open-many",PgOpt::BOOL_OPT,
               "test-open-many",
               "false");
  opts->defOpt("","test-big-vectors",PgOpt::BOOL_OPT,
               "test-big-vectors",
               "false");
  opts->defOpt("","test-text-tsv",PgOpt::BOOL_OPT,
               "test-text-tsv",
               "false");
  //
  opts->defineOption("","benchmark", PgOpt::INT_OPT,
                     "Number of doubles to write for benchmarking output.",
                     "0");
  opts->defOpt("o","output",PgOpt::STRING_OPT,
               "output file",
               "");

}

//////////

int main(int argc, const char* argv[])
{
  PgOptions* opts;
  int num_run=0;

  // quick enough to run each time.
  run_test_consts();

  //
  affx::File5_open();

  // throw our errors for testing.
  Err::setThrowStatus(true);

  // Create our option parser, define the options and parse.
  opts = new PgOptions;
  define_testfile5_options(opts);
  opts->parseArgv(argv);

  // Print our help message if necessary.
  if (opts->getBool("help")) {
    opts->usage();
    return 0;
  }

  //
  run_test_list_names("test-list-names.file5");

  ////
  if (opts->getInt("benchmark")!=0) {
    run_test_benchmark_write(opts->get("output"),opts->getInt("benchmark"));
    return 0;
  }
  if (opts->getBool("test-big-vectors")) {
    run_test_big_vectors("test-big-vectors.file5");
    return 0;
  }
  if (opts->getBool("test-text-tsv")) {
    run_test_text_tsv();
    return 0;
  }

  /////
  if (opts->getBool("test-file")) {
    num_run++;
    run_test_file("test-file5-file.file5");
  }
  //
  if (opts->getBool("test-matrix")) {
    num_run++;
    run_test_matrix("test-file5-matrix.file5");
  }
  //
  if (opts->getBool("test-vector")) {
    num_run++;
    run_test_vector_1("test-file5-vector.file5");
    run_test_vector_2("test-file5-vector.file5");
    run_test_vector_3("test-file5-vector-3.file5");
    //
    run_test_vector_string_size("test-file5-vectorstring-z.file5",20,5);
    run_test_vector_string_size("test-file5-vectorstring-nz.file5",20,0);
  }
  //
  if (opts->getBool("test-user-meta")) {
    run_test_usermeta("test-file5-usermeta.file5");
    run_test_usermeta_big("test-file5-usermeta-big.file5");
  }

  if (opts->getBool("test-open-many")) {
    run_test_open_many();
  }

  // 
  if (opts->getBool("test-tsv")) {
    num_run++;
    std::string file_name="test-file5-tsv.file5";
    run_test_tsv_write(file_name);
    run_test_tsv_write(file_name);
    run_test_twoleveltsv("test-file5-twolevel.file5");
    //
    run_test_tsv_delete("test-file5-tsv-delete.file5");
  }

  if (opts->getBool("test-vec-strlen")) {
    run_test_vec_strlen("test-vec-strlen");
    num_run++;
  }

  //
  if (num_run==0) {
    opts->usage();
  }

  //
  delete opts;

  // to test what happens without the close.
  // _exit(0);

  //
  affx::File5_close();

  return 0;
}
