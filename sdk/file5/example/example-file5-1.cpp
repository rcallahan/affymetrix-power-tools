////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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

// @file   example-file5-1
// @author Harley Gorrell
// @brief  An example of file5 usage.

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
#include <string>
#include <stdio.h>
#include <vector>
//
#ifndef WIN32
#include <unistd.h>
#endif


// Just create the file.
void example_create_file(const std::string& file_name)
{
  // a pointer to a file5 object.
  // it isnt bound to a file yet.
  affx::File5_File* file5=new affx::File5_File();
  // Attempt to create a file with it.
  //
  try {
    file5->open(file_name,affx::FILE5_CREATE);
  }
  catch (...) {
    //
  }
  // even if not open, we can close it.
  file5->close();
  delete file5;
}

void example_open_file(const std::string& file_name)
{
  // a pointer to a file5 object.
  // it isnt bound to a file yet.
  affx::File5_File* file5=new affx::File5_File();
  // Attempt to create a file with it.
  //
  try {
    file5->open(file_name,affx::FILE5_OPEN);
  }
  catch (...) {
    //
  }
  // even if not open, we can close it.
  file5->close();
  delete file5;
}

//////////

void example_write_tsv(const std::string& file_name,int num_rows)
{
  // a pointer to a file5 object.
  // it isnt bound to a file yet.
  affx::File5_File* file5=new affx::File5_File();

  // This will be a new file
  file5->open(file_name,affx::FILE5_REPLACE);

  // "/tsv" is the same as "/tsv" because the tsv is the root group.
  affx::File5_Tsv* tsv5=file5->openTsv("tsv",affx::FILE5_REPLACE);

  int col_max=10;

  char colnamebuf[100];
  for (int cidx=0;cidx<col_max;cidx++) {
    sprintf(colnamebuf,"col-%02d",cidx);
    tsv5->defineColumn(0,cidx,colnamebuf,affx::FILE5_DTYPE_INT);
  }

  int val=0;
  for (int n=0;n<num_rows;n++) {
    for (int cidx=0;cidx<col_max;cidx++) {
      // set_i => sets the integer val (cause we dont want to fiddle with casting.)
      tsv5->set_i(0,cidx,val++);
    }
    tsv5->writeLevel(0);
  }

  //
  tsv5->close();
  delete tsv5;

  // even if not open, we can close it.
  file5->close();
  delete file5;
}

//////////


void example_read_tsv(const std::string& file_name,int max_print_rows)
{
  // a pointer to a file5 object.
  // it isnt bound to a file yet.
  affx::File5_File* file5=new affx::File5_File();

  // This will be a new file
  file5->open(file_name,affx::FILE5_OPEN);

  // "/tsv" is the same as "/tsv" because the tsv is the root group.
  affx::File5_Tsv* tsv5=file5->openTsv("tsv");

  // how many columns are there?
  int col_max=tsv5->getColumnCount(0);

  // header rows
  std::string colname;
  for (int cidx=0;cidx<col_max;cidx++) {
    if (cidx!=0) {
      printf("\t");
    }
    tsv5->getColumnName(0,cidx,&colname);
    printf("%s",colname.c_str());
  }
  printf("\n");

  int val;
  tsv5->rewind();
  while (tsv5->nextLevel(0)==affx::FILE5_OK) {
    // limit the number of lines printed.
    if ((max_print_rows!=0)&&(tsv5->lineNum()>=max_print_rows)) {
      break;
    }
    // go across the current row and print the values.
    for (int cidx=0;cidx<col_max;cidx++) {
      if (cidx!=0) {
        printf("\t");
      }
      tsv5->get(0,cidx,&val);
      printf("%d",val);
    }
    printf("\n");
  }

  //
  tsv5->close();
  delete tsv5;

  // even if not open, we can close it.
  file5->close();
  delete file5;
}

//////////

void example_read_tsv_vec(const std::string& file_name)
{
  // this is what we want to fill.
  std::vector<int> vector_to_fill;

  // This is just like the above.
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_OPEN);
  affx::File5_Tsv* tsv5=file5->openTsv("tsv");

  // this holds a ptr
  affx::File5_TsvColumn* colptr;

  // colptr will be null, we dont have that many columns.
  colptr=tsv5->getColumnPtr(0,999);
  assert(colptr==NULL);

  // now get the first column for real.
  colptr=tsv5->getColumnPtr(0,0);
  // Or you could use the the name as the second arg.
  // colptr=tsv5->getColumnPtr(0,"col-00");
  assert(colptr!=NULL);

  // Now for some magic. 
  // Remember that File5_TsvColumns are File5_vectors.
  // We get its size and adjust the vector we are going to fill.
  int col_len=colptr->size();
  vector_to_fill.resize(col_len);
  // now ask for it to be filled.
  // read_vector will try and put as many items as will fit.
  // Note that this might throw an exception of the data types dont match.
  // (we cant check that they match at compile time.)
  colptr->read_vector(0,&vector_to_fill);
  printf("read_tsv_vec: read: 1: %d of %d elements\n",
         int(vector_to_fill.size()),col_len);
  fflush(NULL);
  // Check the length is correct. (We know it should be this len.)
  assert(vector_to_fill.size()==col_len);

  // if the read is short, the vector is resized to what was 
  // able to be read -- in this case one element.
  colptr->read_vector(col_len-1,&vector_to_fill);
  printf("read_tsv_vec: read: %d of %d elements\n",
         int(vector_to_fill.size()), col_len);
  assert(vector_to_fill.size()==1);

  // close it, just like above.
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

//////////

// All our options are defined here to keep main more readable.
void define_options(PgOptions* opts)
{
  opts->setUsage("example-file5 -- A file5 test program.");

  opts->defOpt("h","help",PgOpt::BOOL_OPT,
               "Print this message.",
               "false");
  // file
  opts->defOpt("","create-file",PgOpt::STRING_OPT,
               "run the example file creation.",
               "");
  opts->defOpt("","open-file",PgOpt::STRING_OPT,
               "run the example file creation.",
               "");
  //
  opts->defOpt("","write-tsv",PgOpt::STRING_OPT,
               "",
               "");
  opts->defOpt("","read-tsv",PgOpt::STRING_OPT,
               "",
               "");
  opts->defOpt("","read-tsv-vec",PgOpt::STRING_OPT,
               "",
               "");
  opts->defOpt("","num-rows",PgOpt::INT_OPT,
               "",
               "100");
}

int main(int argc,const char* argv[]) {
  PgOptions* opts;

  // this needs to be called to init the File5 library.
  // (In turn it calls the HDF5 init function.)
  affx::File5_open();

  // throw our errors for testing.
  Err::setThrowStatus(true);

  // Create our option parser, define the options and parse.
  opts = new PgOptions;
  define_options(opts);
  opts->parseArgv(argv);

  // Print our help message if necessary.
  if (opts->getBool("help")) {
    opts->usage();
  }
  // figure out what the user wanted...
  else if (opts->get("create-file")!="") {
    example_create_file(opts->get("create-file"));
  }
  else if (opts->get("open-file")!="") {
    example_open_file(opts->get("open-file"));
  }
  else if (opts->get("write-tsv")!="") {
    example_write_tsv(opts->get("write-tsv"),opts->getInt("num-rows"));
  }
  else if (opts->get("read-tsv")!="") {
    example_read_tsv(opts->get("read-tsv"),opts->getInt("num-rows"));
  }
  else if (opts->get("read-tsv-vec")!="") {
    example_read_tsv_vec(opts->get("read-tsv-vec"));
  }
  else {
    opts->usage();
  }

  // delete the opts so valgrind doenst report them.
  delete opts;

  // for debugging with valgrind calling _exit will exit the program without
  // letting HDF5 clean up its memory.  This can help id memory leaks.
  // _exit(0);

  // Programs which use HDF5 should call this before exiting.
  // this cleans up resources it has allocated.
  affx::File5_close();

  return 0;
}
