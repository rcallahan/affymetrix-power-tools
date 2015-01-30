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

// @file   apt-file5-util
// @author Harley Gorrell
// @brief  A command to Work with APT files in the HDF5 format.

//
#include "file5/File5.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Convert.h"
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


/////

void find_data(const std::string& file_name,
               const std::string& tsv_name,
               const std::string& col_name,
               const std::string& val)
{
  //
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_name,affx::FILE5_OPEN);
  affx::File5_Tsv* tsv5;
  tsv5=file5->openTsv(tsv_name,affx::FILE5_OPEN);

  // print the header
  tsv5->printfHeader();
  if (tsv5->findFirst(0,col_name,val)==affx::FILE5_OK) {
    tsv5->printfRow();
  }

  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;
}

//////////

/// @todo: these functions should be relocated to the library.

void TsvToTsv5_convert(const std::string& file_in,const std::string& file_out,const std::string& tsv_name, bool append=false);
void TsvToTsv5_copy_schema(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out);
void TsvToTsv5_copy_metadata(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out);
void TsvToTsv5_copy_data(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out);

void Tsv5ToTsv_convert(const std::string& file_in,const std::string& tsv5_name,const std::string& file_out);
void Tsv5ToTsv_copy_schema(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out);
void Tsv5ToTsv_copy_metadata(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out);
void Tsv5ToTsv_copy_data(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out);



/////

void TsvToTsv5_convert(const std::string& file_in,
                       const std::string& file_out,
                       const std::string& tsv_name,
                       bool copymetadata,
                       bool append)
{
  printf("### tsv_to_tsv5: '%s' to '%s:%s'...\n",
         file_in.c_str(),
         file_out.c_str(),tsv_name.c_str());
  //
  affx::TsvFile* tsv_in=new affx::TsvFile();
  tsv_in->open(file_in);

  //
  affx::File5_File* file5=new affx::File5_File();
  if(append) {
    file5->open(file_out,affx::FILE5_OPEN|affx::FILE5_CREATE);
  } else {
    file5->open(file_out,affx::FILE5_REPLACE);
  }
  if (!file5->is_open()) {
    Err::errAbort("Cant open '"+file_out+"'");
  }

  //
  affx::File5_Tsv* tsv5_out;
  tsv5_out=file5->openTsv(tsv_name,affx::FILE5_REPLACE);

  //
  TsvToTsv5_copy_schema(tsv_in,tsv5_out);
  if (copymetadata) {
    TsvToTsv5_copy_metadata(tsv_in,tsv5_out);
  }
  TsvToTsv5_copy_data(tsv_in,tsv5_out);

  //
  tsv5_out->close();
  delete tsv5_out;
  // force a core dump after tsv5_out is closed...
  // Any rechable allocations we see should be as a result of it.
  // assert(0);
  file5->close();
  delete file5;
  //
  tsv_in->close();
  delete tsv_in;
}

void TsvToTsv5_copy_schema(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out)
{
  // Figure out what the types of the Tsv5 columns should be.
  //if (tsv_in->deduce_types()!=affx::TSV_OK) {
  //  Err::errAbort("Unable to deduce_types for '"+tsv_in->getFileName()+"'");
  //}
  // what are the maximum lengths of the columns?
  //tsv_in->deduce_sizes();

  // do the deductions first.
  tsv_in->deduce_types();
  tsv_in->deduce_sizes();

  std::string key;
  std::string val_str;
  int val_int;
  //
  affx::TsvFileField* col;
  std::string col_name;
  affx::tsv_type_t col_type;
  affx::File5_dtype_t col_dtype;
  int col_max_size;

  // create the columns
  int clvl_max=tsv_in->getLevelCount();
  for (int clvl=0;clvl<clvl_max;clvl++) {
    int cidx_max=tsv_in->getColumnCount(clvl);
    for (int cidx=0;cidx<cidx_max;cidx++) {
      // grab a ref to this column
      col=tsv_in->clvlcidx2colptr(clvl,cidx);
      APT_ERR_ASSERT(col!=NULL,"internal error.");
      col_name=col->get_name();
      // start with what we guessed.
      col_type=col->get_type();
      // look to see if there is a header which declares the type.
      key=FILE5_TSVMETA_PREFIX"-"+ToStr(clvl)+"-"+ToStr(cidx);
      if (tsv_in->getHeader(key+"-type",val_str)==affx::TSV_OK) {
        col_type=affx::TsvFile::stringToColType(val_str);
        APT_ERR_ASSERT(col_type!=affx::TSV_TYPE_ERR,"unknown TsvFile type.");
      }
      if (tsv_in->getHeader(key+"-len",val_int)==affx::TSV_OK) {
        col->m_max_size=val_int;
      }
      //
      col_dtype=affx::as_file5_dtype(col_type);
      col_max_size=col->get_max_size();
      //
      //printf("col: %d: %d: type=%2d max_size=%d name='%s'\n",clvl,cidx,col_dtype,col_max_size,col_name.c_str());
      //
      tsv5_out->defineColumn(clvl,cidx,col_name,col_dtype,col_max_size);
    }
  }
  // did we define the Tsv5 columns correctly?
  // tsv5_out->dump();
}

void TsvToTsv5_copy_metadata(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out)
{
  std::string key;
  std::string val;

  tsv_in->headersBegin();
  while (tsv_in->headersNext(key,val)==affx::TSV_OK) {
    // not a meta prefix?
    if (key.find(FILE5_TSVMETA_PREFIX)!=0) {
      tsv5_out->addHeader(key,val);
    }
  }
}

void TsvToTsv5_copy_data(affx::TsvFile* tsv_in,affx::File5_Tsv* tsv5_out)
{
  // Now that the columns are defined
  int clvl;
  while (tsv_in->nextLine()==affx::TSV_OK) {
    clvl=tsv_in->lineLevel();
    if (clvl<0) {
      continue;
    }
    //
    //printf(" %3d:",clvl);
    //
    for (int cidx=0;cidx<tsv_in->getColumnCount(clvl);cidx++) {
      std::string val;
      tsv_in->get(clvl,cidx,val);
      //printf(" '%s' ",val.c_str());
      //
      int val_i;
      float val_f;
      double val_d;
      std::string val_string;

      affx::File5_TsvColumn* col_ptr=tsv5_out->getColumnPtr(clvl,cidx);
      assert(col_ptr!=NULL);
      affx::File5_dtype_t dtype=col_ptr->m_dtype;

      switch (dtype) {
      case affx::FILE5_DTYPE_CHAR:
      case affx::FILE5_DTYPE_SHORT:
      case affx::FILE5_DTYPE_INT:
        tsv_in->get(clvl,cidx,val_i);
        tsv5_out->set_i(clvl,cidx,val_i);
        break;
        //
      case affx::FILE5_DTYPE_FLOAT:
        tsv_in->get(clvl,cidx,val_f);
        tsv5_out->set_f(clvl,cidx,val_f);
        break;
      case affx::FILE5_DTYPE_DOUBLE:
        tsv_in->get(clvl,cidx,val_d);
        tsv5_out->set_d(clvl,cidx,val_d);
        break;
        //
      case affx::FILE5_DTYPE_STRING:
        tsv_in->get(clvl,cidx,val_string);
        tsv5_out->set_string(clvl,cidx,val_string);
        break;
      default:
        Err::errAbort("dtype.");
        assert(0);
      }
    }
    //printf("\n");
    tsv5_out->writeLevel(clvl);
  }

  printf("### converted %d rows.\n",tsv_in->lineNum());
}

//////////

void Tsv5ToTsv_convert(const std::string& file_in,
                       const std::string& internal_name,
                       const std::string& file_out,
                       int precision)
{
  printf("### tsv5_to_tsv: '%s' to '%s'...\n",file_in.c_str(),file_out.c_str());

  //
  affx::TsvFile* tsv_out=new affx::TsvFile();
  // make sure that escapes are used.
  tsv_out->m_optEscapeOk=true;
  // and set default precision 
  tsv_out->setPrecision(precision);

  //
  affx::File5_File* file5=new affx::File5_File();
  file5->open(file_in,affx::FILE5_OPEN_RO);
  //
  affx::File5_Tsv* tsv5_in;
  tsv5_in=file5->openTsv(internal_name,affx::FILE5_OPEN);

  // Add a bunch of headers describing the conversion process.
  tsv_out->addHeader(FILE5_TSVMETA_PREFIX,"Headers which begin with this prefix are ignored.");
  tsv_out->addHeader(FILE5_TSVMETA_PREFIX,"They will be stripped when converted into File5_Tsv objects.");
  tsv_out->addHeader(FILE5_TSVMETA_PREFIX"-filename",file_in);
  tsv_out->addHeader(FILE5_TSVMETA_PREFIX"-internalname",internal_name);

  Tsv5ToTsv_copy_schema(tsv5_in,tsv_out);
  Tsv5ToTsv_copy_metadata(tsv5_in,tsv_out);

  tsv_out->writeTsv(file_out);

  Tsv5ToTsv_copy_data(tsv5_in,tsv_out);

  //
  tsv_out->close();
  delete tsv_out;
  //
  tsv5_in->close();
  delete tsv5_in;
  file5->close();
  delete file5;
}

void Tsv5ToTsv_copy_schema(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out)
{
  // create the columns
  int clvl_max=tsv5_in->getLevelCount();
  for (int clvl=0;clvl<clvl_max;clvl++) {
    int cidx_max=tsv5_in->getColumnCount(clvl);
    for (int cidx=0;cidx<cidx_max;cidx++) {
      //
      //affx::File5_TsvColumn* col=tsv5_in->getColumn(clvl,cidx);
      //int col_type=col->get_dtype();
      //affx::tsv_type_t tsv_ctype=tsv_in.get_type(clvl,cidx);
      //affx::File5_dtype_t cdtype=affx::as_file5_dtype(tsv_ctype);
      std::string cname;
      tsv5_in->getColumnName(clvl,cidx,&cname);
      //printf("col: %d: %d: type=%2d name='%s'\n",clvl,cidx,col_type,cname.c_str());
      //
      tsv_out->defineColumn(clvl,cidx,cname);
      //
      affx::File5_TsvColumn* col_ptr=tsv5_in->getColumnPtr(clvl,cidx);
      affx::File5_dtype_t dtype=col_ptr->m_dtype;
      std::string key=FILE5_TSVMETA_PREFIX"-"+ToStr(clvl)+"-"+ToStr(cidx);
      tsv_out->addHeader(key+"-type",affx::file5_dtype_to_string(dtype));
      if (dtype==affx::FILE5_DTYPE_STRING) {
        tsv_out->addHeader(key+"-len",col_ptr->m_dtype_size);
      }
    }
  }
}

void Tsv5ToTsv_copy_metadata(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out)
{
  std::string key;
  std::string val;

  tsv5_in->headersBegin();
  while (tsv5_in->headersNext(&key,&val)==affx::FILE5_OK) {
    // printf ("Tsv5ToTsv: addHeader('%s','%s')\n",key.c_str(),val.c_str());
    tsv_out->addHeader(key,val);
  }
}

void Tsv5ToTsv_copy_data(affx::File5_Tsv* tsv5_in,affx::TsvFile* tsv_out)
{
  tsv5_in->rewind();

  while (tsv5_in->nextLine()==affx::FILE5_OK) {
    int clvl=tsv5_in->lineLevel();
    //
    for (int cidx=0;cidx<tsv5_in->getColumnCount(clvl);cidx++) {
      // copy the correct type of data.
      char val_c;
      int val_i;
      float val_f;
      double val_d;
      std::string val_s;

      affx::File5_TsvColumn* col_ptr=tsv5_in->getColumnPtr(clvl,cidx);
      assert(col_ptr!=NULL);
      affx::File5_dtype_t dtype=col_ptr->m_dtype;

      switch (dtype) {
      case affx::FILE5_DTYPE_CHAR:
        tsv5_in->get(clvl,cidx,&val_c);
        tsv_out->set(clvl,cidx,val_c);
        break;
      case affx::FILE5_DTYPE_SHORT:
        tsv5_in->get(clvl,cidx,&val_s);
        tsv_out->set(clvl,cidx,val_s);
        break;
      case affx::FILE5_DTYPE_INT:
        tsv5_in->get(clvl,cidx,&val_i);
        tsv_out->set(clvl,cidx,val_i);
        break;
      case affx::FILE5_DTYPE_FLOAT:
        tsv5_in->get(clvl,cidx,&val_f);
        tsv_out->set(clvl,cidx,val_f);
        break;
      case affx::FILE5_DTYPE_DOUBLE:
        tsv5_in->get(clvl,cidx,&val_d);
        tsv_out->set(clvl,cidx,val_d);
        break;
      case affx::FILE5_DTYPE_STRING:
        tsv5_in->get(clvl,cidx,&val_s);
        tsv_out->set(clvl,cidx,val_s);
        break;
      default:
        Err::errAbort("Tsv5ToTsv_copy_data: bad type. (dtype="+ToStr(dtype)+")");
        assert(0);
      }
    }
    tsv_out->writeLevel(clvl);
  }
  printf("### converted %d rows\n",tsv5_in->lineNum()+1);
}


void generate_test_file(const std::string& file_name,
                        int num_rows,
                        int num_levels,
                        int num_cols_str,
                        int num_cols_int,
                        int num_cols_double)
{
  printf("Creating test file: '%s' (s=%d,i=%d,d=%d)\n",
         file_name.c_str(),
         num_cols_str,num_cols_int,num_cols_double);
  //
  affx::TsvFile* tsv_out=new affx::TsvFile();
  //
  char col_name[100];
  int cidx;
  //

  for (int clvl=0;clvl<num_levels;clvl++) {
    cidx=0;
    for (int i=0;i<num_cols_str;i++,cidx++) {
      sprintf(col_name,"test-s-%03d-%03d",clvl,cidx);
      tsv_out->defineColumn(clvl,cidx,col_name);
    }
    for (int i=0;i<num_cols_int;i++,cidx++) {
      sprintf(col_name,"test-i-%03d-%03d",clvl,cidx);
      tsv_out->defineColumn(clvl,cidx,col_name);
    }
    for (int i=0;i<num_cols_double;i++,cidx++) {
      sprintf(col_name,"test-d-%03d-%03d",clvl,cidx);
      tsv_out->defineColumn(clvl,cidx,col_name);
    }
  }

  // add some '#%key=val' headers.
  for (int i=0;i<20;i++) {
    char key_buf[100];
    char val_buf[100];
    sprintf(key_buf,"key-%04d",i);
    sprintf(val_buf,"val-%04d",i);
    tsv_out->addHeader(key_buf,val_buf);
  }

  //
  tsv_out->writeTsv(file_name);

  //
  int val=0;
  char val_buf[100];

  for (int i=0;i<num_rows;i++) {
    int cidx=0;
    int clvl=i%num_levels;
    //
    for (int i=0;i<num_cols_str;i++,cidx++) {
      sprintf(val_buf,"str-%06d",val++);
      tsv_out->set(clvl,cidx,val_buf);
    }
    for (int i=0;i<num_cols_int;i++,cidx++) {
      tsv_out->set(clvl,cidx,val++);
    }
    for (int i=0;i<num_cols_double;i++,cidx++) {
      tsv_out->set(clvl,cidx,(val++)+0.01);
    }
    //
    tsv_out->writeLevel(clvl);
  }

  //
  tsv_out->close();
  delete tsv_out;
}

//////////

void define_aptfile5util_options(PgOptions* opts)
{
  opts->setUsage("apt-file5-util -- A utility for working with APT file5 formats.\n"
                 "\n"
                 "EXAMPLES:\n"
                 "\n"
                 "* To convert an A5 TsvReport to text:\n"
                 "    apt-file5-util --to-tsv --internal-name CN5/CN5.calls -o calls.txt  calls.a5\n"
                 "  You can find the 'internal-name' by running 'h5dump -n' on the A5 file.\n"
                 "  The name is the path up to the 'tsv-col-000-000' dataset.\n"
                 "\n"
                 "* To convert one text tsv file to A5:\n"
                 "    apt-file5-util --to-tsv5 --internal-name=CN5/CN5.calls -o calls.a5 calls.txt\n"
                 "\n"
                 "If you have multiple tsv files which have a 'file5-tsv-meta-internalname',\n"
                 "(which is what  --to-tsv will add) then several tsvfiles can be imported\n"
                 "at the same time.\n"
                 "   apt-file5-util --to-tsv5 -o OUTPUT.a5 *.tsv\n"
                 "\n"
                 "   Note: you may want to run h5repack -i INPUT.a5 -o OUTPUT.a5 to reclaim space.\n"
                 "\n"
                 "* To find and print a matching value in an tsv5 dataset:\n"
                 "    apt-file5-util -i test9.tsv5 --internal-name test9 --find-col COLNAME --find-val MATCHVAL\n"
                 "  The column names can be found by converting the file to text format.\n"
                 "\n"
                 "NOTE: '--append' isnt the best name.  It appends the data set to an existing HDF5/A5 file.\n"
                 "      It does the append by REPLACING the dataset in that location.\n"
                 "      If there isnt a dataset by that name, then it really is an append.\n"
                 "\n"
                 "NOTE: After adding or replacing data in an HDF5/A5 file, there might be some extra space\n"
                 "      left unallocated in the file.  To reclaim this space to should run 'h5repack'.\n"
                 "      (This should be done to reference files as they will be copied many times.)\n"
                 "          h5repack -v -i file_from.a5 -o file_to.a5\n"
                 "\n"
                 );

  opts->defOpt("h","help",PgOpt::BOOL_OPT,
               "Print this message.",
               "false");

  opts->defOpt("5","to-tsv5",PgOpt::BOOL_OPT,
               "Convert a tsv file to the tsv5 format.",
               "false");

  opts->defOpt("tsv","to-tsv",PgOpt::BOOL_OPT,
               "Convert a tsv file to the tsv5 format.",
               "false");

  opts->defOpt("","append",PgOpt::BOOL_OPT,
               "Append to an existing a5 file.",
               "false");

  opts->defOpt("","metadata",PgOpt::BOOL_OPT,
               "Copy the metadata ('#%' lines) to the tsv5 file.",
               "true");

  opts->defOpt("o","output",PgOpt::STRING_OPT,
               "The name of the output file.",
               "");

  opts->defOpt("p","precision",PgOpt::INT_OPT,
               "The text output precision..",
               "6");

  opts->defOpt("","internal-name",PgOpt::STRING_OPT,
               "Selects the name of the File5_Tsv "
               "when converting tsv files. IE 'CN5/CN5.calls'.",
               "");
  opts->defOpt("i","file-name",PgOpt::STRING_OPT,
               "The file name to read from.",
               "");

  //
  opts->defOpt("mn","matrix-file",PgOpt::STRING_OPT,
               "Create a matrix in this file.",
               "");
  opts->defOpt("ms","matrix-size",PgOpt::INT_OPT,
               "Create a matrix of this size.",
               "10");
  //
  opts->defOpt("vf","vector-file",PgOpt::STRING_OPT,
               "Create a vector in this file.",
               "");
  opts->defOpt("vn","vector-name",PgOpt::STRING_OPT,
               "Vector name to create or dump",
               "vector1");
  opts->defOpt("dvf","dump-vector-file",PgOpt::STRING_OPT,
               "Dump the vector in this file.",
               "");
  //
  opts->defOpt("gtf","gen-test-file",PgOpt::STRING_OPT,
               "Generate a test file.",
               "");
  opts->defOpt("gtfs","gen-test-file-size",PgOpt::INT_OPT,
               "Size of the test file to generate.",
               "1000");
  opts->defOpt("gtfl","gen-test-file-levels",PgOpt::INT_OPT,
               "Number of levels to generate.  (1 or more.)",
               "1");
  opts->defOpt("gtf-s","gen-test-str-cols",PgOpt::INT_OPT,
               "",
               "10");
  opts->defOpt("gtf-i","gen-test-int-cols",PgOpt::INT_OPT,
               "",
               "10");
  opts->defOpt("gtf-d","gen-test-double-cols",PgOpt::INT_OPT,
               "",
               "10");

  //
  opts->defOpt("","find-val",PgOpt::STRING_OPT,
               "Value to find",
               "");
  opts->defOpt("","find-col",PgOpt::STRING_OPT,
               "Column to search when finding.",
               "");
  //
  opts->defOpt("lc","line-count",PgOpt::BOOL_OPT,
               "Count the lines in a Tsv5 File.",
               "false");
}

int main(int argc, const char* argv[]) {
  try {
    PgOptions* opts;

    // Create our option parser, define the options and parse.
    opts = new PgOptions;
    define_aptfile5util_options(opts);
    opts->parseArgv(argv);

    // Print our help message if necessary.
    if (opts->getBool("help")) {
      opts->usage();
    }

    //
    else if ((opts->get("find-val")!="")||(opts->get("find-col")!="")) {
      if (opts->get("find-val")=="") {
        Err::errAbort("'--find-val' must be set.");
      }
      if (opts->get("find-col")=="") {
        Err::errAbort("'--find-col' must be set.");
      }
      //
      find_data(opts->get("file-name"),opts->get("internal-name"),opts->get("find-col"),opts->get("find-val"));
    }

    //
    else if (opts->get("gen-test-file")!="") {
      std::string file_name=opts->get("gen-test-file");
      int num_rows=opts->getInt("gen-test-file-size");
      int num_levels=opts->getInt("gen-test-file-levels");
      generate_test_file(file_name,
                         num_rows,
                         num_levels,
                         opts->getInt("gen-test-str-cols"),
                         opts->getInt("gen-test-int-cols"),
                         opts->getInt("gen-test-double-cols"));
    }

    //
    else if (opts->getBool("to-tsv5")) {
      bool do_metadata=opts->getBool("metadata");
      bool do_append=opts->getBool("append");
      // have to have an output name.
      std::string file5_out=opts->get("output");
      if (file5_out=="") {
        Err::errAbort("we need an ouput file5.");
      }
      std::string internal_name=opts->get("internal-name");
      if (internal_name!="") {
        // If a name is given, we can only import one file,
        // otherwise we would be overwriting it.
        if (opts->getArgCount()!=1) {
          Err::errAbort("can only have one input tsv when using --internal-name.");
        }
        std::string tsv_file_in=opts->getArg(0);
        TsvToTsv5_convert(tsv_file_in,
                          file5_out,
                          internal_name,
                          do_metadata,
                          do_append);
      }
      else {
        // we have a bunch-o-tsvfiles and we will
        affx::TsvFile tmp_tsv_file;
        std::string tsv_file_in;
        for (int i=0;i<opts->getArgCount();i++) {
          tsv_file_in=opts->getArg(i);
          // peek inside for the internal name.
          tmp_tsv_file.open(tsv_file_in);
          if (tmp_tsv_file.getHeader(FILE5_TSVMETA_PREFIX"-internalname",internal_name)!=affx::TSV_OK) {
            Err::errAbort("Didnt find "FILE5_TSVMETA_PREFIX"-internalname header in '"+tsv_file_in+"'");
          }
          tmp_tsv_file.close();
          //
          TsvToTsv5_convert(tsv_file_in,
                            file5_out,
                            internal_name,
                            opts->getBool("metadata"),
                            true);
        }
      }
    }
    
    //
    else if (opts->getBool("to-tsv")) {
      std::string tsv_name=opts->get("internal-name");
      if (tsv_name=="") {
        printf("Must have --internal-name=NAME with --to-tsv5. Use '--help' for an example.");
      }
      if (opts->get("output")!="") {
        if (opts->getArgCount()==0) {
          Err::errAbort("Must supply one or more input files. Use '--help' for an example.");
        }
        std::string file_in=opts->getArg(0);
        std::string file_out=opts->get("output");
        Tsv5ToTsv_convert(file_in,tsv_name,file_out,opts->getInt("precision"));
      }
      else {
        for (int i=0;i<opts->getArgCount();i++) {
          std::string file_in=opts->getArg(i);
          std::string file_out=file_in+".tsv";
          Tsv5ToTsv_convert(file_in,tsv_name,file_out,opts->getInt("precision"));
        }
      }
    }

    else if (opts->getBool("line-count")) {
      int lineCount;
      if ((opts->get("file-name")=="") || (opts->get("tsv-name")=="")) {
        Err::errAbort("Need file and tsv name.");
      }
      lineCount=affx::File5_Tsv::getFileTsvLineCount(opts->get("file-name"),
                                                     opts->get("tsv-name"));
      printf("%d\n",lineCount);
    }

    // no args, print help.
    else {
      opts->usage();
    }
    
    //
    delete opts;
    return 0;
  }
  catch (...) {
    Verbose::out(1,"Unexpected Error: uncaught exception.");
  }
  return 1;
}
