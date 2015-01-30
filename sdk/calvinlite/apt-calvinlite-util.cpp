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
// CalvinLite/apt-calvinlite-util.cpp ---
//
// $Id: apt-calvinlite-util.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_Object.h"
#include "calvinlite/CL_File.h"
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_string.h"
//
#include "util/AptVersionInfo.h"
//
#include <stdlib.h>
#include <string.h>

//////////

int print_usage(int rv)
{
  printf("apt-calvinlite-util == "
         "(version '" APT_VER_SVN_VERSION "'"
         " on '" APT_VER_COMPILE_DATE "'"
         " by '" APT_VER_COMPILE_USER "')"
         "\n"
         "--help == This help message\n"
         "\n"
         "--set-param CALVINFILE \\\n"
         "            KEY1 VAL1 [KEY_n VAL_n ...]\\\n"
         "            Replaces or adds the key value pairs.\n"
         "            This only changes the top-most (root) GDH.\n"
         " \n"
         " --print-param FILE [keys...]\n"
         " \n"
         " Print the GDH params of the file.  With out keys, prints all keys.\n"
         " With keys, just prints matching keys.\n"
         " \n"
         " The leading 'NUM:NUM' are the indexes of the nested GDH. '0' is the root GDH.\n"
         " \n"
         "--out-dir DIR == The output directory to use when dumping tsvfiles.\n"
         "                 This arg needs to be before the others.\n"
         "\n"
         "--change-param CALVINFILE STRING_FROM1 STRING_TO1 [STRING_FROM2 STRING_TO2...]\n"
         "               Changes the GDH params in a file. \n"
         "               The FROM from string is a full string match.\n"
         "               Output is to 'CALVINFILE.changed'.\n"
         "\n"
         "--change-param-substr CALVINFILE STRING_FROM1 STRING_TO1 [STRING_FROM2 STRING_TO2...]\n"
         "                      Changes the GDH params in a file. \n"
         "                      The FROM from string is a partial string match.\n"
         "                      Output is to 'CALVINFILE.changed'.\n"
         "\n"
         "--copy\n"
         "\n"
         "\n"
         "--dump\n"
         "\n"
         "--debug = Print lots of debugging output and abort as soon as an error is seen.\n"
         "\n"
         "--dump-segs\n"
         "\n"
         "\n"
         "--print-params CALVINFILE [PARAM1 PARAM2...] == print the params of the file\n"
         "\n"
#if CL_WITH_TSVFILE==1
         "--tsv-export CALVINFILE   [DATAGROUP_IDX [DATASET_IDX]]\n"
         "             Dumps data to a tsv text file.\n"
         "             The two indexes are zero-based.\n"
         "             Look at the --dump-segs output for the indexes.\n"
         "\n"
         "--tsv-export-all CALVINFILE1 CALVINFILE2...\n"
         "                 Export all the files.\n"
         "\n"
         "--tsv-import CALVINFILE TSVFILES\n"
         "             Puts the TSVFILES into the CALVINFILE.  The TSVFILES need\n"
         "             special headers like you see from the exported ones.\n"
         "             That is how the program knows where to put the data.\n"
         "             The output goes to CALVINFILE.new.\n"
         "\n"
#endif
         "--version print the version info.\n"
         "\n"
         "\n");
  return rv;
}

//////////

int do_change_file(const char* argv[],bool substr)
{
  std::string filename;
  CL_changevec_t change_vec;

  int argc=0;

  if (argv[argc]==NULL) {
    printf("Filename required.\n");
  }
  filename=argv[argc++];
  if (filename=="") {
    printf("Filename required.\n");
  }

  while (argv[argc]!=NULL) {
    std::string from=argv[argc++];
    if (argv[argc]==NULL) {
      printf("mismatched change list.  Should have matched pairs of FROM and TO.\n");
      return print_usage(1);
    }
    std::string to=argv[argc++];
    //
    change_vec.push_back(std::pair<std::string,std::string>(from,to));
  }
  return CL_File::changeFile(filename,change_vec,substr);
}

int do_setparam_file(const char* argv[])
{
  std::string filename;
  CL_changevec_t change_vec;

  int argc=0;

  if (argv[argc]==NULL) {
    printf("Filename required.\n");
  }
  filename=argv[argc++];
  if (filename=="") {
    printf("Filename required.\n");
  }

  CL_File file;
  if (file.read_File(filename)!=CL_OK) {
    return CL_ERR;
  }

  CL_Gdh* file_gdh=file.getGdh();

  // @todo should accept KEY=VAL too.
  while (argv[argc]!=NULL) {
    std::string key=argv[argc++];
    if (argv[argc]==NULL) {
      printf("mismatched change list.  Should have matched pairs of KEY VAL\n");
      return print_usage(1);
    }
    std::string val=argv[argc++];
    // get rid of the old ones.
    file_gdh->delParam(key);
    file_gdh->addParam(key,val);
  }
  
  file.write_File(filename+".changed");
  return 0;
}


int main(int argc,const char* argv[])
{
  const char** args=argv;
  std::string opt_outdir=".";

  args++;

  if ((*args==NULL)||
      (strcmp(*args,"--help")==0)||(strcmp(*args,"-h")==0)||
      (strcmp(*args,"--version")==0)||(strcmp(*args,"-v")==0)) {
    args++;
    return print_usage(0);
  }

  // perhaps this was invoked as
  // #!/usr/bin/env apt-calvinlite-util
  if (argc==2) {
    CL_Err_t rv=CL_File::tsvImportScript(argv[1]);
    return (rv==CL_OK)?0:1;
  }

  // *args != NULL after this...
  if ((*args==NULL)||
      (strcmp(*args,"--debug")==0)) {
    CL_File_debug_flags=1;
    printf("Setting 'CL_File_debug_flags' to %d\n",CL_File_debug_flags);
    args++;
  }

  //
  if ((strcmp(*args,"--out-dir")==0)||(strcmp(*args,"-o")==0)) {
    args++;
    if (*args!=NULL) {
      opt_outdir=*args;
      args++;
    }
  }

  //
  if (strcmp(*args,"--is-calvin")==0) {
    args++;
    if (*args==NULL) {
      printf("Requires Calvin filename.\n");
      return print_usage(1);
    }
    while (*args!=NULL) {
      CL_Err_t rv=CL_File::isCalvinFormat(*args);
      if (rv==CL_OK) {
        printf("'%s': ok (%d)\n",*args,rv);
      }
      else {
        printf("'%s': err (%d)\n",*args,rv);
      }
      args++;
    }
    return 0;
  }

  //
  if ((strcmp(*args,"--change-param")==0)||(strcmp(*args,"-c")==0)) {
    args++;
    return do_change_file(args,false);
  }
  if (strcmp(*args,"--change-param-substr")==0) {
    args++;
    return do_change_file(args,true);
  }

  //
  if ((strcmp(*args,"--set-param")==0)) {
    args++;
    return do_setparam_file(args);
  }

  //
  if (strcmp(*args,"--copy")==0) {
    args++;
    //
    while (*args!=NULL) {
      std::string from=*args++;
      std::string to=from+".copy";
      CL_File::copyFile(from,to);
    }
    return 0;
  }

  //
  if (strcmp(*args,"--dump")==0) {
    args++;
    while (*args!=NULL) {
      CL_File::dumpFile(*args++);
    }
    return 0;
  }
  if (strcmp(*args,"--dump-segs")==0) {
    args++;
    while (*args!=NULL) {
      CL_File::dumpFileSegs(*args++);
    }
    return 0;
  }

  //
  if (strcmp(*args,"--print-params")==0) {
    args++; // skip opt
    std::string filename=*args++;
    std::vector<std::string> param_vec;
    while (*args!=NULL) {
      param_vec.push_back(*args++);
    }
    if (param_vec.empty()) {
      CL_File::printFileParamsAll(filename);
    }
    else {
      CL_File::printFileParams(filename,param_vec);
    }
    return 0;
  }

#if CL_WITH_TSVFILE==1
  //
  if (strcmp(*args,"--tsv-export")==0) {
    args++;
    if (*args==NULL) {
      printf("Requires Calvin filename.\n");
      return print_usage(1);
    }
    std::string file=*args++;
    int dg_idx=-1;
    int ds_idx=-1;
    if (*args!=NULL) {
      dg_idx=strtol(*args++,NULL,0);
      if (*args!=NULL) {
        ds_idx=strtol(*args++,NULL,0);
      }
    }
    return CL_File::ExportToTsv(file,dg_idx,ds_idx,opt_outdir);
  }

  if (strcmp(*args,"--tsv-export-gdh")==0) {
    args++;
    if (*args==NULL) {
      printf("Requires Calvin filename.\n");
      return print_usage(1);
    }
    CL_File::ExportGdhToTsv(*args++,opt_outdir);
    return 0;
  }

  //
  if (strcmp(*args,"--tsv-export-all")==0) {
    args++;
    if (*args==NULL) {
      printf("Requires Calvin filename.\n");
      return print_usage(1);
    }
    //
    while (*args!=NULL) {
      CL_File::ExportToTsv(*args++,-1,-1,opt_outdir);
    }
    return 0;
  }

  //
  if (strcmp(*args,"--tsv-import")==0) {
    args++;
    if (*args==NULL) {
      printf("Requires Calvin filename.\n");
      return print_usage(1);
    }
    std::string file_in=*args++;
    std::string file_out=file_in+".out";
    //
    std::vector<std::string> tsv_imports;
    for (;*args!=NULL;args++) {
      std::string tmp=*args;
      tsv_imports.push_back(tmp);
    }
    //
    return CL_File::tsvImport(file_in,file_out,tsv_imports);
  }

#endif

  //
  print_usage(1);
  return 1;
}

