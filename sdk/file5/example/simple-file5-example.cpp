////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
/*
 * Simple example of using file5 wrapper around hdf5. Note that for robust code
 * more checking of error return values and write sizes should be done, but they
 * are omitted in this example for clarity.
 */
#include <iostream>
#include "file5/File5.h"

using namespace affx;
using namespace std;

/* Write some example data to the file. */
void  writeSomeData(affx::File5_Tsv *file) {
  /* File format definition - Tell this file what columns will be
     present and what types they are. */
  file->defineColumn(0,0,"probeset", FILE5_DTYPE_STRING, 10);
  file->defineColumn(0,1,"value", FILE5_DTYPE_DOUBLE);
  
  /* Example data. */
  string probesets[] = {"ps1","ps2","ps3"};
  double values[] = {0.1,0.2,0.3};
  
  /* Write data out into file. */
  for(int i = 0; i < 3; i++) {
    file->set_string(0,0,probesets[i]);
    file->set_d(0,1,values[i]);
    file->writeLevel(0);
  }

  // reset to beginning of file.
  file->rewind();

}

/* Just print out a simple string, double file. */
void dumpSimpleFile(affx::File5_Tsv *file) {
  string s;
  float fval;
  int  ival;
  double val;
  int status = 0;
  while (file->nextLevel(0) == FILE5_OK) {
    status = file->get(0, 0, &s);
    if(status != FILE5_OK) {
      exit(1); 
    }
    cout << "As string: " << s << endl;
    status = file->get(0, 0, &fval);
    if(status != FILE5_OK) {
      exit(1); 
    }
    cout << "As float: " << fval << endl;
    status = file->get(0, 1, &val);
    if(status != FILE5_OK) {
      exit(1); 
    }
    status = file->get(0, 1, &ival);
    if(status != FILE5_OK) {
      exit(1); 
    }
    cout << "As int: " << ival << endl;
    cout << s << "\t" << val << endl;
  }
}

/* Demo reading a file */
void openAndReadOnly() {
    affx::File5_File hdf5;
    affx::File5_Group *group = NULL;

    // Open and initialize the hdf5 file 'simple-example.hd5'. 
    int status = hdf5.open("simple-example.hd5", affx::FILE5_OPEN_RO);
    if(status != FILE5_OK) {
      cout << "Bad open" << endl;
      exit(1);
    }

    // Open the /some-data group, like a directory within the file
    group = hdf5.openGroup("some-data", affx::FILE5_OPEN);

    // Open a 'tsv file' within that group. Like opening a tsv file called /some-data/some-file 
    affx::File5_Tsv *file = NULL;
    file  = group->openTsv("some-file", affx::FILE5_OPEN_RO);
  
    // Write out some example data
    cout << "Doing first file" << endl;
    dumpSimpleFile(file);

    // Open another 'tsv file' within that group. Like opening a tsv file called /some-data/some-other-file 
    affx::File5_Tsv *otherFile = NULL;
    otherFile  = group->openTsv("some-other-file", affx::FILE5_OPEN_RO);
  
    // Read some example data
    cout << "Doing second file" << endl;
    dumpSimpleFile(otherFile);

    // Cleanup
    file->close();
    delete file;
    otherFile->close();
    delete otherFile;
    group->close();
    delete group;
    hdf5.close();
}

/* Demo writing file */
void openAndWriteData() {
    affx::File5_File hdf5;
    affx::File5_Group *group = NULL;
    
    // Open and initialize the hdf5 file 'simple-example.hd5'. 
    int status = hdf5.open("simple-example.hd5", affx::FILE5_OPEN_CREATE | affx::FILE5_REPLACE);
    if(status != FILE5_OK) {
      cout << "Bad open" << endl;
      exit(1);
    }

    // Open the /some-data group, like a directory within the file
    group = hdf5.openGroup("some-data", affx::FILE5_OPEN_CREATE);

    // Open a 'tsv file' within that group. Like opening a tsv file called /some-data/some-file 
    affx::File5_Tsv *file = NULL;
    file  = group->openTsv("some-file", affx::FILE5_OPEN_CREATE);
  
    // Write out some example data
    cout << "Doing first file" << endl;
    writeSomeData(file);

    // Open another 'tsv file' within that group. Like opening a tsv file called /some-data/some-other-file 
    affx::File5_Tsv *otherFile = NULL;
    otherFile  = group->openTsv("some-other-file", affx::FILE5_OPEN_CREATE);
  
    // Read some example data
    cout << "Doing second file" << endl;
    writeSomeData(otherFile);

    // Cleanup
    file->close();
    delete file;
    otherFile->close();
    delete otherFile;
    group->close();
    delete group;
    hdf5.close();
}

/* Print usage and quit. */
void usage() {
  cout << "simple-file5-example - creates or reads an example hdf5 file\n"
  "using affymetrix file5 libraries. Try bashing the simple-example.hd5\n"
  "with a hex editor and you'll see that the correct results are produced\n"
  "or an error is reported.\n"
  "usage:\n"
  "   #create an example file\n"
  "   simple-fil5-example write\n"
  "   #read the example file\n"
    "   simple-file5-example read\n";
  exit(1);
}

/* Everybody's favorite function. Example usage*/
int main(int argc, char *argv[]) {
  string read("read");
  string write("write");
  if(argc != 2)
    usage();
  // read data
  if(argv[1] == read) {
    openAndReadOnly();
  }
  // write the data
  else if(argv[1] == write) {
    openAndWriteData();
  }
  else {
    usage();
  }
  return 0;
}
