////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "calvin_files/fusion/src/FusionCELData.h"
//
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
//

#define NUMREADS 20
#define NUMDATA 60000

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;
using namespace std;

void readCelDataSkip(int numFiles, const char *file[]) {
  vector<int> ids;
  for(int x = 0; x < NUMDATA; x++) {
    int index = rand() % 6000000;
    ids.push_back(index);
  }
  cout << "Reading " << ids.size() << " data points." << endl;
  for(int i = 0; i < numFiles; i++) {
    FusionCELData cel;
    cel.SetFileName(file[i]);
    cel.Read();
    for(int x = 0; x < ids.size(); x++) {
      float y = cel.GetIntensity(ids[x]);
      y = x;
    }
    cel.Close();
  }
}
 
void readCelDataAll(int numFiles, const char *file[]) {
  vector<int> ids;
  for(int x = 0; x < NUMDATA; x++) {
    int index = rand() % 6000000;
    ids.push_back(index);
  }
  cout << "Reading all data points. " << endl;
  for(int i = 0; i < numFiles; i++) {
    FusionCELData cel;
    cel.SetFileName(file[i]);
    cel.Read();
    int size = cel.GetRows() * cel.GetCols();
    for(int x = 0; x < size; x++) {
      float y = cel.GetIntensity(x);
      y = x;
    }
    cel.Close();
  }
}

void readCelDataOrdered(int numFiles, const char *file[]) {
  vector<int> ids;
  for(int x = 0; x < NUMDATA; x++) {
    int index = rand() % 6000000;
    ids.push_back(index);
  }
  cout << "Reading ordered data points. " << endl;
  for(int i = 0; i < numFiles; i++) {
    ifstream in;
    in.open(file[i], ios_base::binary);
    in.seekg(1000000);
    int size = NUMDATA * 4;
    char data[size + 2];
    in.read(data, size);
    if(!in.good()) {
      cout << "Read error!" << endl; 
    }
    for(int x = 0; x < size; x++) {
      char y = data[x];
      y = 'c';
    }
  }
}
 

int main(int argc, const char *argv[]) {
  if(argc < 3) {
    cout << "file-stats - open and read a cel file different ways." << endl
         << "usage: " << endl
         << "   file-stats mode celfile.cel" << endl;
    exit(1);
  }
  srand(10);
  cout << "Doing " << argc - 2 << " cel files." << endl;
  if(string(argv[1]) == "skip") {
    cout << "Doing skips. " << endl; 
    readCelDataSkip(argc-2, argv+2);
  } 
  else if(string(argv[1]) == "all") {
    cout << "Doing read all. " << endl; 
    readCelDataAll(argc-2, argv+2);
  } 
  else if(string(argv[1]) == "ordered") {
    cout << "Doing ordered read. " << endl; 
    readCelDataOrdered(argc-2, argv+2);
  } 
  else {
    cout << "don't  recognize mode: '" << argv[1] << endl;
    exit(1);
  }
  return 0;
}


