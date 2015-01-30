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

#include "chipstream/DataStore.h"

using namespace std;

void usage() {
  cout << "cel-to-data-store - Program to put cel files into an hdf5 DataStore representation" << endl
       << "  usage:" << endl
       << "  cel-to-data-store cels.f5 cel1.cel cel2.cel ... celN.cel" << endl;
  exit(1);
}

int main(int argc, char* argv[]) {
  Verbose::out(1, "Making map");
  int numProbes = DiskIntensityMart::getProbeCountFromCel(argv[2]);
  vector<int> map(numProbes);
  for(int i = 0; i < numProbes; i++) {
    map[i] = i;
  }
  Verbose::out(1, "Done.");
  Verbose::out(1, "Making datastore");
  DataStore ds(argv[1]);
  Verbose::out(1, "Done.");
  ds.setProbeOrder(map);
  vector<string> files;
  for(int i = 2; i < argc; i++) {
    files.push_back(argv[i]);
  }
  //  files[0] = "../../regression-data/data/cel/HuEx-1_0-st-v2/huex_wta_cb_A.CEL";
  //  files[1] = "../../regression-data/data/cel/HuEx-1_0-st-v2/huex_wta_cb_C.CEL";
  Verbose::out(1, "Reading cel files.");
  ds.slurpDataFromCelFiles(files);
  Verbose::out(1, "Done.");
  return 0;
}
