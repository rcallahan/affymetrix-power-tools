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
// affy/sdk/bboard/pnode/Apt2Main.h ---
//
// $Id$
//

//
#include <stdio.h>
#include <string>
#include <vector>

/// A standard main which is used to invoke our "runMain".
#define CALL_APT2MAIN(_Class)                   \
  int main(const int argc,const char** argv) {  \
    _Class main_obj;                            \
    return main_obj.internalRunMainArgv(argc,argv);     \
  }

class Apt2Main {
public:
  // simple no arg constructor. (use in CALL_APT2MAIN)
  Apt2Main();

  // this method is used by CALL_APT2MAIN macro to do the invoking of our main.
  // does the per-linux-process setup.
  int internalRunMainArgv(const int argc,const char** argv);
  // this is what people should overload
  // should be pure virtual but isnt for testing.
  virtual int runMain(const std::vector<std::string>& argv);
  // this method wouldnt do any setup, but just run so one "Apt2Main" 
  // could call another without repeating the per-linux-process setup.
  // virtual int run(const std::vector<std::string>& argv);
  // here we can add a bunch of methods for common setup activites.
  //int defineOptions();
  // 
  // releaseAllResources()
  //
  // abort();
  //
  // exit();
};
