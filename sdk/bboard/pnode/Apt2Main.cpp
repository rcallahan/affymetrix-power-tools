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
// affy/sdk/bboard/pnode/Apt2Main.cpp ---
//
// $Id$
//

//
#include "Apt2Main.h"
//
#include <exception>
#include <stdexcept>
#include <stdio.h>

Apt2Main::Apt2Main() {
  // print a message for testing.
  printf("Apt2Main\n");
}

int Apt2Main::internalRunMainArgv(const int argc,const char** argv) {
  // change the args into a something more STL-like.
  std::vector<std::string> argv_vec;
  for (int i=1;i<argc;i++) {
    argv_vec.push_back(argv[i]);
  }
  // do some setup?
  // ... more setup here.
  int rv=0;
  try {
    // run our users subclassed main.
    rv=runMain(argv_vec);
  }
  catch (std::exception& e) {
    // do something with the exception.
    std::string what=e.what();
    printf("An exception occured!!!\n");
    printf("message: '%s'\n",what.c_str());
    // 
    // if Apt2Main had a root blackboard slot, then
    // we would make sure the resources of the BB were freed.
    rv=-1;
  }

  return rv;
}

// this should be pure virtual, but isnt just for testing.
int Apt2Main::runMain(const std::vector<std::string>& argv) {
  printf("Apt2Main::runMain(): hello world\n"
         "this should be over-ridden.\n");
}
