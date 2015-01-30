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
// affy/sdk/bboard/pnode/apt2-main-test1.cpp ---
//
// $Id$
//

// g++ -o apt2-main-test1 apt2-main-test1.cpp Apt2Main.cpp

//
#include "Apt2Main.h"
#include <stdexcept>

// an example of a main
class Apt2Main_Test1 : public Apt2Main {
  int runMain(const std::vector<std::string>& argv) {
    // test the exception stuff
    if (argv.size()==0) {
      throw std::runtime_error("zero args");
    }
    // do something to show we ran.
    for (int i=0;i<argv.size();i++) {
      printf("argv[%3d] : '%s'\n",i,argv[i].c_str());
    }
  }
};

//
CALL_APT2MAIN(Apt2Main_Test1);

