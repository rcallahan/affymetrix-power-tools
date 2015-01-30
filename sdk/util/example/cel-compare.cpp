////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "util/CelCheck.h"
//
#include <cstdlib>
#include <iostream>
//

using namespace std;

int main(int argc, char *argv[]) {
  if(argc != 3) {
    cout << "usage: " << endl;
    cout << "   cc-chp-compare test.chp gold.chp" << endl;
    exit(0);
  }
  vector<string> test, gold;
  test.push_back(argv[1]);
  gold.push_back(argv[2]);
  cout << "Getting ready to check. " << endl;
  CelCheck checker(test, gold, 0.1f);
  string msg;
  if(!checker.check(msg))
    cout << "File not equivalent!!! " << msg << endl;
  else
    cout << "Files equivalent." << endl;
  return 0;
}
