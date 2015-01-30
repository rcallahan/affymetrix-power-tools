////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#ifndef _CANARYMODEL_H
#define _CANARYMODEL_H

//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

class CanaryModel {

public:
  CanaryModel();

  std::string operator()(int k) { return model_names[k]; }

  std::vector<int> operator[](std::string mstr) { return model_map[mstr]; }
  
  unsigned long size() { return model_names.size(); }

  void erase(std::string model);

private:
  std::vector<std::string> model_names;
  std::map<std::string,std::vector<int> > model_map;
};

#endif
