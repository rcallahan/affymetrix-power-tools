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

#include "canary/CanaryModel.h"
//

using namespace std;

CanaryModel::CanaryModel()
  {
  std::string model_name_array[12] = {"model0","model2","model01","model12",
	"model23","model34","model012","model234","model123","model0123",
	"model1234","model01234"};

	model_names.clear();
	for (int k=0; k<12; k++) model_names.push_back(model_name_array[k]);

  std::string mstr;  // to help avoid typos

	// now the tedious work of setting up the map
	mstr = "model0";
	model_map[mstr].push_back(0);

	mstr = "model2";
	model_map[mstr].push_back(2);

	mstr = "model01";
	model_map[mstr].push_back(0);
	model_map[mstr].push_back(1);

	mstr = "model12";
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);

	mstr = "model23";
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);

	mstr = "model34";
	model_map[mstr].push_back(3);
	model_map[mstr].push_back(4);

	mstr = "model012";
	model_map[mstr].push_back(0);
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);

	mstr = "model234";
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);
	model_map[mstr].push_back(4);

	mstr = "model123";
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);

	mstr = "model0123";
	model_map[mstr].push_back(0);
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);

	mstr = "model1234";
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);
	model_map[mstr].push_back(4);

	mstr = "model01234";
	model_map[mstr].push_back(0);
	model_map[mstr].push_back(1);
	model_map[mstr].push_back(2);
	model_map[mstr].push_back(3);
	model_map[mstr].push_back(4);
	}


void CanaryModel::erase(std::string model)
	{
  for (std::vector<std::string>::iterator iter = model_names.begin();
			iter != model_names.end(); iter++)
		{
		if (*iter== model)
			{
			model_map.erase(*iter);
			model_names.erase(iter);
			return;
			}
		}

	return;
	}
