////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
#include "chipstream/FileGenders.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

using namespace std;

FileGenders::FileGenders(const string &genderFile, const vector<string> &celFiles) {
  affx::TsvFile tsv;
  string gender;
  string celFileName;
  map<string, affx::Gender> celToGender;
  map<string, affx::Gender>::iterator iter;
  int (*fp)(int)=std::tolower;

  tsv.open(genderFile);
  tsv.bind(0,"gender", &gender, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"cel_files", &celFileName, affx::TSV_BIND_REQUIRED);
  while(tsv.nextLevel(0) == affx::TSV_OK) {
      affx::Gender typedGender = affx::UnknownGender;
      transform(gender.begin(), gender.end(), gender.begin(), fp);
      if(gender == "male")
        typedGender = affx::Male;
      else if(gender == "1")
        typedGender = affx::Male;
      else if(gender == "female")
        typedGender = affx::Female;
      else if(gender == "2")
        typedGender = affx::Female;
      celToGender[Fs::basename(celFileName)] = typedGender;
  }
  tsv.close();

  m_Genders.resize(celFiles.size());
  for(int i = 0; i<celFiles.size(); i++) {
      iter = celToGender.find(Fs::basename(celFiles[i]));
      if(iter == celToGender.end())
          Err::errAbort("Unable to find cel file " + Fs::basename(celFiles[i]) + " in gender file " + ToStr(genderFile));
      m_Genders[i] = iter->second;
  }

  m_GenderName="user-supplied";
  m_GenderDescription="user supplied genders";
}
