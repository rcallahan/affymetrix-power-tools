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
#include "chipstream/FileInbred.h"
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

using namespace std;
// this is a monkey-see, monkey do copy of gender
// ideally, we have a general sample covariate file with multiple columns
// that allows us to read in whatever covariates and select the ones we need at any given time
// genders and inbreeding are particular special cases
// But that's probably a job for scientific programming
FileInbred::FileInbred(const string &inbredFile, const vector<string> &celFiles) {
  affx::TsvFile tsv;
  double inbred;
  string celFileName;
  map<string, double> celToInbred;
  map<string, double>::iterator iter;
  //int (*fp)(int)=std::tolower;

  tsv.open(inbredFile);
  tsv.bind(0,"inbred_het_penalty", &inbred, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"cel_files", &celFileName, affx::TSV_BIND_REQUIRED);
  while(tsv.nextLevel(0) == affx::TSV_OK) {
      celToInbred[Fs::basename(celFileName)] = inbred;
  }
  tsv.close();

  m_InbredStatus.resize(celFiles.size());
  for(int i = 0; i<celFiles.size(); i++) {
      iter = celToInbred.find(Fs::basename(celFiles[i]));
      if(iter == celToInbred.end())
      {
	      m_InbredStatus[i] = 0; // no inbreeding penalty or bonus for this file, silent failure?
      }
      else
      {
              m_InbredStatus[i] = iter->second;
      }
  }

  m_InbredName="user-supplied";
  m_InbredDescription="user supplied inbreeding status (inbred_het_penalty)";
}
