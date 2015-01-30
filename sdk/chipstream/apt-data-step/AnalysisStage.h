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
#ifndef _ANALYSISSTAGE_H_
#define _ANALYSISSTAGE_H_

#include <map>
#include <string>
#include <vector>

class AnalysisStage {

public: 

  static std::string nameFromSpec(const std::string &spec) {
    std::string name;
    if(spec.find(".") != string::npos) {
      name = spec.substr(0, spec.find("."));
    }
    else {
      name = spec;
    }
    return name;
  }

  static std::vector<AnalysisStage> makeStages(const string &spec, 
                                               std::map<std::string,std::string> &stdMethods) {
    assert(!spec.empty());
    std::string description;
    bool doingStdMethod = false;
    std::vector<std::string> words;

    /* Check for one of our std methods. */
    if(stdMethods.find(spec) != stdMethods.end()) {
      doingStdMethod = true;
      description = stdMethods[spec];
    }
    else {
      description = spec;
    }

    Util::chopString(description, ',', words);
    if(words.size() < 2) {
      Err::errAbort("Must specify at least a pm adjustment and summary type.");
    };

    std::vector<AnalysisStage> stages;
    // The first n-2 steps are chipstreams
    for(size_t i = 0; i < words.size(); i++) {
      AnalysisStage cStage;
      cStage.spec = words[i];
      cStage.name = nameFromSpec(cStage.spec);
      stages.push_back(cStage);
    }

    return stages;
  }

  std::string spec;
  std::string name;
  
};

#endif // _ANALYSISSTAGE_H_
