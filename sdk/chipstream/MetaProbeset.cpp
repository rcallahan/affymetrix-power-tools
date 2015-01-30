////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
#include "chipstream/MetaProbeset.h"
//
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <set>
#include <vector>
//

using namespace std;
using namespace affx;

MetaProbeset::MetaProbeset() {
  name = NULL;
  displayName = NULL;
  freeSubs = true;
}

MetaProbeset::~MetaProbeset() {
  if(freeSubs) {
    delete [] name;
    for(uint32_t i = 0; i < probesets.size(); i++) {
      delete [] probesets[i];
    }
  }
}

void MetaProbeset::readProbesetSubset(TsvFile &tsv, vector<MetaProbeset *> &metaSets) {
    string probeset;
    tsv.bind(0, "probeset_id", &probeset, TSV_BIND_REQUIRED);
    set<string> seen;
    while(tsv.nextLevel(0) == TSV_OK) {
        if(seen.find(probeset) != seen.end())
            Err::errAbort("Probeset '" + probeset + "' already seen in file: " + tsv.getFileName());
        else 
            seen.insert(probeset);
        MetaProbeset *set = new MetaProbeset();
        set->name = Util::cloneString(probeset.c_str());
        set->probesets.resize(1);
        set->probesets[0] = Util::cloneString(probeset.c_str());
        metaSets.push_back(set);
    }
}

void MetaProbeset::readMetaProbesetList(TsvFile &tsv, vector<MetaProbeset *> &metaSets) {
    string probeset, probesetList;
    vector<string> words;

    tsv.bind(0, "probeset_id", &probeset, TSV_BIND_REQUIRED);
    // list of probesets can be either called probeset_list or 
    if(tsv.cname2cidx(0, "probeset_list") != TSV_ERR_NOTFOUND)
        tsv.bind(0, "probeset_list", &probesetList, TSV_BIND_REQUIRED);
    else if(tsv.cname2cidx(0, "probeset_ids") != TSV_ERR_NOTFOUND)
        tsv.bind(0, "probeset_ids", &probesetList, TSV_BIND_REQUIRED);
    else
        Err::errAbort("Meta probeset list must contain 'probeset_list'");

    set<string> seen;
    while(tsv.nextLevel(0) == TSV_OK) {
        if(seen.find(probeset) != seen.end())
            Err::errAbort("Probeset identifier: '" + probeset + "' already seen in file: " + tsv.getFileName());
        else 
            seen.insert(probeset);
        MetaProbeset *set = new MetaProbeset();
        size_t minSize = 0;
        set->name = Util::cloneString(probeset.c_str());
        Util::chopString(probesetList, ' ', words);
        if(words.empty()) 
            Err::errAbort("Didn't get a list of sub-probesets for meta probeset '" + probeset + "'");
        minSize = words[0].size();
        /* Sanity check to make sure we got real ids and not ,,, */
        for(unsigned int i = 0; i < words.size(); i++) {
            minSize = min(words[i].length(), minSize);
        }        
        if(minSize == 0) 
            Err::errAbort("At line: " + ToStr(tsv.lineNumber()) + " got empty id in '" + probeset + "'");
        set->probesets.resize(words.size(), NULL);
        for(unsigned int i = 0; i < words.size(); i++) {
            set->probesets[i] = Util::cloneString(words[i].c_str());
        }
        metaSets.push_back(set);
    }
}


void MetaProbeset::readMetaProbesets(const string fileName, const bool metaFile, vector<MetaProbeset *> &metaSets) {
    metaSets.clear();
    TsvFile tsv;
    // open file and make sure that it has at least a probeset_id column
    if(tsv.open(fileName) != TSV_OK) 
        Err::errAbort("Couldn't open file: '" + ToStr(fileName) + "'");
    Verbose::out(2, "Opening file: " + ToStr(fileName) + " to read.");
    if(tsv.cname2cidx(0, "probeset_id") == TSV_ERR_NOTFOUND) {
        Err::errAbort("File: " + ToStr(fileName) + " must have column named 'probeset_id'. Is this the right file?");
    }

    // if file has additional probeset_list column (or deprcated
    // probeset_ids) then it is a meta probesets file, else it is a
    // subset file.
    if((tsv.cname2cidx(0, "probeset_list") != TSV_ERR_NOTFOUND) ||
        (tsv.cname2cidx(0, "probeset_ids") != TSV_ERR_NOTFOUND)) {
        if(metaFile == false) {
        Err::errAbort("File: '" + fileName + 
                    "' has both probeset_ids and probeset_list columns - should you be using meta-probesets option?");
        }
        readMetaProbesetList(tsv, metaSets);
    }
    else {
        readProbesetSubset(tsv, metaSets);
    }
    Verbose::out(2, "Read " + ToStr(metaSets.size()) + " probesets.");
    tsv.close();
}
