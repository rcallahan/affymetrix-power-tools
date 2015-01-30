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

#ifndef _CHIPSTREAM_METAPROBESET_H_
#define _CHIPSTREAM_METAPROBESET_H_

//
#include "file/TsvFile/TsvFile.h"
#include "util/FrugalVector.h"
//
#include <vector>

/** Very basic class to store the meta probesets (groups of probesets)
    that a user specifies for analysis. */
class MetaProbeset {
    public:
        const char *name;
		const char *displayName;
        FrugalVector<const char *> probesets;
        bool freeSubs;

        MetaProbeset();
        ~MetaProbeset();

        /** 
        * Read in a subset of probes to be processed. Can be probeset ids or names
        * depending on the options.
        * 
        * @param tsv - File to read from.
        * @param metaSets - Meta probeset vector to fill in.
        */
        static void readProbesetSubset(affx::TsvFile &tsv, std::vector<MetaProbeset *> &metaSets);

        /** 
        * Read in a grouping file for probesets which merges probesets into larger
        * probesets called "meta" probesets. Requires at least two columns of data,
        * a probeset_id for the newly created probeset and a probeset_list which 
        * contains a space separated list of probesets to merge.
        * 
        * @param tsv - File to read meta probesets from.
        * @param metaSets - Meta probesets to fill in from file.
        */
        static void readMetaProbesetList(affx::TsvFile &tsv, std::vector<MetaProbeset *> &metaSets);

        /** 
        * Open a file that could be either a meta probesets file or a subset list
        * of probesets and fill in the metaSets vector appropriately.
        * 
        * @param o - Options
        * @param metaSets - Vector to be filled in with resulting probests or meta probests.
        */
        static void readMetaProbesets(const std::string fileName, const bool metaFile, std::vector<MetaProbeset *> &metaSets);
};

#endif /* _CHIPSTREAM_METAPROBESET_H_ */
