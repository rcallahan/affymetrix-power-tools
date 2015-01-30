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

/// @file   apt-cdf-to-spf.cpp
/// @brief  converts a CDF file to an SPF file.
/// @todo   rename this to "apt-cdf-to-spf.cpp"

//
#include "chipstream/ChipLayout.h"
#include "chipstream/EngineUtil.h"
//
#include "util/Err.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//

using namespace std;

void define_cdftospf_options(PgOptions* opts)
{
  opts->setUsage("apt-cdf-to-spf - Program for converting cdf files to spf (simple probe format).\n"
                 "Will also convert pgf/clf files into spf format as well.\n"
                 "\n"
                 "Usage:\n"
                 "   apt-cdf-to-spf --cdf-file file.cdf --spf-file outfile.spf\n"
                 "   apt-cdf-to-spf --pgf-file file.pgf --clf-file file.clf --spf-file outfile.spf");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Display program options short blurb on usage.",
                     "false");
  opts->defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets, user either --cdf-file or --pgf-file and --clf-file",
                     "");
  opts->defineOption("p", "pgf-file", PgOpt::STRING_OPT,
                     "File defining probe sets.",
                     "");
  opts->defineOption("l", "clf-file", PgOpt::STRING_OPT,
                     "File defining x,y <-> probe id conversion.",
                     "");
  opts->defineOption("s", "spf-file", PgOpt::STRING_OPT,
                     "File to be written with probe sets in spf (simple probe format).",
                     "");
  opts->defineOption("","spf-format", PgOpt::INT_OPT,
                     "Format of spf to write. "
                     "Currently only v2 is allowed. (will support v3 & v4 later.)",
                     "2");
}

int main(int argc,const char* argv[]) {
  try {
    PgOptions *opts = NULL;
    opts = new PgOptions();
    define_cdftospf_options(opts);
    opts->parseArgv(argv);
    
    if(opts->getBool("help") || argc == 1) {
        opts->usage(true);
        exit(0);
    }

    std::string spfOut = opts->get("spf-file");
    int spfFormat = opts->getInt("spf-format");

    ChipLayout layout;

    if(spfOut=="")
        Err::errAbort("Must supply --spf-file for output");
    if (opts->get("cdf-file")!="") {
        Verbose::out(1, "Reading cdf file.");
        if(!layout.openCdfAll(opts->get("cdf-file"))) {
            Err::errAbort("Couldn't open layout file: '" + ToStr(opts->get("cdf-file")) + "'");
        }
    }
    else if ((opts->get("pgf-file")!="") && 
            (opts->get("clf-file")!="")) {
        Verbose::out(1, "Reading pgf file.");
        const std::set<const char *, Util::ltstr> probeSetsToLoad;
        std::vector<bool> probeSubset;
        colrow_t  rows = 0, cols = 0;
        int probeCount = 0;
        std::vector<std::string> chipTypes;
        EngineUtil::getPgfChipType(chipTypes, rows, cols, probeCount, 
                                opts->get("pgf-file"), opts->get("clf-file"));
        const std::string chipType;
        if(!layout.openPgf(opts->get("pgf-file"), rows, cols, probeSetsToLoad, probeSubset, chipTypes[0])) {
        Err::errAbort("Couldn't open layout file: '" + ToStr(opts->get("pgf-file")) + "'");
        }
    }
    else {
        Err::errAbort("Must specify either a cdf file or a pgf/clf file pair.");
    }
    Verbose::out(1, "Writing spf file.");
    layout.writeSpfProbeList(spfOut,spfFormat);
    Verbose::out(1, "Done.");
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
