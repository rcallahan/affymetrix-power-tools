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

/// @file   apt-midas.cpp
/// @brief  main for Midas C++ command line program

//
#include "midas/MidasConfigurePgOptions.h"
#include "midas/MidasConfigureRun.h"
#include "midas/MidasCreateDirectory.h"
#include "midas/MidasEngine.h"
//
#include "util/AptVersionInfo.h"
#include "util/Util.h"
//
#include <iostream>
//

using namespace std;


int main (int argc, char* argv[])
{
  try {
    const string version = AptVersionInfo::versionToReport();

    // PgOptions are globals - this constructor alters them, so it can only be invoked once
    PgOptions opts;
    define_midas_opts(&opts);
    opts.parseArgv(argv);
    
    // optionally display usage message
    if (opts.getBool("help") || argc == 1)
    {
        opts.usage();
        cout << "version: " << version << endl;
        return 0;
    }
    // optionally display version
    if (opts.getBool("version"))
    {
        cout << "version: " << version << endl;
        return 0;
    }
    
    // check output options
    bool wantPvalues = opts.getBool("pvalues");
    const bool wantFstats = opts.getBool("fstats");
    const bool wantNormalized = opts.getBool("normalized");
    
    // require at least one output
    if (! wantPvalues && ! wantFstats && ! wantNormalized)
    {
        Err::errAbort("No outputs selected - please choose at least one");
    }
    // check requested output directory, if any - create
    // one if it's not already present
    std::string outDir = opts.get("out-dir");
    if (outDir != "")
    {
        std::string msg = midasCreateDirectory (outDir);
        // report fatal error if any
        if (msg!="")
            Err::errAbort(ToStr(msg));
    }
    try
    {
        const string execVersion = version;
        // create object to configure, run apt-midas
        midasConfigureRun configureRun ( opts.get("cel-files"), opts.get("genedata"),
                                        opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                        wantPvalues, wantFstats, wantNormalized, opts.getDouble("stabilize"),
                                        opts.commandLine(), execVersion, opts.getBool("no-logtrans"), opts.getBool("keep-path") );
        // configure step may return a non-fatal warning message
        std::string* msg = configureRun.configure();
        if (msg)
            Err::errAbort(ToStr(*msg));
        configureRun.run();
    }
    // report errors
    catch (exception& e)
    {
        Err::errAbort(e.what());
    }
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
