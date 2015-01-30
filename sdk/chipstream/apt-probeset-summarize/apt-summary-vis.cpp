////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/**
/// @file   apt-summary-vis.cpp
/// @brief  Main for creating IGB egr files from summary and annotation files.
*/

#include "chipstream/apt-probeset-summarize/SummaryVis.h"
#include "util/AptVersionInfo.h"
//

using namespace std;


int main (int argc,const char* argv[])
{
  try {
    const string version = AptVersionInfo::versionToReport();

    try
    {
        summaryVis summary (argc, argv, version);
        summary.run();
    }
    
    catch (exception& e)
    {
        // Verbose::out (1, e.what()); // errAbort() writes the error message.
        return 1;
    }
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
