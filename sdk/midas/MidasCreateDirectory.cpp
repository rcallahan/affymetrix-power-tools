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

/// @file   MidasCreateDirectory.cpp
/// @brief  Midas directory creation utility.

#include "midas/MidasCreateDirectory.h"
//
#include "util/Convert.h"
//
#include <sys/stat.h>
#include <sys/types.h>
//

#ifdef WIN32
#include <direct.h>
#endif

using namespace std;

// create output directory if not already present
std::string midasCreateDirectory (const std::string& outDirectory)
{
#ifdef WIN32
  // windows: use ISO C++ conformants; POSIX functions are deprecated
  struct _stat outDirStat;
  int statError = _stat (outDirectory.c_str(), &outDirStat);
#else
  struct stat outDirStat;
  int statError = stat (outDirectory.c_str(), &outDirStat);
#endif

  // error return from fstat - requested output directory
  // doesn't exist, so we create one
  std::string msg = "";
  if (statError != 0)
  {
  #ifdef WIN32
    // windows: use ISO C++ _mkdir; umask is deprecated
    int mkdirError = _mkdir (outDirectory.c_str());
  #else
    mode_t Umask = umask (0);
    // set directory execute permissions based on current umask
    mode_t Mode = S_IRWXU;	// user (this code) can read, write, execute
    Mode |= (~Umask & (S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
    if (Mode & (S_IRGRP | S_IWGRP)) // group, other access is optional
      Mode |= S_IXGRP;
    if (Mode & (S_IROTH | S_IWOTH))
      Mode |= S_IXOTH;
    int mkdirError = mkdir (outDirectory.c_str(), Mode);
    // restore caller's umask
    umask (Umask);
  #endif
    if (mkdirError != 0)
    {
      // report error
      msg = "Unable to create directory " + outDirectory;
      return msg;
    }
  }
#ifndef WIN32
  // requested file exists - just check if it's a directory
  else if (! S_ISDIR (outDirStat.st_mode) )

  {
    // report error
    msg = "Requested output directory " + outDirectory + " is not a directory";
    return msg;
  }
#endif
  // no error found
  return msg;
}
