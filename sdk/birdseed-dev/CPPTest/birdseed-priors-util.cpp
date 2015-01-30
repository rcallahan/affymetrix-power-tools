////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

#include "birdseed-dev/PriorsReader.h"
//
#include "util/AptVersionInfo.h"
#include <cstring>
#include <string.h>
#include <string>
//

//
void
usage()
{
  printf("usage\n");
}

int
main(int argc,char* argv[])
{
  if (argc<=1) {
    printf("self-test: '%s'\n",argv[0]);
  }
  //
  if (strcmp(argv[1],"-h")==0) {
    usage();
    return 0;
  }
  
  if (argc==3) {
    if (strcmp(argv[1],"-tb")==0) {
      printf("text-to-binary: '%s' -> '%s'\n",argv[2],argv[3]);
    }
    else if (strcmp(argv[1],"-bt")==0) {
      printf("binary-to-text: '%s' -> '%s'\n",argv[2],argv[3]);
    } 
    else {
      printf("unknown option: '%s'\n",argv[1]);
      return 1;
    }
  }

  //
  usage();
  return 1;
}
