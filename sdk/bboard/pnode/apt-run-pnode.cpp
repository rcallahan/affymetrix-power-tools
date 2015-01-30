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
//
// ~/Affy-projects/apt2/trunk/apt-run-pile.cpp ---
//
// $Id: apt-run-pnode.cpp,v 1.4 2009-11-05 20:41:39 harley Exp $
//

// These should be collected in '#include "apt-main-includes.h"'
#include "bboard/Bboard.h"
#include "bboard/pnode/Pnode.h"
#include "bboard/pnode/PnodeFactory.h"

//
#include <string>

//////////

// This program is an example of running one Pnode.
// * reads in a BB
// * runs one Pnode
// * Writes out the BB.

// It could be used to chain together a series of Pnodes
// in a shell script.
//
//               (in)       (Pnode)        (out)
// apt-run-pnode data-1.bb  PN_Normalize   data-2.bb
// apt-run-pnode data-2.bb  PN_GenderCall  data-3.bb
// apt-run-pnode data-3.bb  PN_Estimate    data-4.bb
// apt-run-pnode data-4.bb  PN_WriteCalvin data-6.bb

int main(int argc,const char** argv)
{
  if (argc!=4) {
    return APT_ERR;
  }

  // The three args.
  std::string bb_i_filename=argv[1];
  std::string pnode_name=argv[2];
  std::string bb_o_filename=argv[3];
  
  // our BB
  Bboard* bb=new Bboard("main");
  // our Pnode
  Pnode* pnode=PnodeFactory::fromString(pnode_name);

  // now do the three steps: in, process, out.
  printf("== reading BB from '%s'...\n",bb_i_filename.c_str());
  DaoUtil::readFromFile(bb_i_filename,bb);
  printf("== read:\n");
  bb->dump();

  printf("== running node '%s'...\n",pnode_name.c_str());
  pnode->doRun(bb);
  printf("== output:\n");
  bb->dump();

  printf("== writing BB to '%s'...\n",bb_o_filename.c_str());
  DaoUtil::writeToFile(bb_o_filename,bb);

  // and clean up.
  delete bb;
  delete pnode;
  //
  return APT_OK;
}
