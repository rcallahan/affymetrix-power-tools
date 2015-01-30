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
//  ~/Affy/apt2/trunk/Pile.h ---
// 
//  $Id: Pnode.h,v 1.1 2009-10-28 18:25:58 harley Exp $
// 

#ifndef _PNODE_H_
#define _PNODE_H_

//
#include "bboard/Apt2Types.h"
#include "bboard/Bboard.h"
//
#include <string>
#include <vector>

class Pnode {
public:
  // name for debugging uses. (not used by anyone else.)
  std::string m_name;

  // child nodes.
  std::vector<Pnode*> m_pnode_vec;

  // if this is set, the default "doRun" will dump the BB
  // to this filename at the end of "doRun"...
  std::string m_dump_pre_filename;
  std::string m_dump_post_filename;

  Pnode();
  Pnode(const std::string& name);

  virtual ~Pnode();
  virtual void init();
  virtual void doDeleteChildren();

  // This adds a Pnode as a child node.
  // It also takes ownership of the node and will delete them
  // when this Pnode is deleted.
  Pnode* addPnode(Pnode* node);

  // These functions do the step for the node and their children.
  virtual AptErr_t doInit(Bboard* bb);
  virtual AptErr_t doClear(Bboard* bb);
  // If you do something funky, you may want to override this one.
  virtual AptErr_t doRun(Bboard* bb);
  virtual AptErr_t doDump(Bboard* bb,const std::string& prefix);

  // These functions do this step for this node only.
  // These are the ones you will normally override in normal ("PN_") use.

  // setup.
  virtual AptErr_t doInitNode(Bboard* bb);
  // teardown.
  virtual AptErr_t doClearNode(Bboard* bb);
  /// before running children
  virtual AptErr_t doRunNodePre(Bboard* bb);
  /// also before running childen, intended for leaf nodes. (It has a better name.)
  virtual AptErr_t doRunNode(Bboard* bb);
  /// after the children have been run.
  virtual AptErr_t doRunNodePost(Bboard* bb);
  /// This should dump the state of the node to stdout.
  virtual AptErr_t doDumpNode(Bboard* bb,const std::string& prefix);

  // These are stock functions to invoke a Step in the child nodes.
  // (So you dont have to write a loop.)
  virtual AptErr_t doInitChildren(Bboard* bb);
  virtual AptErr_t doClearChildren(Bboard* bb);
  virtual AptErr_t doRunChildren(Bboard* bb);
  virtual AptErr_t doDumpChildren(Bboard* bb,const std::string& prefix);
};

#endif // _PNODE_H_
