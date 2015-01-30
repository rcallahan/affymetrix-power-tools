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
// ~/Affy/apt2/trunk/Pnode.cpp ---
// 
// $Id: Pnode.cpp,v 1.3 2009-11-05 20:41:39 harley Exp $
// 

//
#include "Pnode.h"

//////////

Pnode::Pnode()
{
  init();
}
Pnode::Pnode(const std::string& name)
{
  m_name=name;
}
Pnode::~Pnode()
{
  doDeleteChildren();
}

void Pnode::init() {
  // this should be the setup action.
}

//////////

//
Pnode* Pnode::addPnode(Pnode* pnode)
{
  m_pnode_vec.push_back(pnode);
  return pnode;
}

void Pnode::doDeleteChildren() {
  for (int i=0;i<m_pnode_vec.size();i++) {
    delete m_pnode_vec[i];
  }
  m_pnode_vec.clear();
}

//////////

//
AptErr_t Pnode::doInitNode(Bboard* bb) {
  // no-op
  return APT_OK;
}
AptErr_t Pnode::doClearNode(Bboard* bb) {
  // no-op
  return APT_OK;
}
AptErr_t Pnode::doDumpNode(Bboard* bb,const std::string& prefix) {
  // we should dump more.
  printf("dump: %s %s\n",prefix.c_str(),m_name.c_str());
  return APT_OK;
}

//////////

AptErr_t Pnode::doRunNodePre(Bboard* bb) {
  // no-op
  return APT_OK;
}
AptErr_t Pnode::doRunNode(Bboard* bb) {
  // no-op
  return APT_OK;
}
AptErr_t Pnode::doRunNodePost(Bboard* bb) {
  // no-op
  return APT_OK;
}

//////////

//
AptErr_t Pnode::doInitChildren(Bboard* bb) {
  for (int i=0;i<m_pnode_vec.size();i++) {
    m_pnode_vec[i]->doInit(bb);
  }
  //
  return APT_OK;
}
AptErr_t Pnode::doClearChildren(Bboard* bb) {
  for (int i=0;i<m_pnode_vec.size();i++) {
    m_pnode_vec[i]->doClear(bb);
  }
  //
  return APT_OK;
}
AptErr_t Pnode::doRunChildren(Bboard* bb) {
  for (int i=0;i<m_pnode_vec.size();i++) {
    m_pnode_vec[i]->doRun(bb);
  }
  //
  return APT_OK;
}
AptErr_t Pnode::doDumpChildren(Bboard* bb,const std::string& prefix) {
  for (int i=0;i<m_pnode_vec.size();i++) {
    m_pnode_vec[i]->doDump(bb,prefix);
  }
  //
  return APT_OK;
}

//
AptErr_t Pnode::doInit(Bboard* bb) {
  doInitNode(bb);
  doInitChildren(bb);
  //
  return APT_OK;
}
AptErr_t Pnode::doClear(Bboard* bb) {
  doClearNode(bb);
  doClearChildren(bb);
  //
  return APT_OK;
}

// This is the default series of actions in the correct order.
// Most "PN_" nodes wont need to change this, they will just want to
// supply their own: doRunNodePre(), doRunNode() or doRunNodePost()
// methods.
AptErr_t Pnode::doRun(Bboard* bb) {
  //
  doRunNodePre(bb);
  doRunNode(bb);
  doRunChildren(bb);
  doRunNodePost(bb);

  // we might want to have the BB dumped after the node has run.
  if (m_dump_post_filename!="") {
    DaoUtil::writeToFile(m_dump_post_filename,bb);
  }
  //
  return APT_OK;
}

AptErr_t Pnode::doDump(Bboard* bb,const std::string& prefix) {
  doDumpNode(bb,prefix);
  doDumpChildren(bb,prefix+" "+m_name);
  //
  return APT_OK;
}
