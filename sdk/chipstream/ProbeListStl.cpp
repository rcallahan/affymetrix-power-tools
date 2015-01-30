////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "chipstream/ProbeListStl.h"
//
#include "util/Err.h"

////////// PROBE

ProbeListStlProbe::ProbeListStlProbe()
{
  clear();
}

ProbeListStlProbe::ProbeListStlProbe(const ProbeListStlProbe& p)
{
  // no need to clear, as all are set.
  m_probeId   = p.m_probeId;
  m_probeGc   = p.m_probeGc;
  m_probeApid = p.m_probeApid;
}

ProbeListStlProbe::ProbeListStlProbe(int id)
{
  clear();
  m_probeId=id;
}

ProbeListStlProbe::ProbeListStlProbe(int id,int gc,int apid)
{
  m_probeId   = id;
  m_probeGc   = gc;
  m_probeApid = apid;
}

void
ProbeListStlProbe::clear()
{
  // "-1" is unset.
  m_probeId   = ProbeList::UNSET;
  m_probeGc   = ProbeList::UNSET;
  m_probeApid = ProbeList::UNSET;
}

int
ProbeListStlProbe::id()
{
  return m_probeId;
}

int
ProbeListStlProbe::gc()
{
  return m_probeGc;
}

int
ProbeListStlProbe::apid()
{
  return m_probeApid;
}

bool
ProbeListStlProbe::isNull()
{
  return m_probeId==ProbeList::UNSET;
}

////////// BLOCK

ProbeListStlBlock::ProbeListStlBlock() {
  clear();
}

ProbeListStlBlock::ProbeListStlBlock(int bAnn,int bAllele, int bContext, int bChannel, int bRepType)
{
  m_blockSize=0;
  m_blockAnnotation=bAnn;
  m_blockAllele=bAllele;
  m_blockContext=bContext;
  m_blockChannel=bChannel;
  m_blockRepType=bRepType;
}

ProbeListStlBlock::ProbeListStlBlock(const ProbeListStlBlock& block) {
  m_blockSize=block.m_blockSize;
  m_blockAnnotation=block.m_blockAnnotation;
  m_blockAllele=block.m_blockAllele;
  m_blockContext=block.m_blockContext;
  m_blockChannel=block.m_blockChannel;
  m_blockRepType=block.m_blockRepType;
}

void
ProbeListStlBlock::clear()
{
  m_blockSize=0;
  //
  m_blockAnnotation=ProbeList::UNSET;
  m_blockAllele=ProbeList::UNSET;
  m_blockContext=ProbeList::UNSET;
  m_blockChannel=ProbeList::UNSET;
  m_blockRepType=ProbeList::UNSET;
}

bool ProbeListStlBlock::matchACC(int qAllele,int qContext,int qChannel) const
{
  if ((qAllele==m_blockAllele) &&
      (qContext==m_blockContext) &&
      (qChannel==m_blockChannel)) {
    return true;
  }
  return false;
}

////////// LIST

ProbeListStl::ProbeListStl() {
  clear();
}

ProbeListStl::ProbeListStl(const ProbeListPacked& plp)
{
  // use operator= to make a copy.
  *this=plp;
}

ProbeListStl::ProbeListStl(const ProbeListStl& pls)
{
  m_name=pls.m_name;
  m_max_name_len=pls.m_max_name_len;
  m_type=pls.m_type;
  m_num_match=pls.m_num_match;
  m_block_vec=pls.m_block_vec;
  m_probe_vec=pls.m_probe_vec;
}

ProbeListStl::ProbeListStl(int block_cnt,int probe_cnt)
{
  clear();
  resize(block_cnt,probe_cnt);
}
ProbeListStl::ProbeListStl(int block_cnt,int probe_cnt,const std::string& name)
{
  clear();
  resize(block_cnt,probe_cnt);
  m_name=name;
}

ProbeListStl::ProbeListStl(const std::string& name)
{
  clear();
  m_name=name;
}

ProbeListStl& ProbeListStl::operator=(const ProbeListPacked& plp)
{
  int b_cnt=plp.block_cnt();
  int p_cnt=plp.probe_cnt();

  //printf("operator= ----------\n");

  clear();
  reserve(b_cnt,p_cnt);

  //
  set_name(plp.get_name_string());
  set_type(plp.get_type());
  set_numMatch(plp.get_numMatch());
  //
  for (int b=0;b<b_cnt;b++) {
    push_block(plp.get_blockAnn(b),
               plp.get_blockAllele(b),
               plp.get_blockContext(b),
               plp.get_blockChannel(b),
               plp.get_blockRepType(b));
    set_blockSize(b,plp.get_blockSize(b));
  }
  // set the Apid of each probe as we copy it.
  // we know they are sequental from the block start.
  // this will save looking it up later.
  int apid=plp.getApidStart();
  APT_ERR_ASSERT(apid>=0,"");

  for (int p=0;p<p_cnt;p++) {
    push_probe(plp.get_probeId(p),plp.get_probeGc(p),apid++);
  }

  //printf("operator=\n");
  //plp.dump();
  //dump();

  return *this;
}

void ProbeListStl::resize(int block_cnt,int probe_cnt) {
  m_block_vec.resize(block_cnt);
  m_probe_vec.resize(probe_cnt);
}
void ProbeListStl::reserve(int block_cnt,int probe_cnt) {
  m_block_vec.reserve(block_cnt);
  m_probe_vec.reserve(probe_cnt);
}

void ProbeListStl::clear() {
  m_name="";
  m_max_name_len=1024; // the most we want to store.
  m_type=ProbeList::UNSET;
  m_num_match=0;
  //
  m_block_vec.resize(0);
  m_probe_vec.resize(0);
}

//
int ProbeListStl::get_type() const {
  return m_type;
}
void ProbeListStl::set_type(int type) {
  m_type=type;
}

// accessors
int ProbeListStl::block_cnt() const {
  return m_block_vec.size();
}
int ProbeListStl::probe_cnt() const {
  return m_probe_vec.size();
}
//
int ProbeListStl::get_numMatch() const {
  APT_ERR_ASSERT(m_num_match>=0,"internal error.");
  return m_num_match;
}
void ProbeListStl::set_numMatch(int numMatch) {
  /// @todo this is for debugging.  Maybe we could have 0 matches.
  APT_ERR_ASSERT(numMatch>0,"internal error.");
  m_num_match=numMatch;
}
//
int ProbeListStl::get_blockSize(int i) const {
  return m_block_vec[i].m_blockSize;
}
void ProbeListStl::set_blockSize(int i,int val) {
  m_block_vec[i].m_blockSize=val;
}
//
short ProbeListStl::get_blockAnn(int i) const {
  return m_block_vec[i].m_blockAnnotation;
}
void ProbeListStl::set_blockAnn(int i,short val) {
  m_block_vec[i].m_blockAnnotation=val;
}
//
int ProbeListStl::get_blockAllele(int i) const {
  return m_block_vec[i].m_blockAllele;
}
void ProbeListStl::set_blockAllele(int i,int val) {
  m_block_vec[i].m_blockAllele=val;
}
//
int ProbeListStl::get_blockContext(int i) const {
  return m_block_vec[i].m_blockContext;
}
void ProbeListStl::set_blockContext(int i,int val) {
  m_block_vec[i].m_blockContext=val;
}
//
int ProbeListStl::get_blockChannel(int i) const {
  return m_block_vec[i].m_blockChannel;
}
void ProbeListStl::set_blockChannel(int i,char val) {
  m_block_vec[i].m_blockChannel=val;
}
//
int ProbeListStl::get_blockRepType(int i) const {
  return m_block_vec[i].m_blockRepType;
}
void ProbeListStl::set_blockRepType(int i,char val) {
  m_block_vec[i].m_blockRepType=val;
}
//
int ProbeListStl::get_probeId(int i) const {
  return m_probe_vec[i].m_probeId;
}
void ProbeListStl::set_probeId(int i,int val) {
  m_probe_vec[i].m_probeId=val;
}
//
int8_t ProbeListStl::get_probeGc(int i) const {
  return m_probe_vec[i].m_probeGc;
}
void  ProbeListStl::set_probeGc(int i,int8_t val) {
  m_probe_vec[i].m_probeGc=val;
}

//
int ProbeListStl::get_probeApid(int i) const {
  return m_probe_vec[i].m_probeApid;
}
void  ProbeListStl::set_probeApid(int i,int val) {
  m_probe_vec[i].m_probeApid=val;
}

//
void ProbeListStl::set_name(const std::string& name) {
  m_name=name;
}
// this shouldnt be used.
const char* ProbeListStl::get_name() const {
  return m_name.c_str();
};
// use this instead.
const char* ProbeListStl::get_name_cstr() const {
  return m_name.c_str();
}
//
std::string ProbeListStl::get_name_string() const {
  return m_name;
}
//
int ProbeListStl::max_name_len() const {
  return m_max_name_len;
}

//
void ProbeListStl::push_block(const ProbeListStlBlock& b)
{
  m_block_vec.push_back(b);
}
void ProbeListStl::push_block(const ProbeListStlBlock* b)
{
  APT_ERR_ASSERT(b!=NULL,"Cant push back a null block.");
  m_block_vec.push_back(*b);
}
void ProbeListStl::push_block(int bAnn,int bAllele, int bContext, int bChannel, int bRepType)
{
  m_block_vec.push_back(ProbeListStlBlock(bAnn,bAllele,bContext,bChannel,bRepType));
}

/// @todo blocks should contain probes.
/// we need the pls ref as we need it to get the probes.
/// 
/*
void ProbeListStl::push_blockandprobesfrom(const ProbeListStl& pls_from,int bidx)
{
  m_block_vec.push_back(pls_from.getBlockPtr(bidx));
  int pbase=bidx*m_block_vec[bidx].m_blockSize();
  for (int pidx=0;pls_from.m_block_vec[bidx].m_blockSize();pidx++) {
    push_probe(pls_from.get_probeId(pbase+pidx));
  }
}
*/
//
void ProbeListStl::push_probe(const ProbeListStlProbe& p)
{
  m_probe_vec.push_back(p);
}
void ProbeListStl::push_probe(int id)
{
  m_probe_vec.push_back(ProbeListStlProbe(id));
}
void ProbeListStl::push_probe(int id,int gc,int apid)
{
  m_probe_vec.push_back(ProbeListStlProbe(id,gc,apid));
}

//
void ProbeListStl::push_probelist(const ProbeListStl& pl)
{
  for (int b=0;b<pl.block_cnt();b++) {
    m_block_vec.push_back(pl.m_block_vec[b]);
  }
  for (int p=0;p<pl.probe_cnt();p++) {
    m_probe_vec.push_back(pl.m_probe_vec[p]);
  }
}

const ProbeListStlBlock* ProbeListStl::getBlockPtr(int b)
{
  if ((b<0)||(b>=m_block_vec.size())) {
    return NULL;
  }
  return &m_block_vec[b];
}

//
int ProbeListStl::findBlockMatch(int qAllele,int qContext,int qChannel) const
{
  for (int b=0;b<block_cnt();b++) {
    if ((qAllele==get_blockAllele(b)) &&
        (qContext==get_blockContext(b)) &&
        (qChannel==get_blockChannel(b))) {
      return b;
    }
  }
  return -1;
}

bool ProbeListStl::hasDuplicateBlocks() const
{
  for (int b=0;b<block_cnt();b++) {
    int idx=findBlockMatch(get_blockAllele(b),get_blockContext(b),get_blockChannel(b));
    if (idx!=b) {
      return true;
    }
  }
  return false;
}
    
#define PLS_DUMP "#PLS: "
// Handy method to display ourselves.
void ProbeListStl::dump() const {

  int b_cnt=block_cnt();
  int p_cnt=probe_cnt();
  printf(PLS_DUMP "%p '%s' b_cnt=%d  p_cnt=%d\n",
         this,m_name.c_str(),b_cnt,p_cnt);
  //
  if (b_cnt>0) {
    printf(PLS_DUMP "Size____ Anno____ Allele__ Context_ Channel_ RepType_\n");
    for (int b=0;b<b_cnt;b++) {
      printf(PLS_DUMP "%8d %8d %8d %8d %8d %8d\n",
             get_blockSize(b),
             get_blockAnn(b),
             get_blockAllele(b),
             get_blockContext(b),
             get_blockChannel(b),
             get_blockRepType(b));
    }
  }
  //
  if (p_cnt>0) {
    printf(PLS_DUMP "ID______ GC APID____\n");
    for (int p=0;p<p_cnt;p++) {
      printf(PLS_DUMP "%8d %2d %8d\n",get_probeId(p),get_probeGc(p),get_probeApid(p));
    }
  }
}
