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

//
#include "chipstream/ProbeListFactory.h"
//
#include "chipstream/ProbeListStl.h"
//
#include "file/CDFFileData.h"
#include "util/Err.h"
#include "util/Verbose.h"
//
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

// If you need some output
// #define PROBELIST_DEBUG 1

// If the arch cant do unaligned loads and stores we use these fuctions.
// These are in big-endian format.
#ifdef PROBELIST_NEEDS_UNALIGNED

unsigned int
probelist_load32b(unsigned int* ptr)
{
  unsigned char* cptr=(unsigned char*)ptr;
  unsigned int val;
  val=(cptr[0]<<24)+(cptr[1]<<16)+(cptr[2]<<8)+cptr[3];
  return val;
}
unsigned short
probelist_load16b(unsigned short* ptr)
{
  unsigned char* cptr=(unsigned char*)ptr;
  unsigned short val;
  val=(cptr[0]<<8)+cptr[1];
  return val;
}
void
probelist_store32b(unsigned int* ptr,unsigned int val)
{
  unsigned char* cptr=(unsigned char*)ptr;
  cptr[0]=(val>>24);
  cptr[1]=(val>>16);
  cptr[2]=(val>> 8);
  cptr[3]=(val    );
}
void
probelist_store16b(unsigned short* ptr,unsigned short val)
{
  char* cptr=(char*)ptr;
  cptr[0]=(val>> 8);
  cptr[1]=(val    );
}

#endif

///// ProbeListPacked

ProbeListPacked::ProbeListPacked()
{
  m_headptr=NULL;
}

ProbeListPacked::ProbeListPacked(ProbeList_Head* ptr)
{
  m_headptr=ptr;
}

bool
ProbeListPacked::isNull() const
{
  // return m_headptr==NULL;
  if (m_headptr==NULL) {
    return true;
  } else {
    return false;
  }
}

int
ProbeListPacked::byte_size(int b_cnt,int p_cnt,int name_len)
{
  return sizeof(ProbeList_Head)+
    b_cnt*sizeof(ProbeList_E_Block)+
    p_cnt*sizeof(ProbeList_E_Probe)+
    name_len;
}

int
ProbeListPacked::byte_size() const
{
  return sizeof(ProbeList_Head)+
    block_cnt()*sizeof(ProbeList_E_Block)+
    probe_cnt()*sizeof(ProbeList_E_Probe)+
    max_name_len();
}

//
int
ProbeListPacked::block_cnt() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD32B(&m_headptr->m_block_cnt);
}
int
ProbeListPacked::probe_cnt() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD32B(&m_headptr->m_probe_cnt);
}
int
ProbeListPacked::max_name_len() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD16B(&m_headptr->m_name_len);
}
//
int
ProbeListPacked::get_numMatch() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD32B(&m_headptr->m_numMatch);
}
void
ProbeListPacked::set_numMatch(int numMatch)
{
  assert(m_headptr!=NULL);
  PROBELIST_STORE32B(&m_headptr->m_numMatch,numMatch);
}

int
ProbeListPacked::pidxToBidx(int pidx) const
{
  // in bounds?
  int p_cnt=probe_cnt();
  if ((pidx<0)||(pidx>=p_cnt)) {
    return -1;
  }

  int b_cnt=block_cnt();
  int rem=pidx;

  for (int b=0;b<b_cnt;b++) {
    int b_size=get_blockSize(b);
    if (b_size==0) {
      continue;
    }
    // in this block?
    if (rem<b_size) {
      return b;
    }
    // not in this block...
    rem-=b_size;
  }
  // the probe is out of bounds.
  return -1;
}

//
bool
ProbeListPacked::isPm(int p_idx) const
{
  int num_match=get_numMatch();
  int p_max=probe_cnt();

  // check for sanity.
  assert((num_match==1)||(num_match==2));
  assert((0<=p_idx)&&(p_idx<p_max));
  // APT_ERR_ASSERT and especially ToStr() are big performance hits right here.

//   APT_ERR_ASSERT(((num_match==1)||(num_match==2)),"num_match must be 1 or 2. " + ToStr(num_match));
//   APT_ERR_ASSERT(((0<=p_idx)&&(p_idx<p_max)),"p_idx out of bounds.");

  // by definition, if there is one match, then it has to be a PM.
  if (num_match==1) {
    return true;
  }

  // we have PM and MM probes...
  if (num_match==2) {
    // Like pidxToBidx, but we unroll it here to keep some state.
    // step though the blocks, to find the block and then test to see
    // if it is in the first half of the block.
    int p_rem=p_idx; // probe remainder
    int b_cnt=block_cnt();
    for (int b=0;b<b_cnt;b++) {
      int b_size=get_blockSize(b);
      // no probes in a zero sized block
      if (b_size==0) {
        continue;
      }
      // is the probe index in this block?
      if (p_rem<b_size) {
        // in the first half is a PM;
        if (p_rem<(b_size/2)) {
          return true;
        }
        else {
          return false;
        }
      }
      // not in this block, deduct the size...
      p_rem-=b_size;
    }
  }

  APT_ERR_ABORT("Shouldnt get here.");
  return false;
}

int
ProbeListPacked::get_type() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD32B(&m_headptr->m_type);
}

void
ProbeListPacked::set_type(int type)
{
  assert(m_headptr!=NULL);
  PROBELIST_STORE32B(&m_headptr->m_type,type);
}
int
ProbeListPacked::getApidStart() const
{
  assert(m_headptr!=NULL);
  return PROBELIST_LOAD32B(&m_headptr->m_apid_start);
}
void
ProbeListPacked::setApidStart(int start)
{
  assert(m_headptr!=NULL);
  PROBELIST_STORE32B(&m_headptr->m_apid_start,start);
}

// Block getters & setters
ProbeList_E_Block*
ProbeListPacked::get_addr_E_block(unsigned int idx) const
{
  assert(m_headptr!=NULL);
  assert(idx<m_headptr->m_block_cnt);
  //
  char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+
    (idx*sizeof(ProbeList_E_Block));
  return (ProbeList_E_Block*)ptr;
}

int
ProbeListPacked::get_blockSize(int i) const
{
  return PROBELIST_LOAD32B(&get_addr_E_block(i)->m_blockSize);
}
void
ProbeListPacked::set_blockSize(int i,int val)
{
  PROBELIST_STORE32B(&get_addr_E_block(i)->m_blockSize,val);
}
//
short
ProbeListPacked::get_blockAnn(int i) const
{
  return PROBELIST_LOAD16B(&get_addr_E_block(i)->m_blockAnn);
}
void
ProbeListPacked::set_blockAnn(int i,short val)
{
  PROBELIST_STORE16B(&get_addr_E_block(i)->m_blockAnn,val);
}
//
int
ProbeListPacked::get_blockAllele(int i) const
{
  return PROBELIST_LOAD32B(&get_addr_E_block(i)->m_blockAllele);
}
void
ProbeListPacked::set_blockAllele(int i,int val)
{
  PROBELIST_STORE32B(&get_addr_E_block(i)->m_blockAllele,val);
}
//
int
ProbeListPacked::get_blockContext(int i) const
{
  return PROBELIST_LOAD32B(&get_addr_E_block(i)->m_blockContext);
}
void
ProbeListPacked::set_blockContext(int i,int val)
{
  PROBELIST_STORE32B(&get_addr_E_block(i)->m_blockContext,val);
}
//
int
ProbeListPacked::get_blockChannel(int i) const
{
  return PROBELIST_LOAD8B(&get_addr_E_block(i)->m_blockChannel);
}
void
ProbeListPacked::set_blockChannel(int i,int val)
{
  // test all the channels are "<0" to check sign extention.
  // APT_ERR_ASSERT(val<0,"channel");
  PROBELIST_STORE8B(&get_addr_E_block(i)->m_blockChannel,val);
}
//
int
ProbeListPacked::get_blockRepType(int i) const
{
  return PROBELIST_LOAD8B(&get_addr_E_block(i)->m_blockRepType);
}
void
ProbeListPacked::set_blockRepType(int i,int val)
{
  PROBELIST_STORE8B(&get_addr_E_block(i)->m_blockRepType,val);
}

/////

/// @todo this should be somehow unified with the ProbeListStl version.
int ProbeListPacked::findBlockMatch(int qAllele,int qContext,int qChannel) const
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

bool ProbeListPacked::hasDuplicateBlocks() const
{
  for (int b=0;b<block_cnt();b++) {
    int idx=findBlockMatch(get_blockAllele(b),get_blockContext(b),get_blockChannel(b));
    if (idx!=b) {
      return true;
    }
  }
  return false;
}

/////

/*
// Probe getters & setters
ProbeList_E_Probe* const
ProbeListPacked::get_addr_E_probe(unsigned int idx)
{
  APT_ERR_ASSERT(m_headptr==NULL,
                 "Head pointer is NULL.");
  APT_ERR_ASSERT(idx<probe_cnt(),
                 "Probeset: '" + get_name_string() + "' - wrong number of probes."
                 "(Max: " + ToStr(probe_cnt() - 1) + " Got: " + ToStr(idx)_")");

  char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+
    (block_cnt()*sizeof(ProbeList_E_Block))+
    (idx*sizeof(ProbeList_E_Probe));
  return (ProbeList_E_Probe*)ptr;
}
*/

// Probe getters & setters
ProbeList_E_Probe*
ProbeListPacked::get_addr_E_probe(unsigned int idx) const
{
  assert(m_headptr!=NULL);
  assert(idx<probe_cnt());
  // APT_ERR_ASSERT, get_name_string(), and especially ToStr() are big
  // performance hits right here.

//   APT_ERR_ASSERT(m_headptr!=NULL,
//                  "Head pointer is NULL.");
//   APT_ERR_ASSERT(idx<probe_cnt(),
//                  "Probeset: '" + get_name_string() + "' - wrong number of probes."
//                  "(Max: " + ToStr(probe_cnt() - 1) + " Got: " + ToStr(idx)+")");

  char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+
    (block_cnt()*sizeof(ProbeList_E_Block))+
    (idx*sizeof(ProbeList_E_Probe));
  return (ProbeList_E_Probe*)ptr;
}

int
ProbeListPacked::get_probeId(int i) const
{
  return PROBELIST_LOAD32B(&get_addr_E_probe(i)->m_probeId);
}
void
ProbeListPacked::set_probeId(int i,int val)
{
  PROBELIST_STORE32B(&get_addr_E_probe(i)->m_probeId,val);
}
void
ProbeListPacked::get_probeIds(std::vector<int>& probeids) const
{
  int p_cnt=probe_cnt();
  probeids.resize(p_cnt);

  if (p_cnt==0) {
    return;
  }

  // get the address of the first probe.
  ProbeList_E_Probe* p_ptr=get_addr_E_probe(0);
  for (int i=0;i<p_cnt;i++) {
    probeids[i]=p_ptr->m_probeId;
    p_ptr++;
  }
}

int
ProbeListPacked::get_probeIdxForBlock(int blockIdx) const
{
  APT_ERR_ASSERT(blockIdx<block_cnt(),"blockIdx out of range.");

  int probeIdx=0;
  for (int bi=0;bi<blockIdx;bi++) {
    probeIdx+=get_blockSize(bi);
  }
  return probeIdx;
}

void
ProbeListPacked::get_probeIdsForBlock(int blockIdx,std::vector<int>& probeids) const
{
  int p_cnt=get_blockSize(blockIdx);
  probeids.resize(p_cnt);

  if (p_cnt==0) {
    return;
  }

  // get the address of the first probe.
  int p_idx=get_probeIdxForBlock(blockIdx);
  ProbeList_E_Probe* p_ptr=get_addr_E_probe(p_idx);
  for (int i=0;i<p_cnt;i++) {
    probeids[i]=p_ptr->m_probeId;
    p_ptr++;
  }
}

int
ProbeListPacked::get_probeApid(int i) const
{
  if ((i>=0)&&(i<m_headptr->m_probe_cnt)) {
    return getApidStart()+i;
  }
  return -1;
}
int8_t
ProbeListPacked::get_probeGc(int i) const
{
  return PROBELIST_LOAD8B(&get_addr_E_probe(i)->m_probeGc);
}
void
ProbeListPacked::set_probeGc(int i,int8_t val)
{
  PROBELIST_STORE8B(&get_addr_E_probe(i)->m_probeGc,val);
}

// Name getters & setters;
char*
ProbeListPacked::get_addr_E_name()
{
  assert(m_headptr!=NULL);
  char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+
    (block_cnt()*sizeof(ProbeList_E_Block))+
    (probe_cnt()*sizeof(ProbeList_E_Probe));
  return ptr;
}

// Name getters & setters;
const char*
ProbeListPacked::get_addr_E_name() const
{
  assert(m_headptr!=NULL);
  char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+
    (block_cnt()*sizeof(ProbeList_E_Block))+
    (probe_cnt()*sizeof(ProbeList_E_Probe));
  return ptr;
}

const char*
ProbeListPacked::get_name_cstr() const
{
  assert(m_headptr!=NULL);
  int name_len=PROBELIST_LOAD16B(&m_headptr->m_name_len);
  if (name_len==0) {
    return NULL;
  }
  return get_addr_E_name();
}

std::string
ProbeListPacked::get_name_string() const
{
  std::string name;
  name=get_addr_E_name();
  return name;
}

void
ProbeListPacked::set_name(const std::string &name)
{
  assert(m_headptr!=NULL);
  int   name_len=PROBELIST_LOAD16B(&m_headptr->m_name_len);

  if (name_len>0) { // there is space
    char* name_ptr=get_addr_E_name();
    if (name_len>1) { // more than just a null?
      // trunc the name to fit.
      strncpy(name_ptr,name.c_str(),name_len-1);
    }
    // strncpy might not put the final '0'
    name_ptr[name_len-1]=0;
  }
}

int
ProbeListPacked::findProbeApid(int qPid,int qAllele,int qContext,int qChannel) const
{
  int pidx_start;
  int pidx_end;
  // "pidx" is the internal offset in this ProbeListPacked
  pidx_start=0;
  for (int b=0;b<m_headptr->m_block_cnt;b++) {
    // does this block match?
    if ((get_blockAllele(b)==qAllele) &&
        (get_blockContext(b)==qContext) &&
        (get_blockChannel(b)==qChannel))
    {
      // hit on block, se
      pidx_end=pidx_start+get_blockSize(b);
      for (int pidx=pidx_start;pidx<pidx_end;pidx++) {
        // should be renamed to "get_probeid(pidx)"
        if (get_probeId(pidx)==qPid) {
          return getApidStart()+pidx;
        }
      }
      // not found in this block
      pidx_start+=get_blockSize(b);
    }
    else {
      pidx_start+=get_blockSize(b);
    }
  }
  // not found
  return -1;
}


#define PLP_DUMP "#PLP: "
void
ProbeListPacked::dump() const
{
  if (isNull()) {
    printf(PLP_DUMP " %p: isNull",this);
    return;
  }

  int b_cnt=block_cnt();
  int p_cnt=probe_cnt();
  int num_match=get_numMatch();

  printf(PLP_DUMP "%p: '%s'  b_cnt=%d  p_cnt=%d  num_match=%d\n",
         this,get_name_cstr(),b_cnt,p_cnt,num_match);

  if (b_cnt>0) {
    printf(PLP_DUMP "Size____ Anno____ Allele_ Context_ Channel_ RepType_\n");
    for (int b=0;b<b_cnt;b++) {
      printf(PLP_DUMP "%8d %8d %8d %8d %8d %8d\n",
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
    printf(PLP_DUMP "ID______ GC APID____ ____\n");
    for (int p=0;p<p_cnt;p++) {
      printf(PLP_DUMP "%8d %2d %8d %4s\n",
             get_probeId(p),get_probeGc(p),get_probeApid(p),
             (isPm(p)?"PM":"  mm"));
    }
  }
}


//////////////////////////////

///// ProbeListFactory

// *poof!* A new heap is born!
ProbeListFactory::ProbeListFactory() : m_numCols(0), m_numRows(0), m_numChannels(1)
{
  m_probeid2apid_vecptr=NULL;
  //
  delete_regions(); // also clears counters
  set_default_region_size(1024*1024); // default size
}
//
ProbeListFactory::~ProbeListFactory()
{
  delete_regions();
}

//
void
ProbeListFactory::set_default_region_size(int size)
{
  m_default_region_size=size;
}

//
void
ProbeListFactory::new_region()
{
  new_region(m_default_region_size);
}

void
ProbeListFactory::new_region(int region_size)
{
#ifdef PROBELIST_DEBUG
  printf("alloc_region(%d,%d)\n",m_ridx,region_size);
#endif
  ProbeListFactory_Region new_region;
  m_region.push_back(new_region);
  //
  m_region[m_ridx].m_start_ptr=(char*)malloc(region_size);
  assert(m_region[m_ridx].m_start_ptr!=NULL);
  memset(m_region[m_ridx].m_start_ptr,0,region_size);
  m_region[m_ridx].m_fill_ptr=m_region[m_ridx].m_start_ptr;
  m_region[m_ridx].m_end_ptr=m_region[m_ridx].m_start_ptr+region_size;
}

// Get rid of all our memory.
void
ProbeListFactory::delete_regions()
{
  // clear the memory.
  clear();

  // free the memory...
  for (unsigned int i=0;i<m_region.size();i++) {
    free(m_region[i].m_start_ptr);
  }
  // ...and reset the vec
  m_region.clear();
}

// zero out our info.
void
ProbeListFactory::clear()
{
#ifdef PROBELIST_DEBUG
  printf("heap.clear()\n");
#endif
  //
  m_ridx=0;
  m_allocs=0;
  m_probecnt_total=0;
  //
  m_probelist_vec.clear();
  //
  m_name2idx_vec.clear();
  m_name2idx_issorted=false;
  m_name2idx_sortcnt=0;
  //
  for (size_t i=0;i<m_region.size();i++) {
    m_region[i].m_fill_ptr=m_region[i].m_start_ptr;
    // zero out the region for future use.
    unsigned int len=m_region[i].m_end_ptr-m_region[i].m_start_ptr;
    memset(m_region[i].m_start_ptr,0,len);
  }

  delete m_probeid2apid_vecptr;
  m_probeid2apid_vecptr=NULL;
}

//
char*
ProbeListFactory::alloc_ProbeList_bytes(int size)
{
  m_allocs++;
  // Work along the regions until we find space
  while (1) {
    //printf("alloc_bytes(%d)\n",size);
    // have an active region?
    if (m_region.size()<=m_ridx) {
      // make it at least big enough for this allocation
      if (m_default_region_size<size) {
        new_region(size);
      }
      else {
        new_region(m_default_region_size);
      }
    }

    char* new_fill_ptr;
    char* list_ptr;

    // able to alloc from the current region?
    list_ptr=m_region[m_ridx].m_fill_ptr;
    new_fill_ptr=list_ptr+size;
    if (new_fill_ptr<=m_region[m_ridx].m_end_ptr) {
      m_region[m_ridx].m_fill_ptr=new_fill_ptr;
      // for debugging, wack the region just allocated with
      // an easy to spot value
      // memset(list_ptr,-1,size);
      return list_ptr;
    }

    // Nope! move on to the next region.
    m_ridx++;
  }
  // not reached && keep gcc happy
  assert(0);
  return NULL;
}

//
void ProbeListFactory::reserve(int size)
{
  m_probelist_vec.reserve(size);
  m_name2idx_vec.reserve(size);
}

void
ProbeListFactory::name2idx_dirty()
{
  m_name2idx_issorted=false;
}

//////////

/// Finish adding the packed probelist to the factory
void
ProbeListFactory::add_probelist_finish(ProbeListPacked& pl)
{
  // cant add a null probelist
  APT_ERR_ASSERT(!pl.isNull(),"pl.isNull()");
  // set the start of the list and update the total...
  pl.setApidStart(m_probecnt_total);
  m_probecnt_total+=pl.probe_cnt();
  // also add the index to the name2idx vec...
  m_name2idx_vec.push_back(m_probelist_vec.size());
  // add it
  m_probelist_vec.push_back(pl);
  // ...which is no longer sorted.
  name2idx_dirty();
}

//
ProbeListPacked
ProbeListFactory::add_ProbeList(int block_cnt,int probe_cnt,int name_len)
{
  ProbeListPacked pl;
  //
  /// @todo change all asserts to APT_ERR_ASSERT()
  assert((unsigned int)block_cnt<PROBELIST_BLOCK_CNT_MAX);
  assert((unsigned int)probe_cnt<PROBELIST_PROBE_CNT_MAX);
  assert((unsigned int)name_len<PROBELIST_NAME_LEN_MAX);
  //
  int required_size=pl.byte_size(block_cnt,probe_cnt,name_len);
  pl.m_headptr=(ProbeList_Head*)alloc_ProbeList_bytes(required_size);
  assert(pl.m_headptr!=NULL);
  //
  PROBELIST_STORE32B(&pl.m_headptr->m_block_cnt,block_cnt);
  PROBELIST_STORE32B(&pl.m_headptr->m_probe_cnt,probe_cnt);
  PROBELIST_STORE16B(&pl.m_headptr->m_name_len,name_len);
  //
  add_probelist_finish(pl);
  return pl;
};

//
ProbeListPacked
ProbeListFactory::add_ProbeList(int block_cnt,int probe_cnt, const std::string &name)
{
  // +1 for the null
  ProbeListPacked pl=add_ProbeList(block_cnt,probe_cnt,name.size()+1);
  pl.set_name(name);
  //
  return pl;
}

// Be sure to keep ProbeListStl::opreator= updated.
ProbeListPacked
ProbeListFactory::add_ProbeList(const ProbeListStl& pls)
{
  // Refuse to add probelists which have duplicates.
  APT_ERR_ASSERT(pls.hasDuplicateBlocks()==false,"internal error");

  //
  ProbeListPacked plp;

  // this adds it to m_probelist_vec so we dont.
  plp=add_ProbeList(pls.block_cnt(),pls.probe_cnt(),pls.get_name_string());
  //
  plp.set_type(pls.get_type());
  plp.set_numMatch(pls.get_numMatch());
  //
  for (int i=0;i<pls.block_cnt();i++) {
    plp.set_blockSize(i,pls.get_blockSize(i));
    plp.set_blockAllele(i,pls.get_blockAllele(i));
    plp.set_blockContext(i,pls.get_blockContext(i));
    plp.set_blockChannel(i,pls.get_blockChannel(i));
    plp.set_blockRepType(i,pls.get_blockRepType(i));
  }
  for (int i=0;i<pls.probe_cnt();i++) {
    plp.set_probeId(i,pls.get_probeId(i));
    plp.set_probeGc(i,pls.get_probeGc(i));
    // we dont store the universal probe as ProbeListFactory computes that.
  }
  // compare the two items by eye.
  //pls.dump();
  //plp.dump();
  return plp;
}

double ProbeListFactory::safe_div(double a,double b)
{
  if (b==0.0) {
    return 0.0;
  }
  return (a/b);
}

void
ProbeListFactory::dump_regions() const
{
  unsigned int space=0;
  unsigned int used=0;
  unsigned int total_space=0;
  unsigned int total_used=0;
  double percent;
  double avg_size;

  printf("== ProbeListFactory Usage: (this=%p) allocs=%u\n",this,m_allocs);
  printf("== Reg | Ptr_____________ | Size____ | Used____ | %%____ |\n");
  for (unsigned int r=0;r<m_region.size();r++) {
    space=m_region[r].m_end_ptr -m_region[r].m_start_ptr;
    used =m_region[r].m_fill_ptr-m_region[r].m_start_ptr;
    total_space+=space;
    total_used+=used;
    percent=safe_div(used,space)*100.0;
    printf("== %3d | %16p | %8u | %8u | %5.2f |\n",
           r,m_region[r].m_start_ptr,space,used,percent);
  }
  percent=safe_div(total_used,total_space)*100.0;
  printf("== Total:                 | %8u | %8u | %5.2f |\n",total_space,total_used,percent);
  avg_size=safe_div(total_used,m_allocs);
  printf("== avg size=%6.1f\n",avg_size);

}

// @todo change back
// ProbeListPacked
void
ProbeListFactory::pListFromPSet(const ProbeSet &ps) {
  unsigned char nGroups = max((int)ps.numGroups,1);
  bool hasMM = hasNonPm(ps);
  string name;
  int pmCount = countPmProbes(&ps);
  int probeCount = pmCount;

  if(hasMM) {
    probeCount += pmCount;
  }

  // are we using a name or an id?
  if(ps.name != NULL && !Util::sameString(ps.name, ""))
    name = ps.name;
  else {
    Err::errAbort("ProbeListFactory::pListFromPSet() - No name for probeset.");
  }
  //
  ProbeListPacked pList = add_ProbeList(nGroups, probeCount, name );
  // let people know we have mismatches.
  if(hasMM)
    pList.set_numMatch(2);
  else
    pList.set_numMatch(1);
  //  pList.set_block(0,0);
  pList.set_type(ps.psType);
  // fill the probes with null_probe indicator to trip ourselves up if we try to use them
  for(int i = 0; i < probeCount; i++)
    pList.set_probeId(i, ProbeList::null_probe);
  int pCount = 0;
  vector<int> atomProbes(ps.atoms.size(), 0);
  // fill in the probes.
  for(uint32_t atomIx = 0; atomIx < ps.atoms.size(); atomIx++) {
    Atom *a = ps.atoms[atomIx];
    for(uint32_t probeIx = 0; probeIx < a->probes.size(); ++probeIx) {
      Probe *current = a->probes[probeIx];
      Probe *next = NULL;
      atomProbes[atomIx]++;
      if(probeIx +1 < a->probes.size())
        next = a->probes[probeIx + 1];
      // if it is an pm-mm set then put pm in front and mm in matching position.
      if(Probe::isPm(*current) && next != NULL && Probe::isMm(*next)) {
        pList.set_probeId(pCount, current->id);
        pList.set_probeGc(pCount, current->gcCount);
        pList.set_probeId(pCount + pmCount, next->id);
        pList.set_probeGc(pCount + pmCount, next->gcCount);
        probeIx++; // advance probeIx count as we did two.
      }
      // if it is a lone mm then put a "null_probe" in for the pm probe.
      if(!Probe::isPm(*current)) {
        pList.set_probeId(pCount, ProbeList::null_probe);
        pList.set_probeGc(pCount, -1);
        pList.set_probeId(pCount + pmCount, current->id);
        pList.set_probeGc(pCount + pmCount, current->gcCount);
      }
      // lone pm probe
      else {
        pList.set_probeId(pCount, current->id);
        pList.set_probeGc(pCount, current->gcCount);
      }
      pCount++;
    }
  }

  // the first atom is taken as to have allele and context for its block.
  int atomIdx=0;
  for(int groupIx=0; groupIx < ps.numGroups; groupIx++) {
    Atom* a=ps.atoms[atomIdx];
    atomIdx+=ps.atomsPerGroup[groupIx];
    //
    pList.set_blockAllele(groupIx,a->allele_code);
    pList.set_blockContext(groupIx,a->context_code);
    pList.set_blockChannel(groupIx,a->channel_code);
    pList.set_blockRepType(groupIx,a->replication_type);
  }

  /* Setup the block offsets into the probe array. */
  int last = 0;
  int current = 0;
  for(int groupIx = 0; groupIx < ps.numGroups; groupIx++) {
    current += ps.atomsPerGroup[groupIx];
    int probeAtomCount=0;
    for(; last < current; last++) {
      probeAtomCount += atomProbes[last];
    }
    pList.set_blockSize(groupIx, probeAtomCount);
    pList.set_blockAnn(groupIx, ProbeList::NoAnnotation);
  }
  if(ps.numGroups == 4 &&
     (ps.psType == ProbeSet::GenoType)) {
    pList.set_blockAnn(0, ProbeList::A_AlleleForward);
    pList.set_blockAnn(1, ProbeList::B_AlleleForward);
    pList.set_blockAnn(2, ProbeList::A_AlleleReverse);
    pList.set_blockAnn(3, ProbeList::B_AlleleReverse);
  }
  else if(ps.numGroups == 2 &&
          (ps.psType == ProbeSet::GenoType)) {
    pList.set_blockAnn(0, ProbeList::A_Allele);
    pList.set_blockAnn(1, ProbeList::B_Allele);
  }

  //return pList;
}

bool ProbeListFactory::maybeInsertProbe(Atom *a, int position, int id, unsigned char gcCount, unsigned char type,int apid) {
  if (id != ProbeList::null_probe) {
    Probe *p = new Probe();
    p->id = id;
    p->type = type;
    p->gcCount = gcCount;
    p->m_apid=apid;
    //
    a->probes[position] = p;
    return true;
  }
  return false;
}

/// @todo this should be "ProbeSet::operator=()"
ProbeSet *ProbeListFactory::asProbeSet(const ProbeListPacked& pl) {
  ProbeSet *ps = new ProbeSet();
  ps->name = Util::cloneString(pl.get_name_cstr());
  ps->psType = (ProbeSet::ProbeSetType)pl.get_type();
  ps->numGroups = pl.block_cnt();
  ps->atomsPerGroup.resize(ps->numGroups, 0);
  bool hasMM = pl.get_numMatch() == 2;
  int pmCount = pl.probe_cnt() / pl.get_numMatch();

  ps->atoms.resize(pl.block_cnt()); // make an atom for each block...
  int currentProbeStart = 0;
  for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
    Atom *a = new Atom();
    ps->atoms[atomIx] = a;

    a->id = atomIx;
    a->allele_code=pl.get_blockAllele(atomIx);
    a->context_code=pl.get_blockContext(atomIx);
    a->channel_code=pl.get_blockChannel(atomIx);
    a->replication_type=(Atom::ReplicationType)pl.get_blockRepType(atomIx);

    int notNullCount = 0;
    int atomProbeEnd = pl.get_blockSize(atomIx) + currentProbeStart;

    for(int i = currentProbeStart; i < atomProbeEnd; i++) {
      if(pl.get_probeId(i) != ProbeList::null_probe) {
        notNullCount++;
      }
      if(hasMM && pl.get_probeId(i+pmCount) != ProbeList::null_probe) {
        notNullCount++;
      }
    }
    a->probes.resize(notNullCount);

    ps->atomsPerGroup[atomIx] = 1; // can do this as there is just one atom per block.
    int probesInserted = 0;
    for(int i = currentProbeStart; i < atomProbeEnd; i++) {
      int probeIndex = i;
      int pId = pl.get_probeId(probeIndex);
      // regular pm-mm pair
      ///@todo I do not think we keep track of st vs at probes -- so we should really
      ///      have a Probe::PM and Probe::MM and set that here
      if (pId != ProbeList::null_probe && hasMM) {
        if (maybeInsertProbe(a, probesInserted, pId, pl.get_probeGc(probeIndex), Probe::PMST,pl.get_probeApid(probeIndex)))
          probesInserted++;
        int mmId = pl.get_probeId(i+pmCount);
        if (maybeInsertProbe(a, probesInserted, mmId, pl.get_probeGc(probeIndex+pmCount), Probe::MMST,pl.get_probeApid(probeIndex+pmCount)))
          probesInserted++;
      }
      else if (pId == ProbeList::null_probe && hasMM) {
        int mmId = pl.get_probeId(i+pmCount);
        if (maybeInsertProbe(a, probesInserted, mmId, pl.get_probeGc(probeIndex+pmCount), Probe::MMST,pl.get_probeApid(probeIndex+pmCount)))
          probesInserted++;
      }
      else if (pId != ProbeList::null_probe) {
        if (maybeInsertProbe(a, probesInserted, pId, pl.get_probeGc(probeIndex), Probe::PMST,pl.get_probeApid(probeIndex)))
          probesInserted++;
      }
      else {
        Err::errAbort("Null probe with no mismatch in probeset: " + ToStr(ps->name));
      }
      if(probesInserted > a->probes.size()) {
        Err::errAbort("ProbeListFactory::asProbeSet() - " + ToStr(ps->name) + "More probes than space for.");
      }
    }
    currentProbeStart += atomProbeEnd - currentProbeStart;
  }
  return ps;
}

bool ProbeListFactory::hasNonPm(const ProbeSet &ps) {
  for(uint32_t atomIx = 0; atomIx < ps.atoms.size(); ++atomIx) {
    Atom *a = ps.atoms[atomIx];
    for(uint32_t probeIx = 0; probeIx < a->probes.size(); ++probeIx) {
      if(!Probe::isPm(*(a->probes[probeIx]))) {
        return true;
      }
    }
  }
  return false;
}

unsigned int ProbeListFactory::countPmProbes(const ProbeSet *ps) {
  unsigned int atomPmCount = 0;
  for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
    Atom &atom = *(ps->atoms[atomIx]);
    for(int probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
      Probe *p = atom.probes[probeIx];
      Probe *next = NULL;
      if(probeIx + 1 < atom.probes.size()) {
        next = atom.probes[probeIx + 1];
      }
      if(p->type != Probe::MMST && p->type != Probe::MMAT) {
        atomPmCount++;
        if(next != NULL && Probe::isMm(*next)) {
          probeIx++;
        }
      }
      else if(!Probe::isPm(*p)) {
        atomPmCount++;
      }
    }
  }
  return atomPmCount;
}

//////////

void
ProbeListFactory::dump() const
{
  dump_regions();
  dump_probelists();
}

void
ProbeListFactory::dump_probelist(int i) const
{
  if ((0<=i)||(i<m_probelist_vec.size())) {
    m_probelist_vec[i].dump();
  }
  else {
    printf("# %3d : out of bounds. (0-%d)\n",i,int(m_probelist_vec.size()));
  }
}

void
ProbeListFactory::dump_probelists() const
{
  for (int i=0;i<m_probelist_vec.size();i++) {
    printf(" %3d :",i);
    m_probelist_vec[i].dump();
  }
}

//////////

void read_spf_v2(ProbeListFactory& plf,affx::SpfFile& spfFile)
{
  ProbeListStl pl;
  //
  std::string i_name;
  int i_type;
  int i_num_blocks;
  int i_num_probes;
  // @todo this should be used.
  // int i_block_sizes;
  int i_num_match;
  std::vector<int> i_vector;
  //int pi;

  //
  while (spfFile.nextLevel(0)==affx::TSV_OK) {
    //
    pl.clear();
    //
    spfFile.get(0,spfFile.m_name_cidx,i_name);
    pl.set_name(i_name);
    //
    //printf("%s\n",pl.get_name_cstr());
    //
    i_type=-1;
    spfFile.get(0,spfFile.m_type_cidx,i_type);
    pl.set_type(i_type);
    //
    i_num_blocks=-1;
    i_num_probes=-1;
    spfFile.get(0,spfFile.m_v2_num_blocks_cidx,i_num_blocks);
    spfFile.get(0,spfFile.m_v2_num_probes_cidx,i_num_probes);
    //
    pl.resize(0,0);
    pl.resize(i_num_blocks,i_num_probes);
    // @todo: for now we assume and check all the blocks are the same size.
    spfFile.get(0,spfFile.m_v2_block_sizes_cidx,&i_vector);
    //for (int i=0;i<i_vector.size();i++) {
    //  APT_ERR_ASSERT(i_vector[0]==i_vector[i],"block size mismatch. dataline="+ToStr(spfFile.lineNumber()));
    //}
    for (int i=0;i<i_num_blocks;i++) {
      pl.set_blockSize(i,i_vector[i]);
    }
    //
    i_num_match=-1;
    if (spfFile.m_num_match_cidx>=0) {
      spfFile.get(0,spfFile.m_num_match_cidx,i_num_match);
    }
    else {
      i_num_match=1;
    }
    pl.set_numMatch(i_num_match);
    //
    spfFile.get(0,spfFile.m_v2_block_annotations_cidx,&i_vector);
    for (int i=0;i<i_num_blocks;i++) {
      pl.set_blockAnn(i,i_vector[i]);
    }
    // required annotations
    spfFile.get(0,spfFile.m_v2_block_annotations_cidx,&i_vector);
    for (int i=0;i<i_num_blocks;i++) {
      pl.set_blockAnn(i,i_vector[i]);
    }
    // required probes
    spfFile.get(0,spfFile.m_v2_probes_cidx,&i_vector);
    APT_ERR_ASSERT((i_vector.size()==i_num_probes),"probe count mismatch");
    for (int i=0;i<i_vector.size();i++) {
      pl.set_probeId(i,TO_ZEROBASED(i_vector[i]));
      pl.set_probeApid(i,-1);
    }
    // optional columns...
    //
    if (spfFile.m_v2_block_alleles_cidx>=0) {
      spfFile.get(0,spfFile.m_v2_block_alleles_cidx,&i_vector);
      for (int i=0;i<i_vector.size();i++) {
        pl.set_blockAllele(i,i_vector[i]);
      }
    }
    if (spfFile.m_v2_block_contexts_cidx>=0) {
      spfFile.get(0,spfFile.m_v2_block_contexts_cidx,&i_vector);
      for (int i=0;i<i_vector.size();i++) {
        pl.set_blockContext(i,i_vector[i]);
      }
    }
    if (spfFile.m_v2_block_channels_cidx>=0) {
      spfFile.get(0,spfFile.m_v2_block_channels_cidx,&i_vector);
      for (int i=0;i<i_vector.size();i++) {
        pl.set_blockChannel(i,i_vector[i]);
      }
    }
    if (spfFile.m_v2_block_rep_types_cidx>=0) {
      spfFile.get(0,spfFile.m_v2_block_rep_types_cidx,&i_vector);
      for (int i=0;i<i_vector.size();i++) {
        pl.set_blockRepType(i,i_vector[i]);
      }
    }
    //
    plf.add_ProbeList(pl);
  } // while
}

void read_spf_v3(ProbeListFactory& plf,affx::SpfFile& spfFile)
{
  ProbeListStl pl;
  //
  std::string i_name;
  int i_type;
  int i_num_match;
  int i_annotation_code;
  int i_allele_code;
  int i_context_code;
  int i_channel_code;
  int i_rep_type = 0;
  int i_probe_id;

  //
  while (spfFile.nextLevel(0)==affx::TSV_OK) {
    pl.clear();
    //
    spfFile.get(0,spfFile.m_name_cidx,i_name);
    pl.set_name(i_name);
    //
    i_type=-1;
    spfFile.get(0,spfFile.m_type_cidx,i_type);
    pl.set_type(i_type);
    //
    i_num_match=-1;
    if (spfFile.m_num_match_cidx>=0) {
      spfFile.get(0,spfFile.m_num_match_cidx,i_num_match);
    }
    else {
      i_num_match=1;
    }
    pl.set_numMatch(i_num_match);
    //
    while (spfFile.nextLevel(1)==affx::TSV_OK) {
      spfFile.get(1, spfFile.m_v3_allele_code_cidx, i_allele_code);
      //
      while (spfFile.nextLevel(2)==affx::TSV_OK) {
        spfFile.get(2,spfFile.m_v3_annotation_code_cidx,i_annotation_code);
        spfFile.get(2,spfFile.m_v3_context_code_cidx,i_context_code);
        spfFile.get(2,spfFile.m_v3_channel_code_cidx,i_channel_code);
        //
        pl.push_block(i_annotation_code,i_allele_code,i_context_code,i_channel_code,i_rep_type);
        //
        while (spfFile.nextLevel(3)==affx::TSV_OK) {
          spfFile.get(3,spfFile.m_v3_probe_id_cidx,i_probe_id);
          pl.push_probe(TO_ZEROBASED(i_probe_id));
        }
      }
    }
    //
    pl.set_numMatch(pl.probe_cnt()/pl.block_cnt());
    //pl.dump();
    plf.add_ProbeList(pl);
  } // while
}

void read_spf_v4(ProbeListFactory& plf,affx::SpfFile& spfFile)
{
  ProbeListStl pl;
  //
  std::string i_name;
  int i_type;
  int i_num_match;
  int i_annotation_code;
  int i_allele_code;
  int i_context_code;
  int i_channel_code;
  int i_rep_type;
  int i_probe_id;

  //
  while (spfFile.nextLevel(0)==affx::TSV_OK) {
    pl.clear();
    //
    spfFile.get(0,spfFile.m_name_cidx,i_name);
    pl.set_name(i_name);
    //
    i_type=-1;
    spfFile.get(0,spfFile.m_type_cidx,i_type);
    pl.set_type(i_type);
    //
    i_num_match=-1;
    if (spfFile.m_num_match_cidx>=0) {
      spfFile.get(0,spfFile.m_num_match_cidx,i_num_match);
    }
    else {
      i_num_match=1;
    }
    pl.set_numMatch(i_num_match);
    //
    while (spfFile.nextLevel(1)==affx::TSV_OK) {
      spfFile.get(1,spfFile.m_v4_annotation_code_cidx,i_annotation_code);
      spfFile.get(1,spfFile.m_v4_allele_code_cidx,i_allele_code);
      spfFile.get(1,spfFile.m_v4_context_code_cidx,i_context_code);
      spfFile.get(1,spfFile.m_v4_channel_code_cidx,i_channel_code);
      spfFile.get(1,spfFile.m_v4_rep_type_cidx,i_rep_type);
      //
      pl.push_block(i_annotation_code,i_allele_code,i_context_code,i_channel_code,i_rep_type);
      //
      while (spfFile.nextLevel(2)==affx::TSV_OK) {
        spfFile.get(2,spfFile.m_v4_probe_id_cidx,i_probe_id);
        pl.push_probe(TO_ZEROBASED(i_probe_id));
      }
    }
    pl.set_numMatch(pl.probe_cnt()/pl.block_cnt());
    //pl.dump();
    plf.add_ProbeList(pl);
  } // while
}

void read_spf(ProbeListFactory& plf,affx::SpfFile& spfFile)
{
  // try and help out if this is unset.
  if (spfFile.m_spf_format==0) {
    spfFile.determineFormat();
  }
  //
  if (spfFile.m_spf_format==2) {
    read_spf_v2(plf,spfFile);
  }
  else if (spfFile.m_spf_format==3) {
    read_spf_v3(plf,spfFile);
  }
  else if (spfFile.m_spf_format==4) {
    read_spf_v4(plf,spfFile);
  }
  else {
    APT_ERR_ABORT("Unable to determine SpfFile format.");
  }
}

void ProbeListFactory::readSpfFile(const std::string& spfName)
{
  affx::SpfFile spfFile;
  spfFile.openSpf(spfName);
  // copy over some info.
  m_numCols=spfFile.m_header_numCols;
  m_numRows=spfFile.m_header_numRows;
  m_numChannels=spfFile.m_header_numChannels;
  m_chipTypes=spfFile.m_header_chipTypes;
  // m_header_num_probesets is there, but we dont care about it.
  //
  readSpfFile(spfFile);
  spfFile.close();
}
void ProbeListFactory::readSpfFile(affx::SpfFile& spfFile)
{
  read_spf(*this,spfFile);
}

//////////

void ProbeListFactory::writeToSpfFile_v2(affx::SpfFile& spfFile) const
{
  ProbeListPacked pl;
  std::vector<int> v_block_sizes;
  std::vector<int> v_block_annotations;
  std::vector<int> v_block_alleles;
  std::vector<int> v_block_contexts;
  std::vector<int> v_block_channels;
  std::vector<int> v_block_rep_types;
  std::vector<int> v_probes;

  // name type num_blocks block_sizes block_annotations num_match num_probes probes
  // block_alleles block_contexts block_channels

  for (int pl_i=0;pl_i<m_probelist_vec.size();pl_i++) {
    pl=m_probelist_vec[pl_i];
    //
    v_block_sizes.clear();
    v_block_annotations.clear();
    v_block_alleles.clear();
    v_block_contexts.clear();
    v_block_channels.clear();
    v_block_rep_types.clear();
    v_probes.clear();
    //
    spfFile.set(0,spfFile.m_name_cidx,pl.get_name_string());
    spfFile.set(0,spfFile.m_type_cidx,pl.get_type());
    spfFile.set(0,spfFile.m_v2_num_blocks_cidx,pl.block_cnt());
    spfFile.set(0,spfFile.m_num_match_cidx,pl.get_numMatch());

    // construct the comma seperated block info.
    for (int i=0;i<pl.block_cnt();i++) {
      v_block_sizes.push_back(pl.get_blockSize(i));
      v_block_annotations.push_back(pl.get_blockAnn(i));
      //
      v_block_alleles.push_back(pl.get_blockAllele(i));
      v_block_contexts.push_back(pl.get_blockContext(i));
      v_block_channels.push_back(pl.get_blockChannel(i));
      v_block_rep_types.push_back(pl.get_blockRepType(i));
    }
    spfFile.set(0,spfFile.m_v2_block_sizes_cidx,v_block_sizes);
    spfFile.set(0,spfFile.m_v2_block_annotations_cidx,v_block_annotations);
    //
    if (spfFile.m_v2_block_alleles_cidx>=0) {
      spfFile.set(0,spfFile.m_v2_block_alleles_cidx,v_block_alleles);
    }
    if (spfFile.m_v2_block_contexts_cidx>=0) {
      spfFile.set(0,spfFile.m_v2_block_contexts_cidx,v_block_contexts);
    }
    if (spfFile.m_v2_block_channels_cidx>=0) {
      spfFile.set(0,spfFile.m_v2_block_channels_cidx,v_block_channels);
    }
    if (spfFile.m_v2_block_rep_types_cidx>=0) {
      spfFile.set(0,spfFile.m_v2_block_rep_types_cidx,v_block_rep_types);
    }
    //
    for (int i=0;i<pl.probe_cnt();i++) {
      v_probes.push_back(FROM_ZEROBASED(pl.get_probeId(i)));
    }
    spfFile.set(0,spfFile.m_v2_num_probes_cidx,(uint64_t)v_probes.size());
    spfFile.set(0,spfFile.m_v2_probes_cidx,v_probes);
    //
    spfFile.writeLevel(0);
  }
}

void ProbeListFactory::writeProbeListToSpfFile_v3(affx::SpfFile& spfFile, ProbeListPacked pl)
{
  int probes_per_block;
  int p_i;
  int block_allele_code;
  int last_allele_code=-1;

  //
  spfFile.set(0,spfFile.m_name_cidx,pl.get_name_string());
  spfFile.set(0,spfFile.m_type_cidx,pl.get_type());
  spfFile.set(0,spfFile.m_num_match_cidx,pl.get_numMatch());
  spfFile.writeLevel(0);
  //
  probes_per_block=((pl.block_cnt()>0)?(pl.probe_cnt()/pl.block_cnt()):0);
  p_i=0;
  //
  for (int b_i = 0; b_i < pl.block_cnt(); b_i++) {
      block_allele_code = pl.get_blockAllele(b_i);
      if ( (b_i==0) || (last_allele_code != block_allele_code) ) {
          spfFile.set(1, spfFile.m_v3_allele_lbl_cidx, "allele");
          spfFile.set(1, spfFile.m_v3_allele_code_cidx, block_allele_code);
          spfFile.writeLevel(1);
          last_allele_code=block_allele_code;
      }
      //
      spfFile.set(2, spfFile.m_v3_context_lbl_cidx, "context");
      spfFile.set(2, spfFile.m_v3_annotation_code_cidx, pl.get_blockAnn(b_i));
      spfFile.set(2, spfFile.m_v3_context_code_cidx, pl.get_blockContext(b_i));
      if (spfFile.m_v3_channel_code_cidx>=0) {
          spfFile.set(2, spfFile.m_v3_channel_code_cidx, pl.get_blockContext(b_i));
      }
      spfFile.writeLevel(2);
      //
      for (int i=0;i<pl.get_blockSize(b_i);i++) {
          // FROM_ZEROBASED evals its arg twice.
          int probeid=pl.get_probeId(p_i++);
          spfFile.set(3, spfFile.m_v3_probe_id_cidx, FROM_ZEROBASED(probeid));
          spfFile.writeLevel(3);
      }
  }
}

void ProbeListFactory::writeToSpfFile_v3(affx::SpfFile& spfFile) const
{
  ProbeListPacked pl;

  for (int pl_i=0;pl_i<m_probelist_vec.size();pl_i++) {
      pl=m_probelist_vec[pl_i];
      writeProbeListToSpfFile_v3(spfFile, pl);
  }

}

void ProbeListFactory::writeProbeListToSpfFile_v4(affx::SpfFile &spfFile, ProbeListPacked pl)
{
  int probes_per_block;
  int p_i;
  spfFile.set(0,spfFile.m_name_cidx,pl.get_name_string());
  spfFile.set(0,spfFile.m_type_cidx,pl.get_type());
  spfFile.set(0,spfFile.m_num_match_cidx,pl.get_numMatch());
  spfFile.writeLevel(0);
  //
  probes_per_block=((pl.block_cnt()>0)?(pl.probe_cnt()/pl.block_cnt()):0);
  p_i=0;
  for (int b_i=0;b_i<pl.block_cnt();b_i++) {
    spfFile.set(1,spfFile.m_v4_block_lbl_cidx,"block");
    spfFile.set(1,spfFile.m_v4_annotation_code_cidx,pl.get_blockAnn(b_i));
    spfFile.set(1,spfFile.m_v4_allele_code_cidx,pl.get_blockAllele(b_i));
    spfFile.set(1,spfFile.m_v4_context_code_cidx,pl.get_blockContext(b_i));
    spfFile.set(1,spfFile.m_v4_channel_code_cidx,pl.get_blockChannel(b_i));
    spfFile.set(1,spfFile.m_v4_rep_type_cidx,pl.get_blockRepType(b_i));
    spfFile.writeLevel(1);
    //
    for (int i=0;i<probes_per_block;i++) {
      // FROM_ZEROBASED evals its arg twice.
      int probeid=pl.get_probeId(p_i++);
      spfFile.set(2,spfFile.m_v4_probe_id_cidx,FROM_ZEROBASED(probeid));
      spfFile.writeLevel(2);
    }
  }
}

void ProbeListFactory::writeToSpfFile_v4(affx::SpfFile& spfFile) const
{
  ProbeListPacked pl;


  for (int pl_i=0;pl_i<m_probelist_vec.size();pl_i++) {
    pl=m_probelist_vec[pl_i];
    //
    writeProbeListToSpfFile_v4(spfFile, pl);
  }
}


void ProbeListFactory::writeToSpfFile(affx::SpfFile& spfFile) const
{
  //
  if (spfFile.m_spf_format==2) {
    writeToSpfFile_v2(spfFile);
  }
  else if (spfFile.m_spf_format==3) {
    writeToSpfFile_v3(spfFile);
  }
  else if (spfFile.m_spf_format==4) {
    writeToSpfFile_v4(spfFile);
  }
  else {
    APT_ERR_ABORT("bad format");
  }
}

//////////

void ProbeListFactory::writeSpfFile(const std::string& spfName,
                                    int spfFormat,
                                    int with_allelecontext,
                                    int with_channel) const
{
  affx::SpfFile spfFile;
  // set some info which will be in the headers.
  spfFile.m_header_chipTypes=m_chipTypes;
  spfFile.m_header_numCols=m_numCols;
  spfFile.m_header_numRows=m_numRows;
  spfFile.m_header_numChannels=m_numChannels;
  spfFile.m_header_numProbesets=getProbeSetCount();
  // these will be used to define the format.
  spfFile.m_has_allele_info=with_allelecontext;
  spfFile.m_has_context_info=with_allelecontext;
  spfFile.m_has_channel_info=with_channel;
  // start writing.
  spfFile.writeSpf(spfName,spfFormat);
  //
  writeToSpfFile(spfFile);
  //
  spfFile.close();
}

void ProbeListFactory::writeSpfFile(const std::string& spfName,int spfFormat) const
{
  writeSpfFile(spfName,spfFormat,1,1);
}

//////////

/// This this the comparison function for sorting the names.
/// It dereferences the 'plf' to get the names to pass to strcmp.
class ProbeListFactory_NameIdxCmp {
public:
  const ProbeListFactory* m_plf;
  ProbeListFactory_NameIdxCmp(const ProbeListFactory* plf) {
    m_plf=plf;
  }
  bool operator()(int a_idx,int b_idx) {
    const char* a_name=m_plf->getProbeListNameCstr(a_idx);
    const char* b_name=m_plf->getProbeListNameCstr(b_idx);
    int rv=strcmp(a_name,b_name);
    if (rv<0) {
      return true;
    }
    return false;
  };
};

void ProbeListFactory::sortName2Idx() const
{
  // Because the index is built on demand, the message below
  // tends to show up in unexpected places -- like in the middle
  // of a progress meter on processing probesets
  // You can call "ensureSortedName2Idx()" yourself if you want to control
  // when this message is sent.
  Verbose::out(3,"Building and sorting probeset index.");
  m_name2idx_sortcnt++;

  // rebuild the index
  m_name2idx_vec.clear();
  m_name2idx_vec.reserve(m_probelist_vec.size());
  for (int i=0;i<m_probelist_vec.size();i++) {
    m_name2idx_vec.push_back(i);
  }
  //
  APT_ERR_ASSERT(m_name2idx_vec.size()==m_probelist_vec.size(),
                 "internal error. "
                 "(sort_cnt="+ToStr(m_name2idx_sortcnt)+
                 " m_probelist_vec="+ToStr(m_probelist_vec.size())+
                 " m_name2idx_vec="+ToStr(m_name2idx_vec.size()));
  //
  ProbeListFactory_NameIdxCmp comparer(this);
  sort(m_name2idx_vec.begin(),m_name2idx_vec.end(),comparer);
  m_name2idx_issorted=true;
}

void ProbeListFactory::assertName2IdxSorted(int verbose) const
{
  APT_ERR_ASSERT(m_name2idx_vec.size()==m_probelist_vec.size(),"internal error!");

  int m_name2idx_vec_size=m_name2idx_vec.size();

  if (m_name2idx_vec_size>0) {
    const char* this_name=NULL;
    const char* last_name=m_probelist_vec[m_name2idx_vec[0]].get_name_cstr();
    //
    if (verbose==1) {
      printf("%d : '%s'\n",0,last_name);
    }
    for (int i=1;i<m_name2idx_vec_size;i++) {
      this_name=m_probelist_vec[m_name2idx_vec[i]].get_name_cstr();
      if (verbose==1) {
        printf("%d : '%s'\n",i,this_name);
      }
      APT_ERR_ASSERT(strcmp(last_name,this_name)<=0,"Unsorted vector.");
      last_name=this_name;
    }
  }
}

void ProbeListFactory::ensureSortedName2Idx() const
{
  if (m_name2idx_issorted==false) {
    sortName2Idx();
    //
    assertName2IdxSorted(0); // debugging.
  }
  
  APT_ERR_ASSERT(m_name2idx_vec.size()==m_probelist_vec.size(),
                 "internal error. "
                 " sort_cnt="+ToStr(m_name2idx_sortcnt)+
                 " m_probelist_vec="+ToStr(m_probelist_vec.size())+
                 " m_name2idx_vec="+ToStr(m_name2idx_vec.size()));
}

/// This changes the order of the probes in m_probelist_vec to be
/// in name order.
void ProbeListFactory::sortFactoryByName()
{
  ensureSortedName2Idx();

  int m_probelist_vec_size=m_probelist_vec.size();

  std::vector<ProbeListPacked> tmp_vec;
  tmp_vec.reserve(m_probelist_vec_size);

  for (int i=0;i<m_probelist_vec_size;i++) {
    int idx=m_name2idx_vec[i];
    APT_ERR_ASSERT(((idx>=0)&&(idx<m_probelist_vec_size)),"internal error!");
    tmp_vec.push_back(m_probelist_vec[idx]);
  }
  m_probelist_vec=tmp_vec;
}

/////

/*
class ProbeListFactory_NameFinder {
public:
  const ProbeListFactory* m_plf;
  std::string m_name;
  const char* m_name_cstr;

  ProbeListFactory_NameFinder(const ProbeListFactory* plf,const std::string& name) {
    m_plf=plf;
    setSearchName(name);
  };
  void setSearchName(const std::string& name) {
    m_name=name;
    m_name_cstr=m_name.c_str();
  }
  // return true if a_idx is < then m_name_cstr
  bool operator()(int a_idx,int pretend_m_name_cstring_was_here) {
    const char* a_name=m_plf->getProbeListNameCstr(a_idx);
    int rv=strcmp(a_name,m_name_cstr);
    //printf("'%s'<=>'%s' => %d\n",a_name,m_name_cstr,rv);
    if (rv<0) {
      return true;
    }
    return false;
  };
};
*/
////

int ProbeListFactory::getProbeListIndexByName(const std::string& nameToFind) const
{
  int name_idx;
  int probelist_idx;
  int idx_max;
  std::vector<int>::iterator i;

  ensureSortedName2Idx();

  APT_ERR_ASSERT(m_name2idx_vec.size()==m_probelist_vec.size(),"internal error.");

  // The MS library breaks this by reversing the comparison order
  // ProbeListFactory_NameFinder finder(this,nameToFind);
  // i=lower_bound(m_name2idx_vec.begin(),m_name2idx_vec.end(),0,finder);
  // name_idx=i-m_name2idx_vec.begin();

  // here is our own binary seach
  {
    const char* nameToFind_cstr=nameToFind.c_str();
    //
    int i_begin=0;
    int i=i_begin;
    int i_end=m_name2idx_vec.size();

    while (i_begin<i_end) {
      i=i_begin+(i_end-i_begin)/2;
      const char* i_name=getProbeListNameCstr(m_name2idx_vec[i]);
      int rv=strcmp(i_name,nameToFind_cstr);
      //
      //printf("%3d--%3d--%3d : '%s'<=>'%s' => %2d\n",i_begin,i,i_end,i_name,nameToFind_cstr,rv);
      //
      if (rv==0) { // equal, but maybe not the first.
        break;
      }
      if (rv<0) { // i_name<nameToFind
         i_begin=i+1;
      }
      else { // bigger
        i_end=i;
      }
    }
    // back up to the first one.
    while ((0<i)&&(strcmp(getProbeListNameCstr(m_name2idx_vec[i-1]),nameToFind_cstr)==0)) {
      i--;
    }
    name_idx=i;
    //printf("name_idx=%d\n",name_idx);
  }

  // turn the name_idx into a probelist_idx
  idx_max=m_probelist_vec.size();
  //Verbose::out(1,"name_idx = " + ToStr(name_idx));
  //Verbose::out(1,"idx_max = " + ToStr(idx_max));
  if ((name_idx>=0) && (name_idx<idx_max)) {
    probelist_idx=m_name2idx_vec[name_idx];
    //Verbose::out(1,"probelist_idx = " + ToStr(probelist_idx));
    //Verbose::out(1,"nameToFind = " + ToStr(nameToFind));
    //Verbose::out(1,"found nane = " + ToStr(getProbeListName(probelist_idx)));
    if (nameToFind==getProbeListName(probelist_idx)) {
      return probelist_idx;
    }
  }
  return -1;
}

///
ProbeListPacked ProbeListFactory::getProbeListByName(const std::string& name) const
{
  return getProbeListAtIndex(getProbeListIndexByName(name));
}
// We dont want to do abort then it is not found as this function is used to
// test for existence as well as getting the value.
//    int idx = getProbeListIndexByName(name);
//    if (idx < 0) {
//        Err::errAbort("Couldn't find probeset with name: " + name);
//    }
//    return getProbeListAtIndex(idx);

/////

ProbeListPacked ProbeListFactory::getProbeListAtIndex(int i) const
{
  if ((i>=0)&&(i<m_probelist_vec.size())) {
    return m_probelist_vec[i];
  }
  //
  //Verbose::out(1,"Did not find probelist at index " + ToStr(i));
  return ProbeListPacked(NULL);
}

/////

std::string ProbeListFactory::getProbeListName(int i) const
{
  APT_ERR_ASSERT(i>=0,"too small!");
  APT_ERR_ASSERT(i<m_probelist_vec.size(),"too big!");
  return m_probelist_vec[i].get_name_string();
}
const char* ProbeListFactory::getProbeListNameCstr(int i) const
{
  APT_ERR_ASSERT(i>=0,"too small!");
  APT_ERR_ASSERT(i<m_probelist_vec.size(),"too big!");
  return m_probelist_vec[i].get_name_cstr();
}

/////

void ProbeListFactory::incApid(int cnt)
{
  m_probecnt_total+=cnt;
}

unsigned int ProbeListFactory::getApidMax() const
{
  return m_probecnt_total;
}


int
ProbeListFactory::getProbeSetCount() const
{
  return m_probelist_vec.size();
}

int
ProbeListFactory::findProbeApidByName(const std::string& qName,int qPid,int qAllele,int qContext,int qChannel) const
{
  int qPlIdx=getProbeListIndexByName(qName);
  if (qPlIdx<0) {
    return -1;
  }
  return findProbeApidByIdx(qPlIdx,qPid,qAllele,qContext,qChannel);
}

int
ProbeListFactory::findProbeApidByIdx(int qPlIdx,int qPid,int qAllele,int qContext,int qChannel) const
{
  if ((qPlIdx<0)||(qPlIdx>=m_probelist_vec.size())) {
    return -1;
  }
  return m_probelist_vec[qPlIdx].findProbeApid(qPid,qAllele,qContext,qChannel);
}

//////////

void ProbeListFactory::setDimensions(int cols,int rows)
{
  m_numCols=cols;
  m_numRows=rows;
}
void ProbeListFactory::getDimensions(int cols,int rows) const
{
  cols=m_numCols;
  rows=m_numRows;
}

void ProbeListFactory::setNumChannels(int num_channels)
{
  m_numChannels=num_channels;
}

int ProbeListFactory::numChannels() const
{
  return m_numChannels;
}


/// @brief     Build an index mapping ProbeId to Apid of _unique_ probeids
///            Having a non-unique mapping is ok, using it is not.
void
ProbeListFactory::buildProbeId2ApidIndex() const
{
  // allocate our index
  if (m_probeid2apid_vecptr==NULL) {
    m_probeid2apid_vecptr=new std::vector<int>();
  }

  // Meanings for (probeid)=>apid.
  //   (In the order the values move through.)
  // -1  => no apid for this ProbeId
  //  0+ => the apid for this probe.
  //        It is unique, as we have only seen it once.
  // -2  => this probeid appears in multiple probesets.
  //        We have seen it more than once.
  //        the apids are unique, but the mapping of (probeid)=>apid is not.
  //        the mapping of (probeset,probeid)=>apid is.

  // start the index with all '-1's.
  int pid_max=m_numCols*m_numRows;
  m_probeid2apid_vecptr->resize(pid_max,-1); //

  // walk the packed probesets to build our index
  int psi_max=getProbeSetCount();
  for (int psidx=0;psidx<psi_max;psidx++) {
    ProbeListPacked psp=getProbeListAtIndex(psidx);
    // now walk the probes of this set.
    // the apids are assigned from the start of the block.
    int pidx_max=psp.probe_cnt();
    int apid=psp.getApidStart();
    for (int pidx=0;pidx<pidx_max;pidx++) {
      int probeid=psp.get_probeId(pidx);
      // ignore null probes.
      if (probeid==ProbeList::null_probe) {
        continue;
      }
      // we dont allow duplicates (Has the index been set already?)
      int oldval=(*m_probeid2apid_vecptr)[probeid];
      if (oldval==-1) { // unseen, enter the mapping.
        (*m_probeid2apid_vecptr)[probeid]=apid;
      }
      else if (oldval>=0) { // seen once already, this isnt unique. Flag and warn.
        Verbose::out(1,"buildProbeId2ApidIndex: probeId "+ToStr(probeid)+" has already been used in another probeset. (apid was:"+ToStr(oldval));
        (*m_probeid2apid_vecptr)[probeid]=-2;
      }
      else if (oldval==-2) { // already warned.
        // do nothing
      }
      else {
        APT_ERR_ABORT("internal error.");
      }
      // the apids are sequental, advance it.
      apid++;
    } // for probe
  } // for probeset
}

int
ProbeListFactory::findProbeApidByProbeId(int qPid) const
{
  // in bounds?
  int pid_max=m_numCols*m_numRows;
  if ((qPid<0)||(qPid>=pid_max)) {
    APT_ERR_ABORT("ProbeId "+ToStr(qPid)+" out of bounds. (0-"+ToStr(pid_max)+")");
  }

  // has the index been built?
  if (m_probeid2apid_vecptr==NULL) {
    buildProbeId2ApidIndex();
  }

  // simply look it up.
  int apid=(*m_probeid2apid_vecptr)[qPid];
  // "-1" is an ok return value. "-2" is not as that means it is ambiguous.
  // wrap in if statement otherwise we have performance overhead of all the ToStr()
  if(apid <= -2) {
    APT_ERR_ASSERT(apid!=-2,"ProbeId "+ToStr(qPid)+" is not unique. (apid="+ToStr(apid)+")");
    APT_ERR_ASSERT(apid>-2,"ProbeId "+ToStr(qPid)+" has unexpected apid value. (apid="+ToStr(apid)+")");
  }
  return apid;
}
