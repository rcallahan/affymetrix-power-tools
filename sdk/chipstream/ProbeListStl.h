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

#ifndef _CHIPSTREAM_PROBELISTSTL_H_
#define _CHIPSTREAM_PROBELISTSTL_H_

// forward declare as the included headers use us.
class ProbeListStl;
class ProbeListStlBlock;
class ProbeListStlProbe;

#include "chipstream/ProbeListFactory.h"
//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/// @file   ProbeListStl.h
/// @brief  A resizeable and modifiable ProbeList object.
///         This object should be keept in sync with ProbeListPacked.


/// This is used as an interchange format

/// @brief     A probelist made from stl containers.
///            This is like ProbeListPacked, but it can be resized and modified.
///            See the SpfFile readers for example useage.

/// @todo Rename this to ProbeListMutable or something better.
/// @todo Put the probes in the blocks or make a bare block object for export.

///
class ProbeListStlProbe {
public:
  /// The 0-based probeid
  int    m_probeId;
  /// The GC count
  char   m_probeGc;
  /// The Apid of this probe.
  int    m_probeApid;

  /// A blank probe
  ProbeListStlProbe();
  /// @brief     construct a copy
  /// @param     p
  ProbeListStlProbe(const ProbeListStlProbe& p);
  /// @brief     Constructors
  /// @param     pid
  ProbeListStlProbe(int pid);
  /// @brief     Constructor
  /// @param     pid
  /// @param     gc
  /// @param     apid
  ProbeListStlProbe(int pid,int gc,int apid);

  /// zeros out the fields.
  void clear();

  /// @brief     probeId
  /// @return    probeId
  int id();
  /// @brief     returns "gc" count.
  /// @return    the gc count.
  int gc();
  /// @brief     analysis probe id
  /// @return    analysis probe id
  int apid();
  /// @brief
  /// @return    true of this is a null probe. (m_probeId==-1)
  bool isNull();
};

class ProbeListStlBlock {
public:
  int     m_blockSize;          ///< the number of probes in this block.
  short   m_blockAnnotation;    ///< the annotation_code for this block
  int     m_blockAllele;        ///< the allele_code for this block
  int     m_blockContext;       ///< the context_code for this block
  int     m_blockChannel;       ///< the channel_code for this block
  int     m_blockRepType;       ///< the replicate_type for this block

  // While we could keep the list of probes in each block,
  // we keep a single vector of probes in the ProbeListStl.
  // this reduces the vector count.
  // This is here as a note of another way of doing it.
  // std::vector <ProbeListStlProbe> m_probeid_vec;

  /// @brief     new default block
  ProbeListStlBlock();
  /// @brief     Make a copy of the block
  /// @param     b
  ProbeListStlBlock(const ProbeListStlBlock& b);
  /// @brief     create and push a block.
  /// @param     bAnn
  /// @param     bAllele
  /// @param     bContext
  /// @param     bChannel
  /// @param     bRepType
  ProbeListStlBlock(int bAnn,int bAllele, int bContext, int bChannel, int bRepType);

  /// @brief     does this block match the Allele, Context and Channel?
  /// @param     qAllele   
  /// @param     qContext  
  /// @param     qChannel  
  /// @return    true if they are all equal.
  bool matchACC(int qAllele,int qContext,int qChannel) const;

  /// @brief Clear out the block to unset values.
  void clear();
};

/// If
// : public ProbeList
/// @brief
class ProbeListStl

{
public:
  std::string m_name;  ///< the name of this probelist
  int m_max_name_len; ///< the longest name which could be stored.
  int m_type; ///< type
  int m_num_match; ///< number of matching probes
  std::vector<ProbeListStlBlock> m_block_vec; ///< collection of blocks in this ProbeList
  std::vector<ProbeListStlProbe> m_probe_vec; ///< collection of probes in this ProbeList

  //
  ProbeListStl();
  ProbeListStl(const ProbeListStl& pl);
  ProbeListStl(const ProbeListPacked& pl);
  ProbeListStl(const std::string& name);
  ProbeListStl(int block_cnt,int probe_cnt);
  ProbeListStl(int block_cnt,int probe_cnt,const std::string& name);

  /// @brief     Make ourselves a copy of a ProbeListPacked object.
  /// @param     plp
  /// @return
  ProbeListStl& operator=(const ProbeListPacked& plp);

  /// @brief     clear the contents of this object.
  void clear();
  /// @brief     set the size of this probelist
  /// @param     block_cnt
  /// @param     probe_cnt
  void resize(int block_cnt,int probe_cnt);
  /// @brief     reserve space in this probelist
  /// @param     block_cnt
  /// @param     probe_cnt
  void reserve(int block_cnt,int probe_cnt);

  /// these are the basic accessors, shared with ProbeListPacked
  /// see ProbeList.h for more info.

  /// These are just like ProbeListPacked.
  //
  int get_type() const;
  void set_type(int type);

  //
  int block_cnt() const;
  int probe_cnt() const;
  int max_name_len() const;
  //
  int get_numMatch() const;
  void set_numMatch(int numMatch);
  //
  int get_blockSize(int b) const;
  void set_blockSize(int b,int val);
  //
  short get_blockAnn(int b) const;
  void set_blockAnn(int b,short val);
  //
  int get_blockAllele(int b) const;
  void set_blockAllele(int b,int val);
  //
  int get_blockContext(int b) const;
  void set_blockContext(int b,int val);
  //
  int get_blockChannel(int b) const;
  void set_blockChannel(int b,char val);
  //
  int get_blockRepType(int b) const;
  void set_blockRepType(int b,char val);

  //
  int get_probeId(int b) const;
  void set_probeId(int b,int val);
  //
  int8_t get_probeGc(int d) const;
  void  set_probeGc(int d,int8_t val);

  int get_probeApid(int d) const;
  void  set_probeApid(int d,int val);

  //
  void set_name(const std::string &name);
  // this shouldnt be used.
  const char* get_name() const;
  // use this instead.
  const char* get_name_cstr() const;
  //
  std::string get_name_string() const;

  //// these accessors are in addition to the ProbeList accessors

  /// @brief
  /// @param     b
  void push_block(const ProbeListStlBlock& b);
  void push_block(const ProbeListStlBlock* b);
  /// @brief
  /// @param     bAnn
  /// @param     bAllele
  /// @param     bContext
  /// @param     bChannel
  /// @param     bRepType
  void push_block(int bAnn, int bAllele, int bContext, int bChannel, int bRepType = 0);

  // this is ugly.
  //void push_blockandprobesfrom(const ProbeListStl& pls,int bidx);

  /// @brief     push a probe onto the end of a probelist
  /// @param     p         probelist probe to push
  void push_probe(const ProbeListStlProbe& p);
  /// @brief     push a probe with just an id
  /// @param     id        id of the probe.
  void push_probe(int id);
  /// @brief     create and push a probe with these three values.
  /// @param     id
  /// @param     det
  /// @param     apid
  void push_probe(int id,int det,int apid);
  /// @brief     push an entire probelist onto this one.
  /// @param     pl
  void push_probelist(const ProbeListStl& pl);

  /// @brief     Get a pointer to the block B
  /// @param     b         
  /// @return    the pointer or NULL
  const ProbeListStlBlock* getBlockPtr(int b);

  /// @brief     Find the index of the block with matchin ACC
  /// @param     qAllele   
  /// @param     qContext  
  /// @param     qChannel  
  /// @return    the index or -1
  int findBlockMatch(int qAllele,int qContext,int qChannel) const;

  /// @brief     Return true if two or more blocks have matching ACCs
  /// @return    
  bool hasDuplicateBlocks() const;

  /// @brief     For debuggging...
  void dump() const;
};

#endif // _CHIPSTREAM_PROBELISTSTL_H_
