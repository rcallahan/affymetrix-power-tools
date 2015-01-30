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

#ifndef _CHIPSTREAM_PROBELIST_H_
#define _CHIPSTREAM_PROBELIST_H_

// commented out as the methods are ifdefed out.
// #include <string>

// ProbeList isnt used, as objects with virtual functions
// have an extra ptr.  (Which doubles the size of ProbeListPacked.)

/// A place to park functions and constants common to
/// ProbeListPacked and ProbeListStl.
class ProbeList {
public:
  /// Annotations a specific block.
  enum BlockAnnotation {
    NoAnnotation,
    A_Allele,
    B_Allele,
    A_AlleleForward,
    B_AlleleForward,
    A_AlleleReverse,
    B_AlleleReverse
  };

  /// Use this symbolic value for an unset variable.
  static const int UNSET = -1;
  /// use this symbolic value for the null probe.
  static const int null_probe = -1;

  // ifdefed out as there are no users of these 
  // ProbeListPacked and ProbeListStl dont inherit from ProbeList.
  // (It would add a pointer for the virtual class for each object.
  // Which would make them too big.)

  // The intent is for both ProbeListPacked and ProbeListStl to have 
  // the same interface.  (As documented below.)

#if 0

  virtual int get_type() const = 0;
  virtual void set_type(int type) = 0;
  
  // accessors
  virtual int block_cnt() const = 0;
  virtual int probe_cnt() const = 0;
  virtual int max_name_len() const = 0;
  //
  virtual int get_numMatch() const = 0;
  virtual void set_numMatch(int numMatch) = 0;
  //
  /// @todo add set_XXX(std::vector<int>);
  virtual int get_blockSize(int b) const = 0;
  virtual void set_blockSize(int b,int val) = 0;
  //
  virtual short get_blockAnn(int b) const = 0;
  virtual void set_blockAnn(int b,short val) = 0;
  //
  virtual int get_blockAllele(int b) const = 0;
  virtual void set_blockAllele(int b,int val) = 0;
  //
  virtual int get_blockContext(int b) const = 0;
  virtual void set_blockContext(int b,int val) = 0;
  //
  virtual int get_blockChannel(int b) const = 0;
  virtual void set_blockChannel(int b,uint8_t val) = 0;

  //
  virtual int get_probe(int b) const = 0;
  virtual void set_probe(int b,int val) = 0;
  //
  virtual int8_t get_probeDet(int d) const = 0;
  virtual void  set_probeDet(int d,int8_t val) = 0;

  //
  virtual void set_name(const std::string &name) = 0;
  // this shouldnt be used.
  virtual const char* get_name() const {
    return get_name_cstr();
  };
  // use this instead.
  virtual const char* get_name_cstr() const = 0;
  //
  virtual std::string get_name_string() const = 0;

#endif
};

#endif // _CHIPSTREAM_PROBELIST_H_
