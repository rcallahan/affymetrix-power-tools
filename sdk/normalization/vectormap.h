////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#ifndef __VECTORMAP_H_
#define __VECTORMAP_H_

//
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <utility>
#include <vector>

/// @file   vectormap.h
/// @brief  An STL container like "map" but with reduced space and
///         a requirement to sort before lookups.

//! The standard STL "map" keeps its elements in a sorted
//! tree.  The norm code only needs to sort the data before
//! use and a vector uses less memory.  The disadvantage to
//! vectormap is the elements are not kept in order.  Before
//! any lookups are done the vector is sorted.  If the
//! inserts and lookups are alternated, much time will be
//! spent sorting.

using namespace std;

////////////////////

template <typename T1> class pair_key_asc {
public:
  /// @brief        Compare pairs by the key
  /// @return       true if pair1.key < pair2.key
  /// @remark       The key is the first value of the pair.
  bool operator()(const T1& pair1, const T1& pair2) const {
    return pair1.first < pair2.first;
  }
};

template <typename T1,typename T2> class pair_cmp_key {
public:
  /// @brief       compare a pair to a key value. 
  /// @return      Returns (-1,0,+1)
  int operator()(const T1& pair, const T2& key) const {
    if (pair.first < key) {
      return +1;
    }
    if (key < pair.first) {
      return -1;
    }
    if (pair.first==key) {
      return 0;
    }
    assert(0);
  }
};

template <typename T1> class pair_cmp_key_bool {
 public:
  /// @brief       A weak ordering.
  /// @return      true if pair1.key < pair2.key
  bool operator()(const T1 & pair, const T1 & key) const {
    if (pair.first < key.first) {
      return true;
    }
    return false;
  }
};
      

//! vectormap == A map with the KeyVal pairs stored in a vector.
//! the vector is sorted when needed.

template <typename Key,typename Val>
class vectormap
{
 public:
  //! the key type of the map
  typedef Key key_t;
  //! the value type of the map
  typedef Val val_t;
  typedef pair<key_t,val_t> keyval_t;
  typedef typename vector<keyval_t>::iterator iterator;
  typedef typename vector<keyval_t>::iterator iterator_t;

  //! All the pairs are stored here.
  vector<keyval_t> vec;
  int is_sorted;

  vectormap() {
    clear();
  };

  /// @brief start iterator of the vector
  iterator_t begin() { return vec.begin(); }
  /// @brief end iterator of the vector
  iterator_t end()   { return vec.end();   }

  /// @brief       Random access to the map contents.
  /// @param       n         the index
  /// @return      the keyval at that index
  /// @remark      Contents might not be sorted
  keyval_t& operator[](size_t n) {
    return vec[n];
  }

  /// @brief       Removes all the data
  void clear() { vec.clear(); is_sorted=1; };
  /// @return      the size of the map
  size_t size() { return vec.size(); }

  /// @brief       adds the (key,val) to the map
  /// @param       key       
  /// @param       val       
  /// @remark      The pair is created for you.
  void insert_kv(key_t key,val_t val) {
    is_sorted=0;
    vec.push_back(pair<key_t,val_t>(key,val));
  }

  /// @brief       adds a pair(key,val) to the map
  /// @param       keyval    a pair of (key,val)
  void insert(keyval_t keyval) {
    is_sorted=0;
    vec.push_back(keyval);
  };
  
  /// @brief       sort the map before use. Sets is_sorted to true.
  void sort_key_asc() {
    //printf("SORT!====(before)\n"); print();
    pair_key_asc<keyval_t> func;
    sort(vec.begin(),vec.end(),func);
    is_sorted=1;
    //printf("SORT!====(after)\n"); print();
  };

  // make sure the vectormap is sorted.
  // However! normalization does its own sorting.
  // so dont sort the vectors 
  void ensure_sorted() {
    if (is_sorted==1) {
       return;
    }
    sort_key_asc();
  };

  iterator_t vlower_bound(const key_t key) {
    pair_cmp_key_bool<keyval_t> cmpfunc;
    //ensure_sorted();
    return lower_bound(begin(),end(),keyval_t(key,val_t()),cmpfunc);
  }
  iterator_t vupper_bound(const key_t key) {
    pair_cmp_key_bool<keyval_t> cmpfunc;
    //ensure_sorted();
    return upper_bound(begin(),end(),keyval_t(key,val_t()),cmpfunc);
  }
  pair <iterator_t,iterator_t> equal_range(const key_t key) {
    //ensure_sorted();
    return make_pair<iterator_t,iterator_t>(vlower_bound(key),vupper_bound(key));
  }
  
  /// @brief       Dumps the contents of the map for debugging.
  void print() {
    printf("--------------------\n");
    for (size_t i=0;i<vec.size();i++) {
      printf(" %3d : %3d => %3d\n",int(i),int(vec[i].first),int(vec[i].second));
    }
  }
};

#endif // __VECTORMAP_H_
