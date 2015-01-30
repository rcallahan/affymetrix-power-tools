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

#include "normalization/CPPTest/VectormapTest.h"
//
#include "normalization/vectormap.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
#include <cstdlib>
#include <ctime>
#include <map>
//



// ==============================
// util functions


// multimaps dont have "+/-" defined for iterators.
// We will just count them ourselves.
size_t mmap_iterator_idx(multimap<int,int>& mmap,multimap<int,int>::iterator & iter) {
  size_t idx=0;
  multimap<int,int>::iterator curr;

  curr=mmap.begin();
  while (curr!=iter) {
    idx++;
    curr++;
  }
  return idx;
}

//
size_t vmap_iterator_idx(vectormap<int,int>& vmap,vectormap<int,int>::iterator & iter) {
  return iter-vmap.begin();
}


//
void test_against_mmap_key(vectormap<int,int>& vmap,multimap<int,int>& mmap, int key) {
  typedef multimap<int,int>::iterator  mmap_iterator_t;
  typedef vectormap<int,int>::iterator vmap_iterator_t;
  mmap_iterator_t mmap_iter;
  vmap_iterator_t vmap_iter;
  size_t mmap_idx;
  size_t vmap_idx;

  //
  vmap.ensure_sorted();

  // look up the data
  mmap_iter=mmap.lower_bound(key);
  mmap_idx=mmap_iterator_idx(mmap,mmap_iter);
  vmap_iter=vmap.vlower_bound(key);
  vmap_idx=vmap_iterator_idx(vmap,vmap_iter);

  // verbose
  //if (mmap_iter==mmap.begin()) printf("mmap_i==begin\n");
  //if (mmap_iter==mmap.end())   printf("mmap_i==end\n");

  // check
  if (mmap_iter==mmap.begin()) 
    CPPUNIT_ASSERT(vmap_iter==vmap.begin());
  if (mmap_iter==mmap.end())   
    CPPUNIT_ASSERT(vmap_iter==vmap.end());

  //printf("LOWER (key=%3d) mmap idx=%3d  vmap idx=%3d\n",key,mmap_idx,vmap_idx);

  // look up the data
  mmap_iter=mmap.upper_bound(key);
  mmap_idx=mmap_iterator_idx(mmap,mmap_iter);
  vmap_iter=vmap.vupper_bound(key);
  vmap_idx=vmap_iterator_idx(vmap,vmap_iter);

  // verbose
  //if (mmap_iter==mmap.begin()) printf("mmap_i==begin\n");
  //if (mmap_iter==mmap.end())   printf("mmap_i==end\n");

  // check
  if (mmap_iter==mmap.begin())
    CPPUNIT_ASSERT(vmap_iter==vmap.begin());
  if (mmap_iter==mmap.end())
    CPPUNIT_ASSERT(vmap_iter==vmap.end());

  //printf("UPPER (key=%3d) mmap idx=%3d  vmap idx=%3d\n",key,mmap_idx,vmap_idx);

  // Compare equal_range
  //pair <mmap_iter,mmap_iter> mmap_range;
  //  pair <multimap<int,int>::iterator,multimap<int,int>::iterator> mmap_range;
  //  mmap_range=mmap.equal_range(key);
  /// @todo fix compare equal_range
  //  printf("equal_range(%3d)=mmap[%3d-%3d]\n",
  //       mmap_iterator_idx(mmap_range.first),mmap_iterator_idx(mmap_range.second));

}


// Test our map to the same data in the STL multimap
void test_against_mmap(vectormap<int,int>& vmap) {
  multimap<int,int> mmap;
  size_t vmap_size;

  // clone the data
  mmap.insert(vmap.begin(),vmap.end());

  //// Now compare the results of the method calls.

  // always these
  test_against_mmap_key(vmap,mmap,-100);
  test_against_mmap_key(vmap,mmap,  +0);
  test_against_mmap_key(vmap,mmap,+100);

  // Test the keys of the map...
  vmap_size=vmap.size();
  if (vmap_size!=0) {
    // ...from just beyond here...
    test_against_mmap_key(vmap,mmap,(vmap[0].first-1));
    // ...to just beyond there...
    test_against_mmap_key(vmap,mmap,(vmap[vmap_size-1].first+1));
    // ... and all points in between!
    for (unsigned int idx=0;idx<vmap_size;idx++) {
      test_against_mmap_key(vmap,mmap,vmap[idx].first);
    }
  }
}

// ==============================


CPPUNIT_TEST_SUITE_REGISTRATION( VectormapTest );

void VectormapTest::setUp()
{
}

void VectormapTest::tearDown()
{
}

void VectormapTest::test_Vectormap()
{
  vectormap<int,int> vmap;

  //printf("VectormapTest::test_Vectormap()...\n");

  // the empty case == {}
  vmap.clear();
  test_against_mmap(vmap);

  // sorted case == {1,2,3,4,5,6,7,8,9}
  vmap.clear();
  for (int i=0;i<10;i++) {
    vmap.insert_kv(i,i);
  }
  test_against_mmap(vmap);

  // duplicate case
  vmap.clear();
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      vmap.insert_kv(j,i);
    }
  }
  test_against_mmap(vmap);

  // big random case
  vmap.clear();
  srand( (unsigned)time( NULL ) );
  for (int i=0;i<300;i++) {
    vmap.insert_kv(rand()%1000,rand()%1000);
  }
  test_against_mmap(vmap);
  
}
