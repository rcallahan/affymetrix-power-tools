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

#include "chipstream/ProbeListFactory.h"
//
#include <cassert>
//

#define print_sizeof(x) printf("== sizeof(%s)=%lu\n",#x,sizeof(x))
#define STRORNULL(x) ((x==NULL)?"(null)":x)

int
main(int argc,char* argv[])
{
  //print_sizeof(ProbeList_Original);
  //
  print_sizeof(ProbeList_E_Probe);
  print_sizeof(ProbeList_E_Block);
  print_sizeof(ProbeList_Head);
  print_sizeof(ProbeList);

  // here we just alloc some bytes and print the pointers

  // size of the test
  int test_size=4;

  // malloc will round to 16b
  if (0) {
    for (int i=0;i<10;i++) {
      char* ptr=(char*)malloc(test_size);
      printf("malloc(%d)=%p\n",test_size,ptr);
    }
  }

  // ProbeList heap does not.
  // make the heap
  ProbeListFactory pl_factory;
  ProbeListPacked plist;

  // try changing size of the regions to see what happens.
  //pl_heap.set_default_region_size(1024*1024); // big
  //pl_heap.set_default_region_size(1024); // small
  //pl_heap.set_default_region_size(10);  // tiny

  pl_factory.dump_regions();

  // and get some byte pointers
  printf("== byte ptr test:\n");
  for (int i=0;i<10;i++) {
    char* ptr=pl_factory.alloc_ProbeList_bytes(test_size);
    printf("  heap_allooc(%d)=%p\n",test_size,ptr);
  }

  pl_factory.dump_regions();
  pl_factory.clear();
  pl_factory.dump_regions();

  for (int i=0;i<10000;i++) {
    pl_factory.alloc_ProbeList_bytes(100);
    pl_factory.alloc_ProbeList_bytes(50);
  }

  pl_factory.dump_regions();
  pl_factory.clear();
  pl_factory.dump_regions();

  // alloc some zero
  printf("== name length test:\n");
  for (int i=0;i<2;i++) {
    //
    plist=pl_factory.add_ProbeList(0,0,0);
    plist.set_name("abc");
    printf("  plist: addr=%p name='%s'\n",plist.m_headptr,STRORNULL(plist.get_name_cstr()));
    //
    plist=pl_factory.add_ProbeList(0,0,1);
    plist.set_name("def");
    printf("  plist: addr=%p name='%s'\n",plist.m_headptr,STRORNULL(plist.get_name_cstr()));
    //
    plist=pl_factory.add_ProbeList(0,0,2);
    plist.set_name("ghi");
    printf("  plist: addr=%p name='%s'\n",plist.m_headptr,STRORNULL(plist.get_name_cstr()));
  }

  printf("== heap clear test:\n");
  // use the heap 4 times, clearing at the end.
  // note the reuse of addresses after the first time.
  // (clearing of the allocations above this loop.)
  for (int x=0;x<10;x++) {
    // now generate some ProbeLists
    for (int i=0;i<100;i++) {
      //ProbeList plist=pl_factory.new_ProbeList(i,4,"foo123");
      plist=pl_factory.add_ProbeList(2,4,"foo123");
      if ((x<3)&&(i<5)) {
        printf("  plist: addr=%p name=%s\n",plist.m_headptr,STRORNULL(plist.get_name_cstr()));
      }
      // put some data in
      for (int b=0;b<plist.block_cnt();b++) {
        plist.set_blockSize(b,b);
      }
      for (int p=0;p<plist.probe_cnt();p++) {
        plist.set_probeId(p,p);
      }
    }
    pl_factory.clear();
  }
  pl_factory.delete_regions();

  return 0;
}
