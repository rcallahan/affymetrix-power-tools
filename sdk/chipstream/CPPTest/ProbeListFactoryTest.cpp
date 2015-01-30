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

/**
 * @file   ProbeListFactoryTest.cpp
 * @author harley
 */

//
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeListStl.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstdio>
#include <cstring>
#include <string>
#include <stdio.h>
#include <vector>
//

/**
 * @class ProbeListFactoryTest
 * @brief
 */

class ProbeListFactoryTest : public CppUnit::TestFixture {
public:
  ProbeListFactoryTest() {}

  CPPUNIT_TEST_SUITE( ProbeListFactoryTest );
  CPPUNIT_TEST( plf_test_1 );
  CPPUNIT_TEST( plf_test_2 );
  CPPUNIT_TEST( plf_test_3 );
  CPPUNIT_TEST( plf_test_4 );
  CPPUNIT_TEST_SUITE_END();

  ///
  void plf_test_1();
  void plf_test_2();
  void plf_test_3();
  void plf_test_4();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( ProbeListFactoryTest );

void
ProbeListFactoryTest::plf_test_1()
{
  ProbeListFactory* plf=new ProbeListFactory();
  CPPUNIT_ASSERT(plf!=NULL);
  delete plf;
}

void
ProbeListFactoryTest::plf_test_2()
{
  ProbeListFactory* plf=new ProbeListFactory();
  CPPUNIT_ASSERT(plf!=NULL);

  //
  ProbeListStl pls;
  char buf[100];
  const int l_num=50;
  int b_num=2;
  int p_num=5;
  int p_total=0;
  
  for (int l=0;l<l_num;l++) {
    pls.clear();
    // allocate the names in order so they make sense
    sprintf(buf,"pls-%02d",l);
    pls.set_name(buf);
    //
    for (int b=0;b<b_num;b++) {
      pls.push_block(0,0,0,b);
      pls.set_blockSize(b,p_num);
      for (int p=0;p<p_num;p++) {
        pls.push_probe(p_total++);
      }
    }
    //pls.dump();
    plf->add_ProbeList(pls);
  }
  // plf->dump_probelists();

  // handy helpers
#define PRINT_LBL(_lbl) { printf("\n=== '%s'\n",_lbl); }
#define PRINT_FIND(_expr,_expected) { \
    int rv=_expr; \
    printf("'%s'==%d\n",#_expr,rv); \
    fflush(NULL); \
  }
  //    CPPUNIT_ASSERT(rv==_expected); 

  //
  PRINT_LBL("postitive");
  PRINT_FIND(plf->findProbeApidByName("pls-00", 0,0,0,0),0);
  PRINT_FIND(plf->findProbeApidByName("pls-00", 4,0,0,0),4);
  PRINT_FIND(plf->findProbeApidByName("pls-00", 5,0,0,1),5);
  PRINT_FIND(plf->findProbeApidByName("pls-00", 9,0,0,1),9);
  //
  PRINT_FIND(plf->findProbeApidByName("pls-01",10,0,0,0),10);
  PRINT_FIND(plf->findProbeApidByName("pls-01",14,0,0,0),14);
  PRINT_FIND(plf->findProbeApidByName("pls-01",15,0,0,1),15);
  PRINT_FIND(plf->findProbeApidByName("pls-01",19,0,0,1),19);
  //
  PRINT_FIND(plf->findProbeApidByName("pls-48",480,0,0,0),480);
  PRINT_FIND(plf->findProbeApidByName("pls-49",490,0,0,0),490);
  //
  PRINT_LBL("negative");
  PRINT_FIND(plf->findProbeApidByName("", 0,0,0,1),-1);
  PRINT_FIND(plf->findProbeApidByName("NoName", 0,0,0,1),-1);
  PRINT_FIND(plf->findProbeApidByName("pls-00", 0,0,0,1),-1);
  PRINT_FIND(plf->findProbeApidByName("pls-00", 0,0,0,1),-1);
  PRINT_FIND(plf->findProbeApidByName("pls-01", 0,0,0,0),-1);
  PRINT_FIND(plf->findProbeApidByName("pls-01", 0,0,0,0),-1);
  PRINT_FIND(plf->findProbeApidByName("pls-99", 0,0,0,1),-1);

  delete plf;
}


#define TEST3_PUSH_PROBES(_pls,_cnt) { for (int i=0;i<_cnt;i++) { _pls.push_probe(pid++); } }

void
ProbeListFactoryTest::plf_test_3()
{
  ProbeListFactory* plf=new ProbeListFactory();
  CPPUNIT_ASSERT(plf!=NULL);

  ProbeListStl pls;
  ProbeListPacked plp;
  int pid=0;

  ////
  plf->clear();
  pls.clear();
  pls.set_name("test3-1");
  pls.set_numMatch(1); // 10 PM
  pls.push_block(0,0,0,0);
  pls.set_blockSize(0,10);
  TEST3_PUSH_PROBES(pls,10);

  

  // add it and get PLP back.
  plf->add_ProbeList(pls);
  plp=plf->getProbeListAtIndex(0);
  //plp.dump();



  //
  for (int i=0;i<plp.probe_cnt();i++) {
    CPPUNIT_ASSERT(plp.isPm(i)==true);
  }

  //// test again
  plf->clear();
  pls.clear();
  pls.set_name("test3-2");
  pls.set_numMatch(2); // 5 PM and 5 MM
  pls.push_block(0,0,0,0);
  pls.set_blockSize(0,10);
  TEST3_PUSH_PROBES(pls,10);

  // add it and get PLP back.
  plf->add_ProbeList(pls);
  plp=plf->getProbeListAtIndex(0);
  //plp.dump();

  //
  for (int i=0;i<plp.probe_cnt();i++) {
    if (i<5) {
      CPPUNIT_ASSERT(plp.isPm(i)==true);
    }
    else {
      CPPUNIT_ASSERT(plp.isPm(i)==false);
    }
  }

  ////
  plf->clear();
  pls.clear();
  pls.set_name("test3-3");
  pls.set_numMatch(2); // 5 PM and 5 MM
  pls.push_block(0,0,0,0);
  pls.set_blockSize(0,10);
  pls.push_block(0,1,0,0);
  pls.set_blockSize(1,10);
  TEST3_PUSH_PROBES(pls,20);

  // add it and get PLP back.
  plf->add_ProbeList(pls);
  plp=plf->getProbeListAtIndex(0);
  //plp.dump();

  //
  for (int i=0;i<plp.probe_cnt();i++) {
    if ((i%10)<5) {
      CPPUNIT_ASSERT(plp.isPm(i)==true);
    }
    else {
      CPPUNIT_ASSERT(plp.isPm(i)==false);
    }
  }

  delete plf;
}

void
ProbeListFactoryTest::plf_test_4()
{
  ProbeListFactory* plf=new ProbeListFactory();
  CPPUNIT_ASSERT(plf!=NULL);

  ProbeListStl pls;
  ProbeListPacked plp;
  std::vector<int> probeIds;
  int pid=0;

  ///
  for (int bsize=0;bsize<10;bsize++) {
    plf->clear();
    pls.clear();
    pls.set_name("test4-0");
    pls.set_numMatch(1); // 10 PM
    pls.push_block(0,0,0,0);
    pls.set_blockSize(0,bsize);
    pid=0;
    TEST3_PUSH_PROBES(pls,bsize);
    //
    plp=plf->add_ProbeList(pls);
    plp.get_probeIds(probeIds);
    CPPUNIT_ASSERT(probeIds.size()==bsize);
    if (bsize>=1) {
      CPPUNIT_ASSERT(probeIds[0]==0);
      CPPUNIT_ASSERT(probeIds[bsize-1]==bsize-1);
    }
  }
  //
  delete plf;
}
