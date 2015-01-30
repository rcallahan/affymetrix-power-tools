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
 * @file   ProbeListStlTest.cpp
 * @author harley
 */

//
#include "chipstream/ProbeListStl.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @class ProbeListStlTest
 * @brief
 */

class ProbeListStlTest : public CppUnit::TestFixture {
public:
  ProbeListStlTest() {}

  CPPUNIT_TEST_SUITE( ProbeListStlTest );
  CPPUNIT_TEST( pls_test_1 );
  CPPUNIT_TEST( pls_test_2 );
  CPPUNIT_TEST( pls_test_3 );
  CPPUNIT_TEST( pls_test_4 );
  CPPUNIT_TEST_SUITE_END();

  ///
  void pls_test_1();
  void pls_test_2();
  void pls_test_3();
  void pls_test_4();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( ProbeListStlTest );

void
ProbeListStlTest::pls_test_1()
{
  ProbeListStl* pls=new ProbeListStl();
  CPPUNIT_ASSERT(pls!=NULL);
  delete pls;
}

void
ProbeListStlTest::pls_test_2()
{
  ProbeListStl* pls=new ProbeListStl();
  CPPUNIT_ASSERT(pls!=NULL);

  const ProbeListStlBlock* bptr=NULL;

  //
  const int b_cnt=5;
  for (int b=0;b<b_cnt;b++) {
    pls->push_block(b,b,b,b);
    CPPUNIT_ASSERT(pls->get_blockAnn(b)    ==b);
    CPPUNIT_ASSERT(pls->get_blockAllele(b) ==b);
    CPPUNIT_ASSERT(pls->get_blockContext(b)==b);
    CPPUNIT_ASSERT(pls->get_blockChannel(b)==b);
    //
    bptr=pls->getBlockPtr(b);
    CPPUNIT_ASSERT(bptr->matchACC(b,b,b)==true);
    CPPUNIT_ASSERT(bptr->matchACC(-1,-1,-1)==false);
  }
  CPPUNIT_ASSERT(pls->block_cnt()==b_cnt);

  pls->resize(0,0);
  CPPUNIT_ASSERT(pls->block_cnt()==0);

  
  delete pls;
}

void
ProbeListStlTest::pls_test_3()
{
  ProbeListStl* pls=new ProbeListStl();
  CPPUNIT_ASSERT(pls!=NULL);

  pls->push_block(0,0,0,0);
  pls->push_block(1,1,1,1);
  pls->push_block(2,2,2,2);
  CPPUNIT_ASSERT(pls->hasDuplicateBlocks()==false);

  pls->push_block(2,2,2,0);
  CPPUNIT_ASSERT(pls->hasDuplicateBlocks()==false);
  pls->push_block(2,2,2,1);
  CPPUNIT_ASSERT(pls->hasDuplicateBlocks()==false);

  pls->push_block(0,0,0,0);
  CPPUNIT_ASSERT(pls->hasDuplicateBlocks()==true);

  pls->resize(pls->block_cnt()-1,pls->probe_cnt());
  CPPUNIT_ASSERT(pls->hasDuplicateBlocks()==false);

  delete pls;
}

void
ProbeListStlTest::pls_test_4()
{
  ProbeListStl* pls=new ProbeListStl();
  CPPUNIT_ASSERT(pls!=NULL);

  CPPUNIT_ASSERT(pls->probe_cnt()==0);

  pls->push_probe(0,0,0);
  CPPUNIT_ASSERT(pls->probe_cnt()==1);
  CPPUNIT_ASSERT(pls->get_probeId(0)==0);

  pls->push_probe(1,1,1);
  CPPUNIT_ASSERT(pls->probe_cnt()==2);
  CPPUNIT_ASSERT(pls->get_probeId(1)==1);

  pls->set_probeId(1,100);
  CPPUNIT_ASSERT(pls->get_probeId(1)==100);
  pls->set_probeApid(1,200);
  CPPUNIT_ASSERT(pls->get_probeApid(1)==200);
  
  delete pls;
}
