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
// ~/Affy/apt2/trunk/test-pgm.cpp ---
//
// $Id: test-bboard.cpp,v 1.8 2009-11-05 20:41:39 harley Exp $
//

//
#include "bboard/Bboard.h"
#include "bboard/BboardObj.h"
#include "bboard/dao/Dao.h"
#include "bboard/pnode/Pnode.h"

//
#include "bboard/extra/BioSpecies.h"
#include "util/AptVersionInfo.h"
#include "util/PgOptions.h"

//
#include <cassert>
#include <stdexcept>
#include <string.h>
#include <stdlib.h>

//////////

void example1_sub1(Bboard* bb)
{
  Bboard bad_bb;
  //bb->set("foo",&bad_bb);
}


void example1_main()
{
  // create a top-level BB; call it 'main';
  Bboard* bb_main=new Bboard("main");

  example1_sub1(bb_main);

  delete bb_main;
}

//////////

// Examples of error condtions.
// If you uncomment the lines marked as "BAD" the program should abort.
void error_examples()
{
  Bboard* bb=new Bboard("bad_example1");

  // Dont allow loops to form. 
  // (easily at least, you can still make them in multiple steps 
  // or BB1 -> BB2 -> BB1
  //bb->set("loop",bb);

  // if we attempt to store a reference of something on the 
  // stack to a BB we will die!
  Bboard bad_stack_bb;
  // bb->set("foo",&bad_bb);

  // With the
  // char* ptr=(char*)malloc(10);
  // bb->set("char*",ptr);
  //
  bb->dump();
  //
  delete bb;
}

//////////

// What follows are a bunch of tests
// which arent that 

void test_daopath_1()
{
  FsPath path("foo.bar");
  path.dump();
}

void test_refvarval_1()
{
  BboardBoxRef ref1;
  ref1->setData(1);
  ref1->setData(1.0);
  ref1->setData("one");
  ref1=NULL;
  ref1->setData(2);
  ref1->setData(2.0);
  ref1->setData("two");
  ref1=NULL;
}

void test_bboard_1()
{
  AptErr_t rv;
  
  Bboard* bb1=new Bboard("one");
  //
  int foo_int;
  bb1->set("foo_int",1);
  bb1->get("foo_int",&foo_int);
  double foo_double;
  bb1->set("foo_double",2.0);
  bb1->get("foo_double",&foo_double);

  // abortOnErr==0 and check this fails.
  rv=bb1->get("foo_double",&foo_int,0);
  assert(rv!=APT_OK);
  
  //
  bb1->set("foo_todelete",3.0);

  bb1->dump();
  std::string s = "holla";
  std::string t;
  bb1->set("one", s);
  bb1->get("one", &t);
  printf("first is: %s second is: %s\n", s.c_str(), t.c_str());
  assert(t == s);
  
  Bboard* bb2=new Bboard("two");
  bb2->copyFrom(bb1);
  bb2->dump();

  //
  bb1->unbind("foo_todelete");
  bb2->unbind("foo_todelete");
  // should be a message about it being deleted.

  bb1->dump();
  bb2->dump();

  //
  printf("### refcnts should be the same.");
  bb2->copyFrom(bb1);
  bb2->dump();
  bb2->copyFrom(bb1);
  bb2->dump();
  bb2->copyFrom(bb1);
  bb2->dump();

  delete bb2;
  delete bb1;
}

void test_bboard_2() {
  Bboard* bb1=new Bboard("bb1");

  //
  bb1->arrayResize(10);
  
  if (0) {
    for (int i=0;i<bb1->arraySize();i++) {
      bb1->set(i,i);
    }
  }
  bb1->arrayResize(0);
  bb1->arrayResize(10);

  if (1) {
    // use daoPath here as 'Dao_File' leaves behind garbarge.
    // (We havent finished it yet.)
    for (int i=0;i<bb1->arraySize();i++) {
      FsPath* df=new FsPath();
      bb1->set(i,df);
    }
  }

  if (1) {
    for (int i=0;i<bb1->arraySize();i++) {
      Bboard* bb=new Bboard();
      bb1->set(i,bb);
    }
  }

  //
  bb1->dump();
  bb1->arrayResize(0);
  bb1->dump();

  delete bb1;
}

//
void test_bboard_3() {
  Bboard* bb1=new Bboard("bb1");
  
  BioSpecies* ptr1=BioSpecies::speciesFromName("human");
  bb1->set("biospecies",ptr1);
  bb1->findBoxPtr("biospecies")->setDataOwned(false);

  //
  delete bb1;
  // delete it ourselves
  delete ptr1;

}

void test_biotypes()
{
  // try adding different biotypes to a blackboard.
  Bboard* bb1=new Bboard("bb1");

  //
  BioSpecies* ptr1=BioSpecies::speciesFromName("human");
  // the long form.
  // dont need to include "BioSpecies.h"
  // bb1->setPtrAndType("biospecies",ptr,BBT_BIOSPECIES);
  // the short form, the type comes from the ptr.
  bb1->set("biospecies",ptr1);

  //
  ChipLayout* ptr2=new ChipLayout();
  bb1->set("layout",ptr2);

  bb1->dump();
  delete bb1;
}

void test_parent()
{
  Bboard* bb1=new Bboard("bb1");
  Bboard* bb2=bb1->makeChildBboard("bb2");

  int val;
  bb1->set("one",1);

  //
  bb1->get("one",&val);
  assert(val==1);
  bb1->get("one",&val);
  assert(val==1);

  //
  bb2->set("one",2);

  //
  bb1->get("one",&val);
  assert(val==1);
  bb2->get("one",&val);
  assert(val==2);
  //
  delete bb1;
}

void test_vector() 
{
  Bboard* bb1=new Bboard("bb1");
  std::vector<int>* vec1;

  vec1=new std::vector<int>();
  bb1->set("vec1",vec1);
  bb1->get("vec1",&vec1);
  vec1->resize(10);
  (*vec1)[0]=100;

  delete bb1;
}

//////////

// A helper function to "test_malloc" below.
void test_malloc_1(Bboard* bb,int idx)
{
  void* mem;
  const int mem_len=100;

  // thow an error on the tenth call for testing.
  if (idx==10) {
    throw std::runtime_error("Trigger idx reached.");
  }

  // allocate some memory from the BB.
  // dont deallocate it as we want the BB to clean it up.
  mem=bb->Malloc(mem_len);
  printf("test_malloc_1: %d: %p\n",idx,mem);
  memset(mem,idx,mem_len);

  // this mem is allocated and freed right here
  // as a test of the "->Free()" call.
  mem=bb->Malloc(mem_len);
  bb->Free(mem);
}

void test_malloc()
{
  Bboard* bb1=new Bboard("bb1");
  Bboard* bb2=bb1->makeChildBboard("child");

  // This try/catch is a test of a BB cleaning up allocations
  // in the face of a nonlocal exit. test_malloc_1 will throw
  // an error on the 10th call.
  try {
    for (int i=0;i<20;i++) {
      test_malloc_1(bb2,i);
    }
  }
  catch (std::exception& err) {
    printf("test_malloc: caught error: '%s'\n",err.what());
  }
  
  // Even with the error, the BB reclaims 
  // the resources of bb1 AND bb2
  delete bb1;
}

// Test base types.
void test_io_1()
{
  Bboard* bb1=new Bboard("bb1");
  
  bb1->set("double100.0",100.0);
  bb1->set("double200.0",200.0);
  bb1->set("float1.0",1.0);
  bb1->set("float2.0",2.0);
  bb1->set("int1",1);
  bb1->set("int2",2);
  bb1->set("int3",3);
  bb1->set("string1","string1string1");
  bb1->set("string2","string2string2");

  DaoUtil::writeToFile("test-io-1-1.bb",bb1);

  Bboard* bb2=new Bboard("bb2");
  DaoUtil::readFromFile("test-io-1-1.bb",bb2);
  bb2->dump();
  DaoUtil::writeToFile("test-io-1-2.bb",bb2);

  //
  delete bb1;
  delete bb2;
}

void test_io_2()
{
  Bboard* bb1_1=new Bboard("bb1_1");
  Bboard* bb1_2=new Bboard("bb1_2");
  
  bb1_1->set("int1",1);
  bb1_1->set("int2",2);
  bb1_1->set("int3",3);
  bb1_1->set("bb1_2",bb1_2);
  bb1_1->dump();

  bb1_2->set("int4",4);
  bb1_2->set("int5",5);
  bb1_2->set("int6",6);
  bb1_2->dump();

  DaoUtil::writeToFile("test-io-2-1.bb",bb1_1);

//  Bboard* bb2_1=new Bboard("bb2_1");
//  bb2->readFromTsvFile("test-io1.bb");
//  bb2->dump();
//  bb2->writeToTsvFile("test-io-2",".bb");

  Bboard* bb2_1=new Bboard("bb2_1");
  DaoUtil::readFromFile("test-io-2-1.bb",bb2_1);

  // does both bb1_1 AND bb1_2;
  delete bb1_1;
  delete bb2_1;
  // does both bb1_1
}

void test_io_3()
{
  PgOptions* pgopts=new PgOptions();
  pgopts->defineOption("", "string-opt", PgOpt::STRING_OPT,
                       "A string option.",
                       "a test string");
  pgopts->defineOption("", "double-opt", PgOpt::DOUBLE_OPT,
                       "A doubple option.",
                       "100.0");
  pgopts->defineOption("", "int-opt", PgOpt::INT_OPT,
                       "A int option.",
                       "200");

  Bboard* bb=new Bboard("bb");
  bb->set("foo","bar");

  bb->set("pgopts",pgopts);
  DaoUtil::writeToFile("test-io-3.bb",bb);

  //
  bb->set("foo","some other value which will be replaced.");

  //
  DaoUtil::readFromFile("test-io-3.bb",bb);

  // bb->dump();

  // dont do this! Remember what is given to the BB,
  // is owned by the BB.
  // delete pgopts; // BAD!
  delete bb;
}

//////////

void test_pnode_1()
{
  Pnode* pnode1=new Pnode("pnode1");
  Bboard* bb1=new Bboard("bb1");

  pnode1->doDump(bb1,"");
  pnode1->doRun(bb1);

  delete bb1;
  delete pnode1;
}

//////////

void test_dao_1()
{
  Dao_File* df=new Dao_File();
  df->create("test-dao1.tsv",DAO_CREATE);
  Dao_Group* dg=df->createGroup("group",0);
  Dao_Table* dt=dg->createTable("table",0);
  
  dt->addHeader("foo","bar");

  dt->defineColumn(0,0,"key"  ,DAO_STRING);
  dt->defineColumn(0,1,"type" ,DAO_STRING);
  dt->defineColumn(0,2,"type" ,DAO_STRING);

  dt->endHeaders();

  dt->close();
  delete dt;
  dg->close();
  delete dg;
  df->close();
  delete df;
}

void test_dao_2()
{
  Bboard* bb1=new Bboard("test_dao_2");

  // @todo fix this when Dao_File has its own refcnting.
  //Dao_File* dao_f=new Dao_File();
  //bb1->set("file",dao_f);

  FsPath* dao_p=new FsPath("foo.bar");
  bb1->set("path",dao_p);

  bb1->dump();
  delete bb1;
}

//////////

int main(int argc,const char** argv)
{
  //
  error_examples();
  //
  example1_main();
  //
  test_refvarval_1();
  test_bboard_1();
  test_bboard_2();
  test_bboard_3();
  //
  //test_biotypes();
  //
  test_parent();
  test_vector();
  test_malloc();
  test_pnode_1();
  //
  test_dao_1();
  test_dao_2();
  //
  test_io_1();
  test_io_2();
  test_io_3();

  //
  return 0;
}
