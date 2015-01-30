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
// ~/Affy-projects/apt2/trunk/PN_TestVecStats.cpp ---
//
// $Id: PN_TestVecStat.cpp,v 1.1 2009-10-28 18:25:58 harley Exp $
//

#include "PN_TestVecStat.h"
//
#include <vector>


// PN_TestVecStat is a collection ("Engine") of two Pnodes.
// This node computes the size and sum of a vector called "vec"
// using the two nodes.
// These Pnodes are independent, so it doesnt make a difference which
// order they are run in.

PN_TestVecStat::PN_TestVecStat() {
  // say hello.
  printf("new PN_TestVecStat()  %p\n",this);
  // build up our Pnode tree.
  // in this case we have two children which do the work.
  // When run we will compute the size...
  addPnode(new PN_TestVecStatSize());
  // ...and the sum of the vector.
  addPnode(new PN_TestVecStatSum());
}

// This gets run before the child nodes.
// It creates some default test data if none is given in the input BB.
// (the BB is empty on input for our testing purposes.)
AptErr_t PN_TestVecStat::doRunNodePre(Bboard* bb)
{
  void* voidptr;
  std::vector<int>* vecptr;

  // Test to see if the BB has data already,
  // if so, dont generate test data.
  // note the check for the name AND type of data.
  if (bb->getPtrIfType("vec",voidptr,BBT_VEC_INT)==APT_OK) {
    return APT_OK;
  }

  // no vector data was found on the blackboard, generate some.

  // poof! a new vec.
  vecptr=new std::vector<int>();
  // which we stick into the backboard, which now owns it.
  // The programer doenst delete it, the BB will on exit, or if we ask it to.
  bb->set("vec",vecptr);
  
  // add "0..9" to the vector.
  for (int i=0;i<10;i++) {
    vecptr->push_back(i);
  }
  return APT_OK;
}

AptErr_t PN_TestVecStat::doRunNodePost(Bboard* bb)
{
  // Now we retreive the results and print them out.
  // the results are left in "vec-size" and "vec-sum"
  AptErr_t rv;
  int tmp_vec_size;
  int tmp_vec_sum;
  int tmp_foo_bar;
  //
  rv=bb->get("vec-size",&tmp_vec_size);
  printf("== PN_TestVecStat(): 'vec-size'=%d  (rv=%d)\n",tmp_vec_size,rv);
  rv=bb->get("vec-sum",&tmp_vec_sum);
  printf("== PN_TestVecStat(): 'vec-sum'=%d  (rv=%d)\n",tmp_vec_sum,rv);
  // this should have an error code...
  rv=bb->get("foo-bar",&tmp_foo_bar);
  printf("== PN_TestVecStat(): 'foo-bar'=%d  (rv=%d) (should be !=0)\n",tmp_foo_bar,rv);
  //
  return APT_OK;
}
//////////

PN_TestVecStatSize::PN_TestVecStatSize() {
  // hello!
  printf("new PN_TestVecStatSize()  %p\n",this);
}

// This is an example of what most "PN_" nodes will do.
// Programmers will write their own "doRunNode" which will carry out the
// computations they want; In this case we figure out the size and store that on the BB.

AptErr_t PN_TestVecStatSize::doRunNode(Bboard* bb)
{
  //
  AptErr_t rv;
  std::vector<int>* vecptr;

  // bb->get() is an overloaded function.  It written to check the type vecptr. 
  // vecptr will be NULL unless the name is found AND the types match.
  rv=bb->get("vec",&vecptr);

  // the vector was found on the BB, do what we needed to do.
  if (vecptr!=NULL) {
    // which is store the size of the vector on the BB.
    // (Not all that exciting, but the journey starts with the first step.)
    // like "->get()", "->set()" is overloaded so that it stores the type of the data
    // If we need the long form, it is "->setPtrAndType(name,ptr,type)".
    bb->set("vec-size",(int)vecptr->size());
  }
  return APT_OK;
}

//////////

PN_TestVecStatSum::PN_TestVecStatSum() {
  // hello!
  printf("new PN_TestVecStatSum()  %p\n",this);
}

AptErr_t PN_TestVecStatSum::doRunNode(Bboard* bb)
{
  // just like the above...
  std::vector<int>* vecptr;
  bb->get("vec",&vecptr);

  // ...except we sum the contents.
  if (vecptr!=NULL) {
    int sum=0;
    for (int i=0;i<vecptr->size();i++) {
      sum+=(*vecptr)[i];
    }
    // and store the results.
    bb->set("vec-sum",sum);
  }
  return APT_OK;
}
