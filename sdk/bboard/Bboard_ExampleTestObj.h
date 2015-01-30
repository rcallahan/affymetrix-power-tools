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
// /nfs/ss11/harley/Apt/work/apt2-1/affy/sdk/bboard/Bboard_ExampleTestObj.h ---
//
// $Id: Bboard_ExampleTestObj.h,v 1.1 2009-10-31 23:24:34 harley Exp $
//

#ifndef _BBOARD_EXAMPLETESTOBJ_H_
#define _BBOARD_EXAMPLETESTOBJ_H_

#include <string>

// To get this class to play with the GC system,
// it needs to be defined where BBT__NEWTYPE is.
class ExampleTestObj {
public:
  std::string m_msg_goodbye;

  ExampleTestObj() {
    printf("TestExampleObject(): hello:\n");
  };
  ExampleTestObj(const std::string& msg_hello) {
    printf("TestExampleObject(): hello: '%s'\n",msg_hello.c_str());
  };
  ExampleTestObj(const std::string& msg_hello, const std::string& msg_goodbye) {
    printf("TestExampleObject(): hello: '%s'\n",msg_hello.c_str());
    m_msg_goodbye=msg_goodbye;
  };
  ~ExampleTestObj() {
    printf("TestExampleObject(): goodbye: '%s'\n",m_msg_goodbye.c_str());
  };
};

#endif // _BBOARD_EXAMPLETESTOBJ_H_
