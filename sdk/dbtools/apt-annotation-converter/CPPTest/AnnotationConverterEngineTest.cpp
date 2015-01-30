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

//
#include "../AnnotationConverterEngine.h"

//
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
#include <string> 
#include <vector> 

using namespace std;

/**
 * @class AnnotationConverterEngineTest
 * @brief cppunit class for testing AnnotationConverterEngineTest functions.
 */

class AnnotationConverterEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(AnnotationConverterEngineTest);   
  CPPUNIT_TEST ( CheckOptionsTest );
  CPPUNIT_TEST_SUITE_END();

public:  
  void CheckOptionsTest();  
};

CPPUNIT_TEST_SUITE_REGISTRATION(AnnotationConverterEngineTest);

// check the options
void AnnotationConverterEngineTest::CheckOptionsTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("AnnotationConverterEngineTest","checkOptionsTest");
	AnnotationConverterEngineTest engine;
	
}
