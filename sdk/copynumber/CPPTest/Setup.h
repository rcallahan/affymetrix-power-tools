////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#ifndef _SETUP_H_
#define _SETUP_H_
/**
 * @file Setup.h
 *
 * @brief This header definition of input files location. 
 */
#define	INPUT  AffxString()+"../../regression-data/data/idata/copynumber-cyto/cppunit/input" 
#define	OUTPUT AffxString()+"./output" 

#include "util/AffxString.h"

class Setup {};  

#define POSITIVE_TEST(expression) \
    {\
        bool doThrow = Err::getThrowStatus(); \
        Err::setThrowStatus(true);       \
	try \
	{\
		Verbose::setLevel(-1); \
		expression; \
		Verbose::setLevel(2);\
		Verbose::out(1, AffxString() + "\tPositive test passed [" + #expression + "]");\
	} \
	catch(...) \
	{\
		Verbose::setLevel(2); \
		Verbose::out(1, AffxString() + "\tFATAL ERROR: Positive test failed (exception thrown) [" + #expression + "]");\
		CPPUNIT_ASSERT(false);\
	} \
        Err::setThrowStatus(doThrow); \
        }

#define NEGATIVE_TEST(expression, ExceptionType) \
    { \
        bool doThrow = Err::getThrowStatus(); \
        Err::setThrowStatus(true);            \
	try \
	{\
		Verbose::setLevel(-1); \
		expression; \
		Verbose::setLevel(2);\
		Verbose::out(1, AffxString() + "\tFATAL ERROR: Negative test failed (exception not thrown) [" + #expression + "]");\
		CPPUNIT_ASSERT(false);\
	} \
	catch(ExceptionType e) \
	{\
		Verbose::setLevel(2); \
		AffxString str(e.what()); \
		int iFindIndex = str.find("ERROR:"); \
		while (iFindIndex != -1) {str = str.substring(iFindIndex + 6); iFindIndex = str.find("ERROR:");} \
                iFindIndex = str.find_first_not_of(" "); \
                if (iFindIndex > 0) {str = str.substring(iFindIndex - 1);} \
		Verbose::out(1, AffxString() + "\tNegative test passed [" + #expression + "]\tMessage: " + str);\
	} \
        Err::setThrowStatus(doThrow); \
    }

#endif
