////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

#include "util/AffxBinaryFile.h"
#include "util/AffxString.h"
#include "util/CPPTest/Setup.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
using namespace std;


/**
 * @class AffxBinaryFileTest
 * @brief cppunit class for testing AffxBinaryFile functions.
 * last change by vliber on 05/26/09
 */

class AffxBinaryFileTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(AffxBinaryFileTest);
    CPPUNIT_TEST(openTest);
	CPPUNIT_TEST(readwriteTest);
	CPPUNIT_TEST(readwriteLittleEndianTest);
	CPPUNIT_TEST(readwriteBigEndianTest);    
    CPPUNIT_TEST_SUITE_END();

    public:  
     void openTest();
     void readwriteTest();
     void readwriteLittleEndianTest();
     void readwriteBigEndianTest();
     void setUp() {
       if ( !Fs::dirExists(OUTPUT) ) {
         Fs::mkdirPath(OUTPUT);
       }
   };
};

CPPUNIT_TEST_SUITE_REGISTRATION(AffxBinaryFileTest );

void AffxBinaryFileTest::openTest()
{
	cout<<endl;
	Verbose::out(1, "***AffxBinaryFileTest testcases***");
    Verbose::out(1, "AffxBinaryFileTest::openTest");
	AffxBinaryFile f1,f2,f3,f4,f5;
	CPPUNIT_ASSERT(f1.open(OUTPUT+"/test2.bin",AffxBinaryFile::SAVE));
	f1.close();
        Fs::rm(OUTPUT+"/test2.bin");
	CPPUNIT_ASSERT(!f2.open(INPUT+"/1lqcdfpsi_winTest1.bcdf",AffxBinaryFile::LOAD));
	f2.close();
	CPPUNIT_ASSERT(f3.open(INPUT+"/1lqcdfpsi_win.bcdf",AffxBinaryFile::LOAD));
	f3.close();
	CPPUNIT_ASSERT(f4.open(OUTPUT+"/1lqcdfpsi_winTest2.bcdf",AffxBinaryFile::APPEND));
	f4.close();
	CPPUNIT_ASSERT(f5.open(INPUT+"/1lqcdfpsi_win.bcdf",AffxBinaryFile::APPEND));
	f5.close();
	Fs::rm(OUTPUT+"/1lqcdfpsi_winTest2.bcdf");
}

void AffxBinaryFileTest::readwriteTest()
{
	  Verbose::out(1, "AffxBinaryFileTest::readwriteTest");
	  AffxBinaryFile f;
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary1.bin",AffxBinaryFile::SAVE));
	  f.writeInt(123456789);
	  f.writeUnsignedInt(1234);
	  f.writeFloat(1234.456f);
	  f.writeShort(13);
      f.writeUnsignedShort(2);
      f.writeChar('c');
	  f.writeUnsignedChar('G');
	  AffxString s="acgtacgt";
	  f.writeString(s);
	  f.writeString("ACGTACGT",8);
	  f.writeString("AaCcGgTt",8);
      f.close();
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary1.bin",AffxBinaryFile::LOAD));
      CPPUNIT_ASSERT(f.readInt()==123456789);
	  f.setOffset(0);
	  CPPUNIT_ASSERT(f.readLittleEndianInt()==123456789);
	  CPPUNIT_ASSERT(f.readUnsignedInt()==1234);
	  CPPUNIT_ASSERT(f.readFloat()==1234.456f);
	  CPPUNIT_ASSERT(f.readShort()==13);
	  CPPUNIT_ASSERT(f.readUnsignedShort()==2);
	  CPPUNIT_ASSERT(f.readChar()==99);
	  CPPUNIT_ASSERT(f.readUnsignedChar()==71);
	  CPPUNIT_ASSERT(f.readString()=="acgtacgt");
	  CPPUNIT_ASSERT(f.readString(8)=="ACGTACGT");
	  char out[9];
	  CPPUNIT_ASSERT(f.readString(out,8)=="AaCcGgTt");
	  CPPUNIT_ASSERT(f.readString(out,7)=="AaCcGgT");
	  f.close();
	  //real file
	  CPPUNIT_ASSERT(f.open(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP",AffxBinaryFile::LOAD));
	  CPPUNIT_ASSERT(f.readUnsignedInt()==65);
	  CPPUNIT_ASSERT(f.readUnsignedInt()==1);
	  CPPUNIT_ASSERT(f.readUnsignedShort()==1164);
	  CPPUNIT_ASSERT(f.readUnsignedShort()==1164);
	  f.setOffset(24);
	  CPPUNIT_ASSERT(f.readString()=="GeneChip.CallGEBaseCall.1");
	  f.setOffset(0);
	  f.advanceOffset(10);
	  CPPUNIT_ASSERT(f.getOffset()==10);
	  CPPUNIT_ASSERT(f.readUnsignedShort()==1164);
	  CPPUNIT_ASSERT(f.getOffset()==12);
	  f.setOffset(885);
      CPPUNIT_ASSERT(f.readInt()==13);
	  f.close();
}
void AffxBinaryFileTest::readwriteLittleEndianTest()
{
	  Verbose::out(1, "AffxBinaryFileTest::readwriteLittleEndianTest");
	  AffxBinaryFile f;
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary2.bin",AffxBinaryFile::SAVE));
	  f.writeLittleEndianInt(123456789);
	  f.writeLittleEndianUnsignedInt(1234);
	  f.writeLittleEndianShort(13);
      f.writeLittleEndianUnsignedShort(2);
	  AffxString s="acgtacgt";
	  f.writeLittleEndianString(s);
	  f.writeLittleEndianFloat(1234.456f);
	  
	  f.close();
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary2.bin",AffxBinaryFile::LOAD));
      CPPUNIT_ASSERT(f.readLittleEndianInt()==123456789);
	  CPPUNIT_ASSERT(f.readUnsignedInt()==1234);
	  CPPUNIT_ASSERT(f.readLittleEndianShort()==13);
	  CPPUNIT_ASSERT(f.readLittleEndianUnsignedShort()==2);
	  CPPUNIT_ASSERT(f.readLittleEndianString()=="acgtacgt");
	  #ifdef WIN32
	  CPPUNIT_ASSERT(f.readLittleEndianFloat()==1234.456f);
      #else
	  //CPPUNIT_ASSERT(f.readLittleEndianFloat()==1234.456f);//todo vliber
	  #endif
	  f.setOffset(4);
	  CPPUNIT_ASSERT(f.readLittleEndianUnsignedInt()==1234);
	  f.close();
	  //real file
	  CPPUNIT_ASSERT(f.open(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP",AffxBinaryFile::LOAD));
	  f.setOffset(0);
	  CPPUNIT_ASSERT_EQUAL(65,f.readLittleEndianInt());
	  CPPUNIT_ASSERT_EQUAL(1,f.readLittleEndianInt());
	  CPPUNIT_ASSERT_EQUAL((short)1164,f.readLittleEndianShort());
	  CPPUNIT_ASSERT_EQUAL((short)1164,f.readLittleEndianShort());
	  CPPUNIT_ASSERT(f.getOffset()==12);
	  f.setOffset(24);
	  CPPUNIT_ASSERT(f.readLittleEndianString()=="GeneChip.CallGEBaseCall.1");
	  f.setOffset(885);
      CPPUNIT_ASSERT(f.readLittleEndianInt()==13);
	  f.close();
}
void AffxBinaryFileTest::readwriteBigEndianTest()
{
	  Verbose::out(1, "AffxBinaryFileTest::readwriteBigEndianTest");
	  AffxBinaryFile f;
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary3.bin",AffxBinaryFile::SAVE));
	  f.writeBigEndianInt(123456789);
	  f.writeBigEndianUnsignedInt(1234);
	  f.writeBigEndianShort(13);
      f.writeBigEndianUnsignedShort(2);
      AffxString s="acgtacgt";
	  f.writeBigEndianString(s);
	  f.writeBigEndianFloat(1234.456f);
	  f.close();
	  CPPUNIT_ASSERT(f.open(OUTPUT+"/testBinary3.bin",AffxBinaryFile::LOAD));
    CPPUNIT_ASSERT(f.readBigEndianInt()==123456789);
	  CPPUNIT_ASSERT(f.readBigEndianUnsignedInt()==1234);
	  CPPUNIT_ASSERT(f.readBigEndianShort()==13);
	  CPPUNIT_ASSERT(f.readBigEndianUnsignedShort()==2);
	  CPPUNIT_ASSERT(f.readBigEndianString()=="acgtacgt");
	  #ifdef WIN32
	  CPPUNIT_ASSERT(f.readBigEndianFloat()==1234.456f);
      #else
	  //CPPUNIT_ASSERT(f.readBigEndianFloat()==1234.456f); //todo vliber
      #endif
	  f.setOffset(16);
	  CPPUNIT_ASSERT(f.readBigEndianString(7)=="acg");
	  f.setOffset(0);
	  f.advanceOffset(16);
	  CPPUNIT_ASSERT(f.readBigEndianString(5)=="a");
	  f.advanceOffset(16);
	  CPPUNIT_ASSERT(f.readBigEndianString(4)=="");
	  f.close();
}





