////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

//

#include "util/Verbose.h"
#include <util/Util.h>
#include <util/md5sum.h>
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
//
#ifndef WIN32
  #include <arpa/inet.h>
#endif

using namespace std;
using namespace affx;

/**
 * @class md5sumTest
 * @brief cppunit class for testing generation of MD5 hash codes (derived from test-md5sum.cpp)
 */
class md5sumTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( md5sumTest );
  //
  CPPUNIT_TEST( testIntegerMD5sum );
  CPPUNIT_TEST( testStringMD5sum );
  CPPUNIT_TEST( testFileMD5sum );
  //
  CPPUNIT_TEST_SUITE_END();

public:
  /** one-time utility function **/
  void generate_md5_reference_results();
  /** test case prototypes */
  void testIntegerMD5sum();
  void testStringMD5sum();
  void testFileMD5sum();
};

void md5sumTest::generate_md5_reference_results()
{
	 Verbose::out(1, "md5sumTest::generate_md5_reference_results");
  // this function generates input files used to validate MD5 hash codes for 
  // integer types validated with md5sum.exe downloaded from activestate.com

  #ifdef WIN32
    #pragma warning( push )
    #pragma warning( disable : 4996 )
  #endif

  FILE *f;

  uint8_t val_uint8 = 0;
  f = fopen("val_uint8_val=0.dat","w+");
  fwrite( (void *)&val_uint8, 1, sizeof(uint8_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint8_val=0.dat
  // 93b885adfe0da089cdf634904fd59f71 *val_uint8_val=0.dat

  val_uint8 = 1;
  f = fopen("val_uint8_val=1.dat","w+");
  fwrite( (void *)&val_uint8, 1, sizeof(uint8_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint8_val=1.dat
  // 55a54008ad1ba589aa210d2629c1df41 *val_uint8_val=1.dat

  val_uint8 = 123;
  f = fopen("val_uint8_val=123.dat","w+");
  fwrite( (void *)&val_uint8, 1, sizeof(uint8_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint8_val=123.dat
  // f95b70fdc3088560732a5ac135644506 *val_uint8_val=123.dat

  uint16_t val_uint16 = 0;
  uint16_t val_uint16_bigendian=htons(val_uint16);
  f = fopen("val_uint16_val=0.dat","w+");
  fwrite( (void *)&val_uint16_bigendian, 1, sizeof(uint16_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint16_val=0.dat
  // c4103f122d27677c9db144cae1394a66 *val_uint16_val=0.dat

  val_uint16 = 1;
  val_uint16_bigendian=htons(val_uint16);
  f = fopen("val_uint16_val=1.dat","w+");
  fwrite( (void *)&val_uint16_bigendian, 1, sizeof(uint16_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint16_val=1.dat
  // 441077cc9e57554dd476bdfb8b8b8102 *val_uint16_val=1.dat

  val_uint16 = 12345;
  val_uint16_bigendian=htons(val_uint16);
  f = fopen("val_uint16_val=12345.dat","w+");
  fwrite( (void *)&val_uint16_bigendian, 1, sizeof(uint16_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint16_val=12345.dat
  // 0a8005f5594bd67041f88c6196192646 *val_uint16_val=12345.dat

  uint32_t val_uint32 = 0;
  uint32_t val_uint32_bigendian=htonl(val_uint32);
  f = fopen("val_uint32_val=0.dat","w+");
  fwrite( (void *)&val_uint32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint32_val=0.dat
  // f1d3ff8443297732862df21dc4e57262 *val_uint32_val=0.dat

  val_uint32 = 1;
  val_uint32_bigendian=htonl(val_uint32);
  f = fopen("val_uint32_val=1.dat","w+");
  fwrite( (void *)&val_uint32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint32_val=1.dat
  // f1450306517624a57eafbbf8ed995985 *val_uint32_val=1.dat

  val_uint32 = 123456;
  val_uint32_bigendian=htonl(val_uint32);
  f = fopen("val_uint32_val=123456.dat","w+");
  fwrite( (void *)&val_uint32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_uint32_val=123456.dat
  // d747074027ce566f5dd8697337091787 *val_uint32_val=123456.dat

  int32_t val_int32 = 1;
  uint32_t val_int32_bigendian=htonl((uint32_t)val_int32);
  f = fopen("val_int32_val=1.dat","w+");
  fwrite( (void *)&val_int32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_int32_val=1.dat
  // f1450306517624a57eafbbf8ed995985 *val_int32_val=1.dat

  val_int32 = -1;
  val_int32_bigendian=htonl((uint32_t)val_int32);
  f = fopen("val_int32_val=-1.dat","w+");
  fwrite( (void *)&val_int32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_int32_val=-1.dat
  // a54f0041a9e15b050f25c463f1db7449 *val_int32_val=-1.dat

  val_int32 = 123456;
  val_int32_bigendian=htonl((uint32_t)val_int32);
  f = fopen("val_int32_val=123456.dat","w+");
  fwrite( (void *)&val_int32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_int32_val=123456.dat
  // d747074027ce566f5dd8697337091787 *val_int32_val=123456.dat

  val_int32 = -123456;
  val_int32_bigendian=htonl((uint32_t)val_int32);
  f = fopen("val_int32_val=-123456.dat","w+");
  fwrite( (void *)&val_int32_bigendian, 1, sizeof(uint32_t), f );
  fclose(f);

  // C:\apt-cppunit-dev\affy\sdk\util\CPPTest>md5sum val_int32_val=-123456.dat
  // 143d3d9511743a43b6f1de3e29b0e040 *val_int32_val=-123456.dat

  #ifdef WIN32
    #pragma warning( pop ) 
  #endif

  return;
}

void md5sumTest::testIntegerMD5sum()
{
	Verbose::out(1, "***md5sumTest testcases***");
  Verbose::out(1, "md5sumTest::testIntegerMD5sum");

  md5sum md5sum;
  std::string sum;

  // unsigned 8-bit integer test cases 

  md5sum.init();
  md5sum.update_nbo( (uint8_t)0 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "93b885adfe0da089cdf634904fd59f71" );

  md5sum.init();
  md5sum.update_nbo( (uint8_t)1 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "55a54008ad1ba589aa210d2629c1df41" );

  md5sum.init();
  md5sum.update_nbo( (uint8_t)123 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "f95b70fdc3088560732a5ac135644506" );

  // unsigned 16-bit integer test cases 

  md5sum.init();
  md5sum.update_nbo( (uint16_t)0 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "c4103f122d27677c9db144cae1394a66" );

  md5sum.init();
  md5sum.update_nbo( (uint16_t)1 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "441077cc9e57554dd476bdfb8b8b8102" );

  md5sum.init();
  md5sum.update_nbo( (uint16_t)12345 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "0a8005f5594bd67041f88c6196192646" );

  // unsigned 32-bit integer test cases 

  md5sum.init();
  md5sum.update_nbo( (uint32_t)0 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "f1d3ff8443297732862df21dc4e57262" );

  md5sum.init();
  md5sum.update_nbo( (uint32_t)1 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "f1450306517624a57eafbbf8ed995985" );

  md5sum.init();
  md5sum.update_nbo( (uint32_t)123456 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "d747074027ce566f5dd8697337091787" );

  // signed 32-bit integer test cases 

  md5sum.init();
  md5sum.update_nbo( (int32_t)0 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "f1d3ff8443297732862df21dc4e57262" );

  md5sum.init();
  md5sum.update_nbo( (int32_t)1 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "f1450306517624a57eafbbf8ed995985" );

  md5sum.init();
  md5sum.update_nbo( (int32_t)-1 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "a54f0041a9e15b050f25c463f1db7449" );

  md5sum.init();
  md5sum.update_nbo( (uint32_t)123456 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "d747074027ce566f5dd8697337091787" );

  md5sum.init();
  md5sum.update_nbo( (uint32_t)-123456 );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "143d3d9511743a43b6f1de3e29b0e040" );

}

void md5sumTest::testStringMD5sum() 
{
  Verbose::out(1, "md5sumTest::testStringMD5sum");

  md5sum md5sum;
  std::string sum;

  md5sum.init();
  md5sum.update( "" );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "d41d8cd98f00b204e9800998ecf8427e" );

  md5sum.init();
  md5sum.update( "a" );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "0cc175b9c0f1b6a831c399e269772661" );

  md5sum.init();
  md5sum.update( "The quick brown fox jumps over the lazy dog" );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "9e107d9d372bb6826bd81d3542a419d6" );

  md5sum.init();
  md5sum.update( "The quick brown fox jumps over the lazy dog." );
  md5sum.final(sum);
  CPPUNIT_ASSERT( sum == "e4d909c290d0fb1ca068ffaddf22cbd0" );
}

void md5sumTest::testFileMD5sum()
{
  std::string sum;
  md5sum md5sum;

  std::string hw_file="input/hello_world.txt";
  CPPUNIT_ASSERT(md5sum.ofFile(hw_file,sum)==0);
  printf("md5sum.ofFile('%s')=='%s'\n",hw_file.c_str(),sum.c_str());
  CPPUNIT_ASSERT(sum=="6f5902ac237024bdd0c176cb93063dc4");

  CPPUNIT_ASSERT(md5sum.ofFile(hw_file,sum,"my-salt")==0);
  printf("md5sum.ofFile('%s', 'my-salt')=='%s'\n",hw_file.c_str(),sum.c_str());
  CPPUNIT_ASSERT(sum=="8a0710454e151557b1744c85e55a8b9f");

  CPPUNIT_ASSERT(md5sum.ofFile("input/this_file_does_not_exist.",sum)!=0);
  CPPUNIT_ASSERT(sum=="");

}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( md5sumTest );

