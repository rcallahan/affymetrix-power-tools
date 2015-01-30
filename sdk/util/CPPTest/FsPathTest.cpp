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

/**
 * @file   FsPathTest.cpp
 * @author
 * @date
 *
 * @brief
 */

//
#define __STDC_LIMIT_MACROS
//
#include "util/FsPath.h"
#include "util/Fs.h"
#include "util/CPPTest/Setup.h"
//
#include <iostream>
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

//
class FsPathTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( FsPathTest );
  //
  CPPUNIT_TEST( test_constructor );
  CPPUNIT_TEST( test_assign_1 );
  CPPUNIT_TEST( test_assign_2 );
  CPPUNIT_TEST( test_internalnames );
  CPPUNIT_TEST( test_ext2fmt );
  CPPUNIT_TEST( test_getPathFor );
  CPPUNIT_TEST( test_pushpop );
  CPPUNIT_TEST( test_dirnames );
  CPPUNIT_TEST( test_mkdir );
  CPPUNIT_TEST( test_iswriteabledir );
  CPPUNIT_TEST( test_touch );
  CPPUNIT_TEST( test_listdir );
  CPPUNIT_TEST( test_samevolume );
  CPPUNIT_TEST( test_freediskspace );
  CPPUNIT_TEST( test_win_mkdirpath );
  CPPUNIT_TEST( test_unix_mkdirpath );
  CPPUNIT_TEST( test_setDirPath );
  //CPPUNIT_TEST( test_ );
  //
  CPPUNIT_TEST_SUITE_END();

  void test_constructor();
  void test_assign_1();
  void test_assign_2();
  void test_internalnames();
  void test_ext2fmt();
  void test_getPathFor();
  void test_pushpop();
  void test_dirnames();
  void test_mkdir();
  void test_iswriteabledir();
  void test_touch();
  void test_listdir();
  void test_samevolume();
  void test_freediskspace();
  void test_win_mkdirpath();
  void test_unix_mkdirpath();
  void test_setDirPath();
  //void test_();
  void test_1();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( FsPathTest );

void FsPathTest::test_constructor()
{
  FsPath p1;
  FsPath p2(p1);
  FsPath p3("");
  //
  CPPUNIT_ASSERT(p1.asString()=="");
  CPPUNIT_ASSERT(p2.asString()=="");
  CPPUNIT_ASSERT(p3.asString()=="");

  p1.setPath("foo");
  FsPath p4(p1);
  CPPUNIT_ASSERT(p1==p4);

  // this is for bamboo machines only... shouldnt have this
#if 0
  std::string absolute_drive_path = "e:/apt/regression-data/probeset-genotype/Mapping250K_Sty/doChpFiles/chp/NA06985_B01_Sty_Plate1.brlmm.chp";

  if ( Fs::fileExists(absolute_drive_path ) && Fs::isReadable(absolute_drive_path))  {
	  
	  FsPath chpPath(absolute_drive_path);
	  CPPUNIT_ASSERT(chpPath.isReadable() );
  }
#endif

}

void FsPathTest::test_assign_1()
{
  FsPath p1;
  CPPUNIT_ASSERT(p1==p1);
  CPPUNIT_ASSERT(p1.empty()==true);

  p1.setPath("foo1.bar1");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="foo1.bar1");
  CPPUNIT_ASSERT(p1.getBaseName()=="foo1.bar1");
  CPPUNIT_ASSERT(p1.empty()==false);

  FsPath p2=p1;
  CPPUNIT_ASSERT(p1==p2);

  p2.setFileExt("bar2");
  //p2.dump();
  CPPUNIT_ASSERT(p2.empty()==false);
  CPPUNIT_ASSERT(!(p1==p2));
  CPPUNIT_ASSERT(p2.getBaseName()=="foo1.bar2");
  CPPUNIT_ASSERT(p2.asString()=="foo1.bar2");

  p2.setFileExt("bar3");
  //p2.dump();
  CPPUNIT_ASSERT(!(p1==p2));
  CPPUNIT_ASSERT(p2.asString()=="foo1.bar3");
  CPPUNIT_ASSERT(p2.getBaseName()=="foo1.bar3");

  p2.clear();
  CPPUNIT_ASSERT(p2.empty()==true);

  p2.setPath("a","b","c");
  CPPUNIT_ASSERT(p2.asString()=="a/b.c");

  //
  p1.setPath("c:/dir/foo1.bar1");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="c:/dir/foo1.bar1");
  p2.clear();
  p2=p1;
  CPPUNIT_ASSERT(p2.getDirName()=="c:/dir");
  CPPUNIT_ASSERT(p2.getBaseName()=="foo1.bar1");
  
}

void FsPathTest::test_assign_2()
{
  FsPath p1;

  p1.clear();
  p1.setPath("directory/");
  CPPUNIT_ASSERT(p1.getDirName()=="directory");
  CPPUNIT_ASSERT(p1.getFileName()=="");
  CPPUNIT_ASSERT(p1.getFileExt()=="");

  p1.clear();
  p1.setDirName("directory");
  CPPUNIT_ASSERT(p1.getDirName()=="directory");
  CPPUNIT_ASSERT(p1.getFileName()=="");
  CPPUNIT_ASSERT(p1.getBaseName()=="");
  CPPUNIT_ASSERT(p1.getFileExt()=="");

  p1.clear();
  p1.setPath("file");
  CPPUNIT_ASSERT(p1.getDirName()=="");
  CPPUNIT_ASSERT(p1.getFileName()=="file");
  CPPUNIT_ASSERT(p1.getFileExt()=="");

  p1.clear();
  p1.setPath(".ext");
  CPPUNIT_ASSERT(p1.getDirName()=="");
  CPPUNIT_ASSERT(p1.getFileName()=="");
  CPPUNIT_ASSERT(p1.getFileExt()=="ext");
  CPPUNIT_ASSERT(p1.asString()==".ext");
}

void FsPathTest::test_internalnames()
{
  FsPath p1;
  // @todo - change to ./foo.bar:group1/table1 for consistency with hdf5
  p1.setPath("foo.bar:./group1/table1.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="foo.bar:./group1/table1.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="foo.bar");

  p1.setPath("dir1/dir2/foo.bar:./group1/table1.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="dir1/dir2/foo.bar:./group1/table1.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="dir1/dir2/foo.bar");

  p1.setPath("/foo.bar:/group2/table2.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="/foo.bar:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="/foo.bar");

  p1.setPath("dir1/dir2/foo.bar:/group2/table2.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="dir1/dir2/foo.bar:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="dir1/dir2/foo.bar");

  p1.setPath("dir1/dir2/:/group2/table2.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="dir1/dir2/:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="dir1/dir2/");

  p1.setPath("foo.bar:/group2/table2.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="foo.bar:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="foo.bar");

  p1.setPath("foo.bar:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asString()=="foo.bar:/group2/table2.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="foo.bar");
  CPPUNIT_ASSERT(p1.getBaseName()=="foo.bar");

  p1.setPath("foo:group2/table2");
  CPPUNIT_ASSERT(p1.asString()=="foo:group2/table2");

  p1.setPath(":table1.baz");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()==":table1.baz");
  CPPUNIT_ASSERT(p1.asUnixPath()=="");
  CPPUNIT_ASSERT(p1.getBaseName()=="");

  //
  p1.setPath("C:foo.bar:group1/table1");
  //std::cout << p1.asString() << "\n";
  CPPUNIT_ASSERT(p1.asString()=="C:foo.bar:group1/table1");
  CPPUNIT_ASSERT(p1.getDirDrive()=="C:");
  CPPUNIT_ASSERT(p1.getBaseName()=="foo.bar");
  CPPUNIT_ASSERT(p1.getInternalGroupName()=="group1");
  CPPUNIT_ASSERT(p1.getInternalTableName()=="table1");

  p1.setPath("C:/foo.bar:group1/table1");
  //std::cout << p1.asString() << "\n";
  CPPUNIT_ASSERT(p1.asString()=="C:/foo.bar:group1/table1");
  CPPUNIT_ASSERT(p1.getDirDrive()=="C:");
  CPPUNIT_ASSERT(p1.getBaseName()=="foo.bar");
  CPPUNIT_ASSERT(p1.getInternalGroupName()=="group1");
  CPPUNIT_ASSERT(p1.getInternalTableName()=="table1");
}

void FsPathTest::test_ext2fmt()
{
  // is the default set loaded?
  CPPUNIT_ASSERT(FsPath::ext2Fmt("tsv")==FsPath::FILEFMT_TSVFILE);
  CPPUNIT_ASSERT(FsPath::ext2Fmt("ref")==FsPath::FILEFMT_FILE5);
  CPPUNIT_ASSERT(FsPath::ext2Fmt("chp")==FsPath::FILEFMT_CALVIN);

  // is an unknown treated as none?
  CPPUNIT_ASSERT(FsPath::ext2Fmt("none")==FsPath::FILEFMT_NONE);

  // can we add and update them?
  FsPath::addExt2Fmt("test1",FsPath::FILEFMT_TSVFILE);
  CPPUNIT_ASSERT(FsPath::ext2Fmt("test1")==FsPath::FILEFMT_TSVFILE);
  FsPath::addExt2Fmt("test1",FsPath::FILEFMT_CALVIN);
  CPPUNIT_ASSERT(FsPath::ext2Fmt("test1")==FsPath::FILEFMT_CALVIN);
}

void FsPathTest::test_getPathFor()
{
  // should work everywhere but windows.
#ifndef _MSC_VER
  FsPath p1=FsPath::getPathFor("HOME");
  //p1.dump();
  CPPUNIT_ASSERT(p1.empty()==false);

  FsPath p2=FsPath::getPathFor("PWD");
  //p2.dump();
  CPPUNIT_ASSERT(p2.empty()==false);

  // TMP isnt set for bambooadmin
  // FsPath p3=FsPath::getPathFor("TMP");
  // //p3.dump();
  // CPPUNIT_ASSERT(p3.empty()==false);
#endif

  FsPath p4=FsPath::getPathFor("THISENVDOESNOTEXIST_123_456");
  CPPUNIT_ASSERT(p4.empty()==true);
}

void FsPathTest::test_pushpop()
{
  FsPath p1;

  //
  p1.setDirName("a/b/c");
  //p1.dump();
  CPPUNIT_ASSERT(p1.canPopDir()==true);
  CPPUNIT_ASSERT(p1.popDir()=="c");
  CPPUNIT_ASSERT(p1.canPopDir()==true);
  CPPUNIT_ASSERT(p1.popDir()=="b");
  CPPUNIT_ASSERT(p1.canPopDir()==true);
  CPPUNIT_ASSERT(p1.popDir()=="a");
  CPPUNIT_ASSERT(p1.canPopDir()==false);
  CPPUNIT_ASSERT(p1.popDir()=="");
  CPPUNIT_ASSERT(p1.canPopDir()==false);
  CPPUNIT_ASSERT(p1.popDir()=="");
  CPPUNIT_ASSERT(p1.empty()==true);

  FsPath abs(std::string("c:\temp"));
  CPPUNIT_ASSERT(abs.isAbsolute()==true);
  p1.setDirName("/a/b");
  CPPUNIT_ASSERT(p1.isAbsolute()==true);
  CPPUNIT_ASSERT(p1.canPopDir()==true);
  CPPUNIT_ASSERT(p1.popDir()=="b");
  CPPUNIT_ASSERT(p1.canPopDir()==true);
  CPPUNIT_ASSERT(p1.popDir()=="a");
  CPPUNIT_ASSERT(p1.canPopDir()==false);
  CPPUNIT_ASSERT(p1.popDir()=="/");
  CPPUNIT_ASSERT(p1.canPopDir()==false);
  CPPUNIT_ASSERT(p1.popDir()=="/");
  CPPUNIT_ASSERT(p1.isAbsolute()==true);

  p1.setDirName("/a/");
  CPPUNIT_ASSERT(p1.isAbsolute()==true);
  CPPUNIT_ASSERT(p1.popDir()=="a");
  CPPUNIT_ASSERT(p1.popDir()=="/");
  CPPUNIT_ASSERT(p1.isAbsolute()==true);

}

void FsPathTest::test_dirnames()
{
  FsPath p1;
  std::vector<std::string> vec;

  vec=p1.getDirNames();
  CPPUNIT_ASSERT(vec.size()==0);

  p1.setDirName("aaa/bbb/ccc");
  CPPUNIT_ASSERT(p1.getBaseName()=="");
  vec=p1.getDirNames();
  CPPUNIT_ASSERT(vec.size()==3);
  CPPUNIT_ASSERT(vec[0]=="aaa");
  CPPUNIT_ASSERT(vec[1]=="bbb");
  CPPUNIT_ASSERT(vec[2]=="ccc");

  p1.setDirName("a/b/c/");
  vec=p1.getDirNames();
  CPPUNIT_ASSERT(vec.size()==3);
  CPPUNIT_ASSERT(vec[0]=="a");
  CPPUNIT_ASSERT(vec[1]=="b");
  CPPUNIT_ASSERT(vec[2]=="c");

  p1.setDirName("/");
  vec=p1.getDirNames();
  CPPUNIT_ASSERT(vec.size()==1);
  CPPUNIT_ASSERT(vec[0]=="/");

  p1.setDirName("/a/b/c/");
  vec=p1.getDirNames();
  CPPUNIT_ASSERT(vec.size()==4);
  CPPUNIT_ASSERT(vec[0]=="/");
  CPPUNIT_ASSERT(vec[1]=="a");
  CPPUNIT_ASSERT(vec[2]=="b");
  CPPUNIT_ASSERT(vec[3]=="c");
}

void FsPathTest::test_mkdir()
{

  FsPath p1;

  // use a name which should be unique.
  // start with an "A" to be at the top of the dir.
  std::string tmp_dir_name="./A-FsPathTest-123-456-cpptest";
  p1.setDirName(tmp_dir_name);
  //p1.dump();

  // dont throw errors
  p1.setErrAbort(false);

  // make sure it isnt there.
  // we dont care if this rmdir has an error.
  p1.rmdir();
  // but it had better not be there.
  CPPUNIT_ASSERT(p1.exists()==false);
  CPPUNIT_ASSERT(p1.dirExists()==false);

  // try some mkdirs.
  p1.mkdir(true);
  CPPUNIT_ASSERT(p1.getErrNum()==APT_OK);
  p1.mkdir(true);
  CPPUNIT_ASSERT(p1.getErrNum()!=APT_OK); // already there
  p1.mkdir(false);
  CPPUNIT_ASSERT(p1.getErrNum()==APT_OK); // ok if already there

  //
  p1.rmdir(true);
  CPPUNIT_ASSERT(p1.getErrNum()==APT_OK);
  p1.rmdir(true);
  CPPUNIT_ASSERT(p1.getErrNum()!=APT_OK); // already gone
  p1.rmdir(false);
  CPPUNIT_ASSERT(p1.getErrNum()==APT_OK); // ok if already gone


//  p1.osRmdir();
//
//  //
//  p1.setDirName("test-mkdir-2");
//  p1.mkdir(false);
//  //
//  p1.setDirName("./test-mkdir-3");
//  p1.mkdir(false);
//
//  //
//  p1.setDirName("test-mkdir-4/a/b/c");
//  p1.mkdirPath(false);
}

void FsPathTest::test_touch()
{
  FsPath p1;
  p1.setErrAbort(false);

  p1.setPath("./A-FsPathTest-123-456.test");
  // dont care about this return value.
  p1.rm();
  // better not exist!
  CPPUNIT_ASSERT(p1.exists()==false);
  CPPUNIT_ASSERT(p1.fileExists()==false);
  CPPUNIT_ASSERT(p1.dirExists()==false);

  // make it exist.
  CPPUNIT_ASSERT(p1.touch()==APT_OK);
  CPPUNIT_ASSERT(p1.exists()==true);
  CPPUNIT_ASSERT(p1.fileExists()==true);
  CPPUNIT_ASSERT(p1.dirExists()==false);

  // make it go away.
  CPPUNIT_ASSERT(p1.rm()==APT_OK);
  CPPUNIT_ASSERT(p1.exists()==false);
  CPPUNIT_ASSERT(p1.fileExists()==false);
  CPPUNIT_ASSERT(p1.dirExists()==false);
}

void FsPathTest::test_iswriteabledir() {
#if 0
  FsPath p1;
  p1.setDirName(".");
  CPPUNIT_ASSERT(p1.osIsWriteableDir()==true);

  //  p1.setDirName("C:\\My Documents");
  p1.setDirName("C:\\Documents and Settings\\hgorre\\My Documents");
  //p1.dump();
  CPPUNIT_ASSERT(p1.osIsWriteableDir()==true);
#endif
}

void FsPathTest::test_listdir()
{
  FsPath p1;
  p1.setErrAbort(false);

  std::vector<std::string> names;
  p1.setDirName(".");

  p1.listDir(names);
  CPPUNIT_ASSERT(names.size()!=0);

  CPPUNIT_ASSERT(p1.listDir(names)==APT_OK);
  CPPUNIT_ASSERT(names.size()!=0);

  p1.setDirName("This/Dir/Does/Not/Exist");
  p1.listDir(names);
  CPPUNIT_ASSERT(names.size()==0);

  CPPUNIT_ASSERT(p1.listDir(names)!=APT_OK);
  CPPUNIT_ASSERT(names.size()==0);
}

void FsPathTest::test_samevolume()
{
  FsPath p1;

  p1.setDirName(".");
  CPPUNIT_ASSERT(p1.isSameVolume(".")==true);

  //CPPUNIT_ASSERT(p1.isSameVolume(p1)==true);
}

//////////

void FsPathTest::test_freediskspace()
{
  FsPath p1;
  p1.setDirName(".");
  CPPUNIT_ASSERT(p1.getFreeDiskSpace()>0);

  p1.setErrAbort(false);
  p1.setDirName("/thispathdoesnotexist");
  CPPUNIT_ASSERT(p1.getFreeDiskSpace()==-1);
}

//////////

void FsPathTest::test_win_mkdirpath()
{
#ifdef _MSC_VER
  FsPath p1;
  p1.setErrAbort(false);

  p1.setDirDrive("C:");
  p1.setDirName("/");
  CPPUNIT_ASSERT(p1.mkdirPath(false)==APT_OK);
  // should it return an error if the directory already exists?
  CPPUNIT_ASSERT(p1.mkdirPath(true)!=APT_OK);
#endif
}

///////////

void FsPathTest::test_unix_mkdirpath()
{
#ifndef _MSC_VER

  FsPath p_foo;
  FsPath p_foo_bar;

  // start with paths ./foo and ./foo/bar
  std::string foo_dir_name="./foo";
  p_foo.setDirName(foo_dir_name);
  std::string foo_bar_dir_name="./foo/bar";
  p_foo_bar.setDirName(foo_bar_dir_name);

  // dont throw errors
  p_foo.setErrAbort(false);
  p_foo_bar.setErrAbort(false);

  // make sure dirs aren't there.
  // we dont care if this rmdir has an error.
  p_foo_bar.rmdir();
  p_foo.rmdir();
  // but it had better not be there.
  CPPUNIT_ASSERT(p_foo.dirExists()==false);

  // try to make both dirs
  CPPUNIT_ASSERT(p_foo_bar.mkdirPath(true) == APT_OK);
  CPPUNIT_ASSERT(p_foo_bar.dirExists());
  CPPUNIT_ASSERT(p_foo.dirExists());

  // try just making last dir in path
  p_foo_bar.rmdir();			// remove ./foo/bar
  CPPUNIT_ASSERT(p_foo.dirExists());	// ./foo still exists
  CPPUNIT_ASSERT(p_foo_bar.dirExists() == false);
  CPPUNIT_ASSERT(p_foo_bar.mkdirPath(true) == APT_OK);
  CPPUNIT_ASSERT(p_foo_bar.getErrNum()==APT_OK);
  CPPUNIT_ASSERT(p_foo_bar.dirExists());

  // try making last dir again
  CPPUNIT_ASSERT(p_foo_bar.mkdirPath(true) != APT_OK); // error if exists
  CPPUNIT_ASSERT(p_foo_bar.mkdirPath(false) == APT_OK); // no error if exists
  
  // check works with errIfExists == false
  p_foo_bar.rmdir();			// remove ./foo/bar
  CPPUNIT_ASSERT(p_foo.dirExists());	// ./foo still exists
  CPPUNIT_ASSERT(p_foo_bar.dirExists() == false);
  CPPUNIT_ASSERT(p_foo_bar.mkdirPath(false) == APT_OK);
  CPPUNIT_ASSERT(p_foo_bar.dirExists());

  // cleanup
  p_foo_bar.rmdir();
  p_foo.rmdir();

#endif
} // test_unix_mkdirpath



//////////
void FsPathTest::test_1()
{
  FsPath p1;
  CPPUNIT_ASSERT(true);
}
void FsPathTest::test_setDirPath()
{
  FsPath p1;
  p1.setDirPath("foo1.bar1");
  //p1.dump();
  
  CPPUNIT_ASSERT(p1.asString()=="foo1.bar1");
  CPPUNIT_ASSERT(p1.getBaseName()=="foo1.bar1");
  CPPUNIT_ASSERT(p1.empty()==false);

  FsPath p2=p1;
  CPPUNIT_ASSERT(p1==p2);
  
  //
  p1.setDirPath("c:/dir/foo1.bar1");
  //p1.dump();
  CPPUNIT_ASSERT(p1.asString()=="c:/dir/foo1.bar1");
  p2.clear();
  p2=p1;
  CPPUNIT_ASSERT(p2.getDirName()=="c:/dir");
  CPPUNIT_ASSERT(p2.getBaseName()=="foo1.bar1");

}

