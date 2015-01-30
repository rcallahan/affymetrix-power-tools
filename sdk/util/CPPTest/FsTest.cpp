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
 * @file   FsTest.cpp
 * @author
 * @date
 *
 * @brief
 */

//
#include "util/Fs.h"
//
#include "util/CPPTest/Setup.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#ifndef _WIN32
#include <sys/stat.h>
#endif


//
class FsTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(FsTest);
  //
  CPPUNIT_TEST(test_unixify);
  CPPUNIT_TEST(test_mkdir);
  CPPUNIT_TEST(test_ensurewriteabledir);
  CPPUNIT_TEST(test_iswriteabledir);
  CPPUNIT_TEST(test_touch);
  CPPUNIT_TEST(test_listdir);
  CPPUNIT_TEST(test_samevolume);
  CPPUNIT_TEST(test_freediskspace);
  CPPUNIT_TEST(test_join);
  CPPUNIT_TEST(test_split);
  CPPUNIT_TEST(test_convertToUncPath);
  CPPUNIT_TEST(test_permissions);
  CPPUNIT_TEST(test_fileRemove);
  CPPUNIT_TEST(test_directoryReadable);
  CPPUNIT_TEST(test_fileRename);
  CPPUNIT_TEST(test_mustOpenToWrite);
  CPPUNIT_TEST(test_carefulClose);
  CPPUNIT_TEST(test_fileCopy);
  //
  CPPUNIT_TEST( test_trailingSlash );
  //
  CPPUNIT_TEST( test_isbinary );
  CPPUNIT_TEST( test_iscalvin );
  CPPUNIT_TEST( test_ishdf5 );
  //
  CPPUNIT_TEST_SUITE_END();

  void test_unixify();
  void test_mkdir();
  void test_ensurewriteabledir();
  void test_base_dir_ext();
  void test_iswriteabledir();
  void test_touch();
  void test_listdir();
  void test_samevolume();
  void test_freediskspace();
  void test_join();
  void test_split();
  void test_convertToUncPath();
  void test_permissions();
  void test_fileRemove();
  void test_directoryReadable();
  void test_mustOpenToWrite();
  void test_carefulClose();
  void test_fileCopy();
  void test_fileRename();
  //
  void test_trailingSlash();
  //
  void test_isbinary();
  void test_iscalvin();
  void test_ishdf5();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FsTest);

//////////

void FsTest::test_unixify()
{
  CPPUNIT_ASSERT(Fs::unixifyPath("") == "");
  CPPUNIT_ASSERT(Fs::unixifyPath("/") == "/");
  CPPUNIT_ASSERT(Fs::unixifyPath("a") == "a");
  CPPUNIT_ASSERT(Fs::unixifyPath("a/b") == "a/b");
  CPPUNIT_ASSERT(Fs::unixifyPath("a\\b") == "a/b");
  CPPUNIT_ASSERT(Fs::unixifyPath("c:/a/b") == "c:/a/b");
  CPPUNIT_ASSERT(Fs::unixifyPath("c:a/b") == "c:a/b");
  CPPUNIT_ASSERT(Fs::unixifyPath("c:a\\b") == "c:a/b");
  CPPUNIT_ASSERT(Fs::unixifyPath("a\\\\b") == "a//b");
  //
  CPPUNIT_ASSERT(Fs::windowifyPath("") == "");
  CPPUNIT_ASSERT(Fs::windowifyPath("/") == "\\");
  CPPUNIT_ASSERT(Fs::windowifyPath("a") == "a");
  CPPUNIT_ASSERT(Fs::windowifyPath("a/b") == "a\\b");
  CPPUNIT_ASSERT(Fs::windowifyPath("a\\b") == "a\\b");
  CPPUNIT_ASSERT(Fs::windowifyPath("a//b") == "a\\\\b");
}

void FsTest::test_join()
{
  // 2 args
  CPPUNIT_ASSERT(Fs::join("", "") == "");
  CPPUNIT_ASSERT(Fs::join("a", "") == "a");
  CPPUNIT_ASSERT(Fs::join("", "b") == "b");
  CPPUNIT_ASSERT(Fs::join("/", "") == "/");
  CPPUNIT_ASSERT(Fs::join("/", "b") == "/b");
  CPPUNIT_ASSERT(Fs::join("/a", "b") == "/a/b");
  //
  CPPUNIT_ASSERT(Fs::join("\\a", "b") == "/a/b");

  //
  CPPUNIT_ASSERT(Fs::join("c:", "abc") == "c:/abc");
  // drive letter example.
  //CPPUNIT_ASSERT(Fs::join("c:"+"abc")=="c:abc");

  // @todo: note that ".."s are preserved not cut out.
  //        add a method to cut out ".."s?

  // 3 args
  CPPUNIT_ASSERT(Fs::join("", "", "") == "");
  //
  CPPUNIT_ASSERT(Fs::join("a", "", "") == "a");
  CPPUNIT_ASSERT(Fs::join("", "b", "") == "b");
  CPPUNIT_ASSERT(Fs::join("", "", "c") == "c");
  CPPUNIT_ASSERT(Fs::join("a", "b", "c") == "a/b/c");
  //
  CPPUNIT_ASSERT(Fs::join("/", "b", "") == "/b");
  CPPUNIT_ASSERT(Fs::join("", "/", "c") == "/c");

  // what should happen here?
  // CPPUNIT_ASSERT(Fs::join("/","/b","")=="/b");
  // CPPUNIT_ASSERT(Fs::join("/a","/b","")=="/a/b");

  // 4 args
  CPPUNIT_ASSERT(Fs::join("", "", "", "") == "");
  //
  CPPUNIT_ASSERT(Fs::join("a", "", "", "") == "a");
  CPPUNIT_ASSERT(Fs::join("", "b", "", "") == "b");
  CPPUNIT_ASSERT(Fs::join("", "", "c", "") == "c");
  CPPUNIT_ASSERT(Fs::join("", "", "", "d") == "d");
  CPPUNIT_ASSERT(Fs::join("a", "b", "c", "d") == "a/b/c/d");
  CPPUNIT_ASSERT(Fs::join("/", "b", "c", "d") == "/b/c/d");

  /// vector args
  std::vector<std::string> parts;
  CPPUNIT_ASSERT(Fs::join(parts) == "");
  //
  parts.clear();
  parts.push_back("a");
  CPPUNIT_ASSERT(Fs::join(parts) == "a");
  parts.push_back("b");
  CPPUNIT_ASSERT(Fs::join(parts) == "a/b");
  //
  parts.clear();
  parts.push_back(".");
  CPPUNIT_ASSERT(Fs::join(parts) == ".");
  parts.push_back("a");
  CPPUNIT_ASSERT(Fs::join(parts) == "./a");
  parts.push_back("b");
  CPPUNIT_ASSERT(Fs::join(parts) == "./a/b");
  //
  parts.clear();
  parts.push_back("/");
  CPPUNIT_ASSERT(Fs::join(parts) == "/");
  parts.push_back("a");
  CPPUNIT_ASSERT(Fs::join(parts) == "/a");
  parts.push_back("b");
  CPPUNIT_ASSERT(Fs::join(parts) == "/a/b");
}

void FsTest::test_split()
{
  std::string drive;
  std::vector<std::string> parts;

  Fs::splitPath("", drive, parts);
  CPPUNIT_ASSERT(drive == "");
  CPPUNIT_ASSERT(parts.size() == 0);

  Fs::splitPath("a:", drive, parts);
  CPPUNIT_ASSERT(drive == "a:");
  CPPUNIT_ASSERT(parts.size() == 0);

  Fs::splitPath("b:/", drive, parts);
  CPPUNIT_ASSERT(drive == "b:");
  CPPUNIT_ASSERT(parts.size() == 1);
  CPPUNIT_ASSERT(parts[0] == "/");

  Fs::splitPath("c:abc", drive, parts);
  CPPUNIT_ASSERT(drive == "c:");
  CPPUNIT_ASSERT(parts.size() == 1);
  CPPUNIT_ASSERT(parts[0] == "abc");

  Fs::splitPath("d:/abc", drive, parts);
  CPPUNIT_ASSERT(drive == "d:");
  CPPUNIT_ASSERT(parts.size() == 2);
  CPPUNIT_ASSERT(parts[0] == "/");
  CPPUNIT_ASSERT(parts[1] == "abc");

  Fs::splitPath("e:\\abc\\def", drive, parts);
  CPPUNIT_ASSERT(drive == "e:");
  CPPUNIT_ASSERT(parts.size() == 3);
  CPPUNIT_ASSERT(parts[0] == "/"); // normalized seperator
  CPPUNIT_ASSERT(parts[1] == "abc");
  CPPUNIT_ASSERT(parts[2] == "def");
}


void FsTest::test_base_dir_ext()
{
  // should be an error
  CPPUNIT_ASSERT(Fs::basename("") == "");
  CPPUNIT_ASSERT(Fs::basename("foo") == "foo");
  CPPUNIT_ASSERT(Fs::basename("foo.bar") == "foo.bar");
  CPPUNIT_ASSERT(Fs::basename("a/b/c/foo.bar") == "foo.bar");
  CPPUNIT_ASSERT(Fs::basename("/foo.bar") == "foo.bar");
  CPPUNIT_ASSERT(Fs::basename(".foo") == ".foo");

  // should be an error.
  //CPPUNIT_ASSERT(Fs::dirname("")=="");
  CPPUNIT_ASSERT(Fs::dirname("/") == "/");
  CPPUNIT_ASSERT(Fs::dirname(".") == ".");
  CPPUNIT_ASSERT(Fs::dirname("foo.bar") == ".");
  CPPUNIT_ASSERT(Fs::dirname("a/b/c") == "a/b");
  CPPUNIT_ASSERT(Fs::dirname("a/b/c/") == "a/b");
  CPPUNIT_ASSERT(Fs::dirname("/a/b/c/") == "/a/b");

  // Test cases preserved from util/CPPTest/UtilTest.cpp

  CPPUNIT_ASSERT(Fs::dirname(INPUT + "\\t.e_s#-$tc\\") == INPUT + "\\t.e_s#-$tc");
  CPPUNIT_ASSERT(Fs::dirname(INPUT + "\\t.e_s#-$tc$") == INPUT);
  CPPUNIT_ASSERT(Fs::dirname(".\\input$t.e_s#-$tc") == ".");
  CPPUNIT_ASSERT(Fs::dirname(INPUT + "/t.e_s#-$tc/") == INPUT + "\\t.e_s#-$tc");
  CPPUNIT_ASSERT(Fs::dirname(INPUT + "/t.e_s#-$tc%") == INPUT);
  CPPUNIT_ASSERT(Fs::dirname("./input@t.e_s#-$tc") == ".");
  CPPUNIT_ASSERT(Fs::dirname("\\\\home\\foo.txt") == "\\home");
  CPPUNIT_ASSERT(Fs::dirname("\\\\home\\bar\\foo.txt") == "\\home\\bar");
  ///@todo this is not really the desired answer. The right answer is probably "D:\\"
  CPPUNIT_ASSERT(Fs::dirname("D:\\foo.txt") == "D:");
  CPPUNIT_ASSERT(Fs::dirname("D:\\home\\foo.txt") == "D:\\home");
  CPPUNIT_ASSERT(Fs::dirname("D:\\home\\bar\\foo.txt") == "D:\\home\\bar");
  CPPUNIT_ASSERT(Fs::dirname(INPUT + "/t.e_s#-$tc/") == INPUT + "/t.e_s#-$tc");
  CPPUNIT_ASSERT(Fs::dirname(INPUT + "/t.e_s#-$tc%") == INPUT);
  CPPUNIT_ASSERT(Fs::dirname("./input@t.e_s#-$tc") == ".");
  CPPUNIT_ASSERT(Fs::dirname("/home/foo.txt") == "/home");
  CPPUNIT_ASSERT(Fs::dirname("/home/bar/foo.txt") == "/home/bar");

  // "extname" is the longest ext! Leftmost dot!
  // to get a partial, use extnameN.
  CPPUNIT_ASSERT(Fs::extname("/") == "");
  CPPUNIT_ASSERT(Fs::extname(".") == ".");
  CPPUNIT_ASSERT(Fs::extname("..") == "..");
  CPPUNIT_ASSERT(Fs::extname(".a") == ".a");
  CPPUNIT_ASSERT(Fs::extname("a..") == "..");
  CPPUNIT_ASSERT(Fs::extname("foo.bar") == ".bar");
  CPPUNIT_ASSERT(Fs::extname("a/b/c.d") == ".d");
  CPPUNIT_ASSERT(Fs::extname("a/b/c") == "");
  CPPUNIT_ASSERT(Fs::extname("/a/b.c/d") == "");
  CPPUNIT_ASSERT(Fs::extname("/a/b.c/d.e") == ".e");
  CPPUNIT_ASSERT(Fs::extname("c.d.e") == ".d.e");
  CPPUNIT_ASSERT(Fs::extname("/a/b/c.d.e") == ".d.e");

  // This is the N rightmost exts.
  // CPPUNIT_ASSERT(Fs::extnameN("/a/b/c.d.e",1)==".e");
  // CPPUNIT_ASSERT(Fs::extnameN("/a/b/c.d.e",2)==".d.e");
  // CPPUNIT_ASSERT(Fs::extnameN("/a/b/c.d.e",9)==".d.e");

  //
  CPPUNIT_ASSERT(Fs::noextname("") == "");
  CPPUNIT_ASSERT(Fs::noextname("/") == "/");
  CPPUNIT_ASSERT(Fs::noextname("/.") == "/");
  CPPUNIT_ASSERT(Fs::noextname("foo") == "foo");
  CPPUNIT_ASSERT(Fs::noextname("foo.bar") == "foo");
  CPPUNIT_ASSERT(Fs::noextname("a/b/c/foo.bar.baz") == "a/b/c/foo");
  CPPUNIT_ASSERT(Fs::noextname("a.b.c/foo.bar.baz") == "a.b.c/foo");

  // extnameN
  CPPUNIT_ASSERT(Fs::noextname1("") == "");
  CPPUNIT_ASSERT(Fs::noextname1("/") == "/");
  CPPUNIT_ASSERT(Fs::noextname1("/.") == "/");
  CPPUNIT_ASSERT(Fs::noextname1("foo") == "foo");
  CPPUNIT_ASSERT(Fs::noextname1("foo.bar") == "foo");
  CPPUNIT_ASSERT(Fs::noextname1("a/b/c/foo.bar.baz") == "a/b/c/foo.bar");
  CPPUNIT_ASSERT(Fs::noextname1("a.b.c/foo.bar.baz") == "a.b.c/foo.bar");

  //
  CPPUNIT_ASSERT(Fs::noextnameN("", 3) == "");
  CPPUNIT_ASSERT(Fs::noextnameN("a", 3) == "a");
  CPPUNIT_ASSERT(Fs::noextnameN("XXXXXXXXXX/aaa.bbb.ccc.ddd.eee", 1) == "XXXXXXXXXX/aaa.bbb.ccc.ddd");
  CPPUNIT_ASSERT(Fs::noextnameN("1.2/3.4/a.b.c.d.e", 1) == "1.2/3.4/a.b.c.d");
  CPPUNIT_ASSERT(Fs::noextnameN("a.b.c", 3) == "a");
  CPPUNIT_ASSERT(Fs::noextnameN("a.b.c.d", 1) == "a.b.c");
  CPPUNIT_ASSERT(Fs::noextnameN("a.b.c.d", 2) == "a.b");
  CPPUNIT_ASSERT(Fs::noextnameN("a.b.c.d", 3) == "a");
  CPPUNIT_ASSERT(Fs::noextnameN("a.b.c.d", 4) == "a");
  CPPUNIT_ASSERT(Fs::noextnameN("1.2.3/a", 3) == "1.2.3/a");
  CPPUNIT_ASSERT(Fs::noextnameN("1.2.3/a.b", 3) == "1.2.3/a");
  CPPUNIT_ASSERT(Fs::noextnameN("1.2.3/a.b.c", 1) == "1.2.3/a.b");
  CPPUNIT_ASSERT(Fs::noextnameN("1.2.3/a...", 1) == "1.2.3/a..");
}

void FsTest::test_directoryReadable()
{
#ifdef WIN32
  CPPUNIT_ASSERT(Fs::isReadableDir(INPUT + "\\"));
  AffxString s1 = INPUT + "\\";
  const char * path = s1.c_str();
  CPPUNIT_ASSERT(Fs::isReadableDir(path));
  std::string path_s = INPUT + "\\";
  CPPUNIT_ASSERT(Fs::isReadableDir(path_s));

  CPPUNIT_ASSERT(Fs::isReadableDir(INPUT + "/chp")); //true
  //CPPUNIT_ASSERT(!Fs::isReadableDir(INPUT + "/chp/")); //todo vliber / return false

  std::string path_s1_w = INPUT + "/calvin";
  CPPUNIT_ASSERT(Fs::isReadableDir(path_s1_w)); //true
  //std::string path_s1 = INPUT + "/calvin/";
  //  CPPUNIT_ASSERT(!Fs::isReadableDir(path_s1)); //todo vliber / return false

  AffxString s2_w = INPUT + "/cel";
  const char* path_c_w = s2_w.c_str();
  CPPUNIT_ASSERT(Fs::isReadableDir(path_c_w)); //true
  //AffxString s2 = INPUT + "/cel/";
  //const char* path_c = s2.c_str();
  //CPPUNIT_ASSERT(!Fs::isReadableDir(path_c)); //todo vliber / return false
  //negative
  //CPPUNIT_ASSERT(!Fs::isReadableDir(".\\input1\\"));//todo vliber / return false
#else
  CPPUNIT_ASSERT(Fs::isReadableDir(INPUT + "/chp/"));
  std::string path_s1 = INPUT + "/calvin/";
  CPPUNIT_ASSERT(Fs::isReadableDir(path_s1));
  AffxString s1 = INPUT + "/cel/";
  const char* path_c = s1.c_str();
  CPPUNIT_ASSERT(Fs::isReadableDir(path_c));
#endif
}

void FsTest::test_mkdir()
{
  // use a name which should be unique.
  // start with an "A" to be at the top of the dir.
  std::string tmp_dir_name = "./A-FsTest-123-456-cpptest";

  // make sure it doesnt exist before the test.
  Fs::rmdir(tmp_dir_name, false);
  CPPUNIT_ASSERT(Fs::dirExists(tmp_dir_name) == false);

  //
  CPPUNIT_ASSERT(Fs::mkdir(tmp_dir_name) == APT_OK);

  //
  CPPUNIT_ASSERT(Fs::isReadable(tmp_dir_name) == true);
  CPPUNIT_ASSERT(Fs::isReadableDir(tmp_dir_name) == true);

  // this shouldn't work - rmdir does not use splitPath to deal with \\'s in path, unlike mkdirPath.
  //
  //POSITIVE_TEST((Fs::mkdirPath(".\\output\\t.e_s#-$tc\\") == APT_OK));
  //POSITIVE_TEST((Fs::rmdir(std::string(".\\output\\t.e_s#-$tc\\")) == APT_OK));
  //
  POSITIVE_TEST((Fs::mkdirPath(OUTPUT + "/t.e_s#-$tcU/") == APT_OK));
  POSITIVE_TEST((Fs::rmdir(OUTPUT + "/t.e_s#-$tcU/") == APT_OK));

  //positive char
  AffxString s1 = INPUT + "/t.e_s#-$tc";
  const char * dirName = s1.c_str();
  CPPUNIT_ASSERT(Fs::mkdirPath(std::string(dirName)) == APT_OK);
  CPPUNIT_ASSERT(Fs::isReadableDir(dirName));
  CPPUNIT_ASSERT(Fs::isWriteableDir(dirName));
  Fs::rmdirPath(dirName, false);

  //negative
#ifndef _MSC_VER
  AffxString s2 = INPUT + "/input1";
  std::string dirName2 = s2.c_str();
  Fs::mkdirPath(dirName2);
  Fs::chmodBasic(dirName2, 0);
  CPPUNIT_ASSERT(!Fs::isWriteableDir(dirName2));
  AffxString s22 = INPUT + "/input1/t.e_s#-$t/";
  std::string dirName22 = s22.c_str();
  //Goal is to receive a message: FATAL ERROR: Error: Util::makeDir() - failed to make directory .\input1\t.e_s#-$t\ -OK
  NEGATIVE_TEST((Fs::mkdirPath(dirName22) == APT_OK), std::exception);
  Fs::chmodBasic(dirName2, 0777);
  Fs::rmdir(dirName2, false);
#endif

  //positive string
  std::string dirName_s = INPUT + "/t.e_s#-$ts";
  CPPUNIT_ASSERT(Fs::mkdirPath(dirName_s) == APT_OK);
  CPPUNIT_ASSERT(Fs::isReadableDir(dirName_s));
  CPPUNIT_ASSERT(Fs::isWriteableDir(dirName_s));
  Fs::rmdir(dirName_s, false);
  //positive literal
  CPPUNIT_ASSERT(Fs::mkdirPath(INPUT + "/t.e_s#-$tl/") == APT_OK);
  CPPUNIT_ASSERT(Fs::isReadableDir(INPUT + "/t.e_s#-$tl"));
  CPPUNIT_ASSERT(Fs::isWriteableDir(INPUT + "/t.e_s#-$tl"));
  AffxString s3 = INPUT + "/t.e_s#-$tl/";
  Fs::rmdir(s3, false);
  //negative
  std::string dirName1_s = "./input1/t.e_s#-$t/";
  //Goal is to receive a message: FATAL ERROR: Error: Util::makeDir() - failed to make directory .\input1\t.e_s#-$t\ -OK
  NEGATIVE_TEST((Fs::mkdir(dirName1_s) == APT_OK), std::exception);

}

void FsTest::test_ensurewriteabledir()
{
  // use a name which should be unique.
  // start with an "A" to be at the top of the dir.
  std::string tmp_dir_name = "a/b/c";

  // make sure it doesnt exist before the test.
  Fs::rmdirPath(tmp_dir_name, false);
  CPPUNIT_ASSERT(Fs::dirExists(tmp_dir_name) == false);

  POSITIVE_TEST(Fs::ensureWriteableDirPath(tmp_dir_name, true) == APT_OK);
  
  CPPUNIT_ASSERT(Fs::rmdirPath(tmp_dir_name, false) == APT_OK);
}

void FsTest::test_touch()
{
  std::string tmp_file_name = "./A-FsTest-123-456.test";
}

void FsTest::test_iswriteabledir()
{
#if 0
  Fs p1;
  p1.setDirName(".");
  CPPUNIT_ASSERT(p1.osIsWriteableDir() == true);

  //  p1.setDirName("C:\\My Documents");
  p1.setDirName("C:\\Documents and Settings\\hgorre\\My Documents");
  p1.dump();
  CPPUNIT_ASSERT(p1.osIsWriteableDir() == true);
#endif
}

void FsTest::test_listdir()
{
  std::vector<std::string> names;

  //
  CPPUNIT_ASSERT(Fs::listDir(".", names, false) == APT_OK);
  CPPUNIT_ASSERT(names.size() != 0);

  CPPUNIT_ASSERT(Fs::listDir("/", names, false) == APT_OK);
  CPPUNIT_ASSERT(names.size() != 0);

  // make sure names really is zeroed out.
  names.push_back("junk");
  CPPUNIT_ASSERT(Fs::listDir("This/Dir/Does/Not/Exist", names, false) != APT_OK);
  CPPUNIT_ASSERT(names.size() == 0);
}

void FsTest::test_samevolume()
{
  AptErr_t rv;
  CPPUNIT_ASSERT(Fs::isSameVolume(".", ".", rv, false) == true);
  CPPUNIT_ASSERT(Fs::getErrNum() == APT_OK);
  CPPUNIT_ASSERT(Fs::isSameVolume("/", "/", rv, false) == true);
  CPPUNIT_ASSERT(Fs::getErrNum() == APT_OK);
}

//////////

void FsTest::test_freediskspace()
{
  CPPUNIT_ASSERT(Fs::getFreeDiskSpace(".") > 0);
  CPPUNIT_ASSERT(Fs::getErrNum() == APT_OK);

  CPPUNIT_ASSERT(Fs::getFreeDiskSpace("/thispathdoesnotexist", false) == -1);
  CPPUNIT_ASSERT(Fs::getErrNum() != APT_OK);
}


void FsTest::test_convertToUncPath()
{
  CPPUNIT_ASSERT(Fs::convertCommandToUnc("") == "");
  CPPUNIT_ASSERT(Fs::convertCommandToUnc(" abc") == " abc");
  CPPUNIT_ASSERT(Fs::convertCommandToUnc("abc ") == "abc ");
  CPPUNIT_ASSERT(Fs::convertCommandToUnc("abc") == "abc");

#ifdef _MSC_VER
  CPPUNIT_ASSERT(Fs::convertToUncPath("\\\\?\\c:\\path") == "\\\\?\\c:\\path");
  // what should happen when the path does not exist?
#else
  CPPUNIT_ASSERT(Fs::convertToUncPath("\\\\?\\c:\\path") == "//?/c:/path");
  CPPUNIT_ASSERT(Fs::convertCommandToUnc("\\a\\b\\c -arg1 -arg2  -arg3") == "/a/b/c -arg1 -arg2  -arg3");
  CPPUNIT_ASSERT(Fs::convertCommandToUnc("/usr/bin/ls -al") == "/usr/bin/ls -al");
#endif

#ifdef _WIN32
  char input1[] = "../../../../../input";
  std::string compare1("..\\..\\..\\..\\..\\input");
  CPPUNIT_ASSERT_EQUAL(Fs::convertToUncPath(input1), compare1);
  CPPUNIT_ASSERT_EQUAL(Fs::convertToUncPath("../../../../../input"), compare1);
  std::string affy_only_network_long_path = \
	  "//ntfs60.ev.affymetrix.com/bioinformatics/apt/apt-regression-data/trunk/idata/windows-max-path-40-char-long-dir-name-1/windows-max-path-40-char-long-dir-name-2/windows-max-path-40-char-long-dir-name-3/windows-max-path-40-char-long-dir-name-4/windows-max-path-40-char-long-dir-name-5/windows-max-path-40-char-long-dir-name-6/windows-max-path-40-char-long-dir-name-7/empty-file.txt";

  affy_only_network_long_path = Fs::Unc(affy_only_network_long_path);

  if ( Fs::dirExists( "//ntfs60.ev.affymetrix.com/bioinformatics") ) {
	  Verbose::out(1,"");
	  Verbose::out(1, std::string("Windows Max Path Network Test: ") + affy_only_network_long_path );
	  CPPUNIT_ASSERT(Fs::fileExists(affy_only_network_long_path ));
  }
  else{ 
	  Verbose::out(1,"");
	  Verbose::out(1, std::string("Windows Max Path Network Test skipped, network path not found: //ntfs60.ev.affymetrix.com/bioinformatics"));
  } 

  
#else
  char input2[] = "..\\..\\..\\..\\..\\input";
  std::string compare2("../../../../../input");
  CPPUNIT_ASSERT_EQUAL(Fs::convertToUncPath(input2), compare2);
  CPPUNIT_ASSERT_EQUAL(Fs::convertToUncPath("..\\..\\..\\..\\..\\input"), compare2);
#endif
  // removing this test - no such checks in code - AK
  //
  //negative
  //Goal is to receive a message: FATAL ERROR: Can't convert C:/apt/affy/sdk/util/CPPTest/input as it contains a ':' character -OK
  //char input3[] = "C:/apt/affy/sdk/util/CPPTest/input";
  //NEGATIVE_TEST(Fs::convertToUncPath(input3), std::exception);

}

void FsTest::test_permissions()
{

  std::cerr << std::endl << "FsTest::test_permissions START" << std::endl;
  std::string permFile = "test_permissions.txt";
  CPPUNIT_ASSERT(Fs::touch(permFile, false) == 0);
  CPPUNIT_ASSERT(Fs::isReadable(permFile));
#ifndef _MSC_VER
  // When both constants are given, they are joined with the bitwise OR operator ( |  ).
  //If write permission is not given, the file is read-only. Note that all files are always readable;
  //it is not possible to give write-only permission. Thus, the modes _S_IWRITE and _S_IREAD | _S_IWRITE are equivalent.
  CPPUNIT_ASSERT(Fs::chmodBasic(permFile, 0) == 0);
  CPPUNIT_ASSERT(!Fs::isReadable(permFile));
  CPPUNIT_ASSERT(Fs::chmodBasic(permFile, 0755) == 0);
  CPPUNIT_ASSERT(Fs::isReadable(permFile));
  CPPUNIT_ASSERT(Fs::isWriteable(permFile));
  CPPUNIT_ASSERT(Fs::chmodBasic(permFile, 0555) == 0);
  CPPUNIT_ASSERT(!Fs::isWriteable(permFile));
#endif
  Fs::rm(permFile);
  std::cerr << "FsTest::test_permissions END" << std::endl;

}

//test fileRemove() rsatin
void FsTest::test_fileRemove()
{
  Verbose::out(1, "Fs::test_fileRemove\n");
  if ( ! Fs::dirExists(OUTPUT) ) {
    Fs::mkdir(OUTPUT);
  }
  // setup for test case
  std::string tmpDir          = OUTPUT + "/output_util_FileIO/";               // temporary directory for testing
  std::string tmpDir_no_slash = OUTPUT + "/output_util_FileIO";                // FIXME: fix bug in chompLastIfSep to avoid this

  if ( ! Fs::dirExists(tmpDir_no_slash)) {
    Fs::mkdirPath(tmpDir_no_slash, false);
  }
  CPPUNIT_ASSERT(Fs::isReadableDir(tmpDir_no_slash));
  CPPUNIT_ASSERT(Fs::isWriteableDir(tmpDir_no_slash));
  std::string srcFile = INPUT + "/testFileOperations.txt";                    // copies of arbitrary text file for testing
  if(!Fs::fileExists(tmpDir + "file1.txt"))                                 // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, tmpDir + "file1.txt"));
  if(!Fs::fileExists(tmpDir + "file2.txt"))                                 // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, tmpDir + "file2.txt"));
  if(!Fs::fileExists(tmpDir + "file3.txt"))                                 // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, tmpDir + "file3.txt"));
  if(!Fs::fileExists(tmpDir + "file4.txt"))                                 // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, tmpDir + "file4.txt"));
  if(!Fs::fileExists(tmpDir + "file5.txt"))                                 // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, tmpDir + "file5.txt"));
  CPPUNIT_ASSERT(Fs::fileExists(tmpDir + "file1.txt"));                     // insure copies of files exist
  CPPUNIT_ASSERT(Fs::fileExists(tmpDir + "file2.txt"));
  CPPUNIT_ASSERT(Fs::fileExists(tmpDir + "file3.txt"));
  CPPUNIT_ASSERT(Fs::fileExists(tmpDir + "file4.txt"));
  CPPUNIT_ASSERT(Fs::fileExists(tmpDir + "file5.txt"));
  // test with std:string input argument with single filename to delete
  CPPUNIT_ASSERT(Fs::rm(tmpDir + "file1.txt", false) == APT_OK);             // delete existent file (positive case)
  CPPUNIT_ASSERT(!Fs::fileExists(tmpDir + "file1.txt"));                    // insure file deleted
  CPPUNIT_ASSERT(Fs::rm(tmpDir + "file0.txt", false) != APT_OK);             // delete non-existent file, boolean status
  POSITIVE_TEST(Fs::rmIfExists(tmpDir + "file0.txt", true));              // delete non-existent file, no exception
#ifdef WIN32
  std::string fileNameInUse = tmpDir + "file2.txt";                         // open for read to lock file blocking deletion
  FILE *f = fopen(fileNameInUse.c_str(), "r");                              // blocks deletion in Windows only
  CPPUNIT_ASSERT(f != NULL);
#else
  std::string fileNameInUse = tmpDir_no_slash;                              // non-empty directory is blocked from deletion
#endif
#ifdef WIN32  //todo rsatin - hangs on Linux and Darwin ---v
  CPPUNIT_ASSERT(Fs::rm(fileNameInUse, false) != APT_OK);                    // delete locked file (negative case)
  CPPUNIT_ASSERT(Fs::fileExists(fileNameInUse));                            // insure file still exists
  NEGATIVE_TEST(Fs::rm(fileNameInUse, true) == APT_OK, std::exception);     // delete locked file, throw exception (neg case)
  CPPUNIT_ASSERT(Fs::fileExists(fileNameInUse));                            // insure file still exists
#else
  printf("\nfileRemoveTest fileNameInUse: %s\n", fileNameInUse.c_str());
#endif        //todo rsatin - hangs on Linux and Darwin ---^
  // test with std::vector input argument with list of filenames to delete
  std::vector<std::string> filesToRemove;
  filesToRemove.insert(filesToRemove.end(), tmpDir + "file3.txt");
  filesToRemove.insert(filesToRemove.end(), tmpDir + "file4.txt");
  filesToRemove.insert(filesToRemove.end(), tmpDir + "file5.txt");
  for(int i = 0; i < filesToRemove.size();  i++) {
    POSITIVE_TEST(Fs::rm(filesToRemove[i]));
  } // delete existent file (positive case)
  CPPUNIT_ASSERT(!Fs::fileExists(tmpDir + "file3.txt"));                    // insure files deleted
  CPPUNIT_ASSERT(!Fs::fileExists(tmpDir + "file4.txt"));
  CPPUNIT_ASSERT(!Fs::fileExists(tmpDir + "file5.txt"));
  for(int i = 0; i < filesToRemove.size();  i++) {
    POSITIVE_TEST(Fs::rmIfExists(filesToRemove[i]) == APT_OK);
  } // delete non-existent files, no exception (neg case)

#ifdef WIN32  //todo rsatin - hangs on Linux ---v
  // test with std::vector input argument with locked file to delete
  filesToRemove.erase(filesToRemove.begin(), filesToRemove.end());
  filesToRemove.insert(filesToRemove.end(), fileNameInUse);
  for(int i = 0; i < filesToRemove.size();  i++) {
    NEGATIVE_TEST((Fs::rm(filesToRemove[i]) == APT_OK), std::exception);
  }   // delete locked file, throw exception (neg case)
  try {
    Verbose::setLevel(-1);
    for(int i = 0; i < filesToRemove.size();  i++) {
      Fs::rm(filesToRemove[i]);
    }
    Verbose::setLevel(3);
  } catch(Except e) {
    Verbose::setLevel(3);
    std::string msg = e.what();
    std::string fname = Fs::basename(fileNameInUse);
    CPPUNIT_ASSERT(msg.find(fname) != std::string::npos);                        // insure locked file listed in error message
  }
#endif        //todo rsatin - hangs on Linux ---^
#ifdef WIN32
  CPPUNIT_ASSERT(fclose(f) != EOF);                                         // close file to unlock
#endif
  CPPUNIT_ASSERT(Fs::rm(tmpDir + "file2.txt", false) == APT_OK);            // delete existent file (positive case)
  CPPUNIT_ASSERT(!Fs::fileExists(tmpDir + "file2.txt"));                    // insure file no longer exists
  // cleanup for test case
  Fs::rmdirPath(tmpDir, false);

}

//test fileRename() rsatin
void FsTest::test_fileRename()
{
  Verbose::out(1, "FsTest::test_fileRename");
  if ( ! Fs::dirExists(OUTPUT) ) {
    Fs::mkdir(OUTPUT);
  }
  // setup for test case
  std::string tmpDir = Fs::join(INPUT, "output_util_FileIO");        // temporary directory for testing
  //if( !Util::directoryWritable(tmpDir) )
  //Util::createDir(tmpDir);
  if ( !Fs::dirExists(tmpDir) ) {
    Fs::mkdirPath(tmpDir);
  }
  CPPUNIT_ASSERT(Fs::isReadableDir(tmpDir));
  CPPUNIT_ASSERT(Fs::isWriteableDir(tmpDir));
  std::string srcFile = Fs::join(INPUT, "testFileOperations.txt");                  // copies of arbitrary text file for testing
  if(!Fs::fileExists(Fs::join(tmpDir,"fileA.txt")))                                            // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, Fs::join(tmpDir,"fileA.txt")));
  if(!Fs::fileExists(Fs::join(tmpDir, "fileB.txt")))                                   // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, Fs::join(tmpDir, "fileB.txt")));
  if(!Fs::fileExists(Fs::join(tmpDir, "fileC.txt")))                                            // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, Fs::join(tmpDir, "fileC.txt")));
  if(!Fs::fileExists(Fs::join(tmpDir,"fileD.txt")))                                            // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::fileCopy(srcFile, Fs::join(tmpDir, "fileD.txt")));
  if(Fs::fileExists(Fs::join(tmpDir ,"renamed.txt")))                                           // insure renamed file not-existent initially
    CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir , "renamed.txt"), false));
  CPPUNIT_ASSERT(Fs::fileExists(Fs::join(tmpDir, "fileA.txt")));                             // insure copy of files exist
  CPPUNIT_ASSERT(Fs::fileExists(Fs::join(tmpDir, "fileB.txt")));                             // insure copy of files exist
  CPPUNIT_ASSERT(Fs::fileExists(Fs::join(tmpDir, "fileC.txt")));                             // insure copy of files exist
  CPPUNIT_ASSERT(!Fs::fileExists(Fs::join(tmpDir, "renamed.txt")));                             // insure renamed file not-existent
  // test with std:string input arguments
  CPPUNIT_ASSERT(Fs::fileRename(Fs::join(tmpDir, "fileA.txt"), Fs::join(tmpDir,"renamed.txt"), false)); // rename existent file (positive case)
  CPPUNIT_ASSERT(Fs::fileExists(Fs::join(tmpDir, "renamed.txt")));                              // insure renamed file exists

  /* BUG rsatin 01may09 ==> fileRename never updates "status", always claims success (and never throws exception) even when failed

    CPPUNIT_ASSERT( !Fs::fileRename(tmpDir+"fileX.txt", tmpDir+"fileY.txt", false) );    // rename non-existing file (neg case)
    CPPUNIT_ASSERT( !Fs::fileExists(tmpDir+"fileY.txt") );                               // insure file still non-existent
    NEGATIVE_TEST( Fs::fileRename(tmpDir+"fileX.txt", tmpDir+"fileY.txt", true), Except );  // neg case with exception result
    CPPUNIT_ASSERT( !Fs::fileExists(tmpDir+"fileY.txt") );                               // insure file still non-existent
    CPPUNIT_ASSERT( !Fs::fileRename(tmpDir+"fileB.txt", tmpDir+"fileC.txt", false) );    // rename overwrites existing file (neg case)
    CPPUNIT_ASSERT(  Fs::fileExists(tmpDir+"fileB.txt") );                               // insure original file still exists
    CPPUNIT_ASSERT(  Fs::fileExists(tmpDir+"fileC.txt") );                               // insure duplicate file still exists
    NEGATIVE_TEST( Fs::fileRename(tmpDir+"fileB.txt", tmpDir+"fileC.txt", true), Except );  // neg case with exception result
    CPPUNIT_ASSERT(  Fs::fileExists(tmpDir+"fileB.txt") );                               // insure original file still exists
    CPPUNIT_ASSERT(  Fs::fileExists(tmpDir+"fileC.txt") );                               // insure duplicate file still exists

    // test attempted rename of locked file, collect results, validate after file unlocked
    std::string fileNameInUse = tmpDir+"fileD.txt";                                        // filename of locked file
    FILE *f = fopen( fileNameInUse.c_str(), "r" );                                         // lock file
    CPPUNIT_ASSERT( f != NULL );
    bool isUnexpectedException = false;
    bool fileRename_status;
    try {
      fileRename_status = Fs::fileRename(fileNameInUse, tmpDir+"renamed2.txt", false); // rename locked file, request boolean status
    } catch( Except e ) {
      isUnexpectedException = true;
    }
    bool isExpectedException = false;
    std::string msgExpectedException;
    bool isExpectedException = false;
    try {
      Fs::fileRename(fileNameInUse, tmpDir+"renamed2.txt", true);                      // rename locked file, request exception if fails
    } catch( Except e ) {
      isExpectedException = true;
      msgExpectedException = e.what();
    }
    CPPUNIT_ASSERT( fclose(f) != EOF );                                                    // unlock file (do it now before assertions can fail)
    f = NULL;
    CPPUNIT_ASSERT( !isUnexpectedException );                                              // insure no unexpected exception
    CPPUNIT_ASSERT( !fileRename_status );                                                  // insure failure status for locked file
    CPPUNIT_ASSERT( isExpectedException );                                                 // insure expected exception when requested
    if( isExpectedException ) {
      CPPUNIT_ASSERT( msgExpectedException.find("fileD.txt") != string::npos );          // insure locked file listed in error message
      CPPUNIT_ASSERT( msgExpectedException.find("renamed2.txt") != string::npos );       // insure target file listed in error message
    }

      // test attempted rename to invalid target filename
    CPPUNIT_ASSERT(  Fs::fileRename(fileNameInUse, "!@#$%^&:?*.txt", false) );            // rename existent file (positive case)
    CPPUNIT_ASSERT(  Fs::fileExists(fileNameInUse) );                                     // insure original file not renamed
      NEGATIVE_TEST( Fs::fileRename(fileNameInUse, "!@#$%^&?*.txt", true), Except );        // neg case with exception result
    CPPUNIT_ASSERT(  Fs::fileExists(fileNameInUse) );                                     // insure original file not renamed
  */
  // cleanup for test case
  CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"renamed.txt"), false) == APT_OK);           // delete files (some should no longer exist)
  //CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"renamed2.txt"), false) == APT_OK );
  //CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"fileA.txt"), false) == APT_OK);
  CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"fileB.txt"), false) == APT_OK);
  CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"fileC.txt"), false) == APT_OK);
  CPPUNIT_ASSERT(Fs::rm(Fs::join(tmpDir,"fileD.txt"), false) == APT_OK);
  Fs::rmdir(tmpDir);
}
//test carefulClose() vliber
void FsTest::test_carefulClose()
{
  Verbose::out(1, "FsTest::test_carefulClose");
  if ( ! Fs::dirExists(OUTPUT) ) {
    Fs::mkdir(OUTPUT);
  }
  std::ofstream out;
  Fs::mustOpenToWrite(out, std::string(OUTPUT + "/GrcSmd.grc"));
  CPPUNIT_ASSERT(out.is_open() == true);
  CPPUNIT_ASSERT(out.good() == true);
  Fs::carefulClose(out);
  CPPUNIT_ASSERT(out.bad() == false);
}

//test fileCopy() vliber
void FsTest::test_fileCopy()
{
  Verbose::out(1, "FsTest::test_fileCopy");
  if ( ! Fs::dirExists(OUTPUT) ) {
    Fs::mkdir(OUTPUT);
  }
  //test unix files
  CPPUNIT_ASSERT(Fs::fileCopy(INPUT + "/1lqcdfpsi_unix.bcdf", OUTPUT + "/out_unix.bcdf") == true);
  AffxString s1 = INPUT + "/1lqcdfpsi_unix.bcdf";
  const char * path1 = s1.c_str();
  std::ifstream ifs1(path1, std::ios::binary);
  AffxString s2 = OUTPUT + "/out_unix.bcdf";
  const char * path2 = s2.c_str();
  std::ifstream ifs2(path2, std::ios::binary);
  std::string line1;
  std::string line2;
  while(!ifs2.eof() || !ifs1.eof()) {
    getline(ifs2, line2);
    getline(ifs1, line1);
    CPPUNIT_ASSERT(line1 == line2);
  }
  ifs1.close();
  ifs2.close();
  // test windows files
  CPPUNIT_ASSERT(Fs::fileCopy(INPUT + "/1lqcdfpsi_win.bcdf", OUTPUT + "/out_win.bcdf") == true);
  AffxString s3 = INPUT + "/1lqcdfpsi_win.bcdf";
  const char * path3 = s3.c_str();
  std::ifstream ifs3(path3, std::ios::binary);
  AffxString s4 = OUTPUT + "/out_win.bcdf";
  const char * path4 = s4.c_str();
  std::ifstream ifs4(path4, std::ios::binary);
  while(!ifs3.eof() || !ifs4.eof()) {
    getline(ifs3, line2);
    getline(ifs4, line1);
    CPPUNIT_ASSERT(line1 == line2);
  }
  ifs3.close();
  ifs4.close();

  //test char* and std::string
  AffxString s5 = INPUT + "/1lqcdfpsi_unix.bcdf";
  const char * pathIn = s5.c_str();
  AffxString s6 = OUTPUT + "/out_unix.bcdf";
  const char * pathOut = s6.c_str();
  CPPUNIT_ASSERT(Fs::fileCopy(std::string(pathIn), std::string(pathOut)) == true);
  std::string pathIn_s = INPUT + "/1lqcdfpsi_win.bcdf";
  std::string pathOut_s = OUTPUT + "/out_win.bcdf";
  CPPUNIT_ASSERT(Fs::fileCopy(pathIn_s, pathOut_s) == true);
}

//test mustOpenToWriteTest() vliber
void FsTest::test_mustOpenToWrite()
{
  Verbose::out(1, "FsTest::test_mustOpenToWrite");
  if ( ! Fs::dirExists(OUTPUT) ) {
    Fs::mkdir(OUTPUT);
  }
  //positive char*
  std::ofstream out;
  AffxString s1 = OUTPUT + "/aA1.2_3-4&5_c.txt";
  const char * path = s1.c_str();
  Fs::mustOpenToWrite(out, std::string(path));
  char a = 'A';
  out.put(a);
  CPPUNIT_ASSERT(out.is_open() == true);
  CPPUNIT_ASSERT(out.good() == true);
  CPPUNIT_ASSERT(out.fail() == false);
  out.close();
  //positive string
  std::string s = OUTPUT + "/aA1.2_3-4&5_s.txt";
  Fs::mustOpenToWrite(out, s);
  char c = 'C';
  out.put(c);
  CPPUNIT_ASSERT(out.is_open() == true);
  CPPUNIT_ASSERT(out.good() == true);
  CPPUNIT_ASSERT(out.fail() == false);
  out.close();
  //positive literal
  Fs::mustOpenToWrite(out, std::string(OUTPUT + "/aA1.2_3-4&5_l.txt"));
  char d = 'D';
  out.put(d);
  CPPUNIT_ASSERT(out.is_open() == true);
  CPPUNIT_ASSERT(out.good() == true);
  CPPUNIT_ASSERT(out.fail() == false);
  out.close();
  //negative
  std::ofstream out1;
  //goal is to receive a message: FATAL ERROR: Couldn't open file: ./inputMissingDirectory/test.txt to write.
  NEGATIVE_TEST(Fs::mustOpenToWrite(out1, std::string("./inputMissingDirectory/test.txt")), std::exception);
  CPPUNIT_ASSERT(out1.is_open() == false);
  CPPUNIT_ASSERT(out1.good() == false);
  CPPUNIT_ASSERT(out1.fail() == true);
  out1.close();
}

void FsTest::test_trailingSlash() {
  // special
  CPPUNIT_ASSERT(Fs::hasTrailingSlash("")==false);
  CPPUNIT_ASSERT(Fs::hasTrailingSlash("/")==true);
  CPPUNIT_ASSERT(Fs::hasTrailingSlash("a")==false);
  CPPUNIT_ASSERT(Fs::hasTrailingSlash("a/")==true);

  // special
  CPPUNIT_ASSERT(Fs::addTrailingSlash("")=="");
  CPPUNIT_ASSERT(Fs::addTrailingSlash("/")=="/");
  //
  CPPUNIT_ASSERT(Fs::addTrailingSlash("a")=="a/");
  CPPUNIT_ASSERT(Fs::addTrailingSlash("a/")=="a/");
  CPPUNIT_ASSERT(Fs::addTrailingSlash("a/b")=="a/b/");
  CPPUNIT_ASSERT(Fs::addTrailingSlash("a/b/")=="a/b/");

  // special
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("")=="");
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("/")=="/");
  //
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("a")=="a");
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("a/")=="a");
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("a/b")=="a/b");
  CPPUNIT_ASSERT(Fs::trimTrailingSlash("a/b/")=="a/b");
}

void FsTest::test_isbinary() {
  //
  CPPUNIT_ASSERT(Fs::isBinaryFile("input/hdf5/empty.hdf5")==true);
  //
  CPPUNIT_ASSERT(Fs::isBinaryFile("")==false);
  CPPUNIT_ASSERT(Fs::isBinaryFile("input/doesnotexist")==false);
  CPPUNIT_ASSERT(Fs::isBinaryFile("input/test.txt")==false);
}

void FsTest::test_iscalvin() {
  //
  CPPUNIT_ASSERT(Fs::isCalvinFile("input/calvin/copynumber/NA06985_GW6_C.CN5_small.CNCHP")==true);
  //
  CPPUNIT_ASSERT(Fs::isCalvinFile("")==false);
  CPPUNIT_ASSERT(Fs::isCalvinFile("input/doesnotexist")==false);
  CPPUNIT_ASSERT(Fs::isCalvinFile("input/test.txt")==false);
}

void FsTest::test_ishdf5()
{
  //
  CPPUNIT_ASSERT(Fs::isHdf5File("input/hdf5/empty.hdf5")==true);
  //
  CPPUNIT_ASSERT(Fs::isHdf5File("")==false);
  CPPUNIT_ASSERT(Fs::isHdf5File("input/hdf5/doesnotexist")==false);
  CPPUNIT_ASSERT(Fs::isHdf5File("input/test.txt")==false);
}
