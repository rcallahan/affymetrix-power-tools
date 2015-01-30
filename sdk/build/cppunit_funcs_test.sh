#!/bin/bash
#
# sdk/build/cppunit_funcs_test ---
#
# $Id$
#

# a example/test of cppunit_funcs.sh

# should figure out the path automaticly.
sdk_root=..
source ${sdk_root}/build/cppunit_funcs.sh

# sets what will be the output file name and inits the library.
cppunit_set_outputfile "cppunit_check.xml"

# Sadly, the XSL file for converting CPPUnit data to JUnitTestResults data
# requires a "::" in the name.  So name all your tests with "BigName::LittleName".

cppunit_test_start "Foo::test_pass_1"
cppunit_test_passed

cppunit_test_start "Foo::test_pass_2"
cppunit_test_passed

cppunit_test_start "Foo::test_fail_3"
cppunit_test_failed ERR_NO_FISH "a failure of porpoise"

cppunit_test_start "Foo::test_pass_4"
cppunit_test_passed

cppunit_test_start "test_fail_5"
cppunit_test_failed ERR_NO_TUNE "You cant tune a fish"

# Calling this function gathers up all the test data, creates the output file
# and deletes the temp files. If you dont call it, you dont get any output.
cppunit_finish

# Local Variables:
# mode: ksh
# End:
