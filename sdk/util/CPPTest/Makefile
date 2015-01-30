#
# cvs:affy/sdk/util/CPPTest/Makefile ---
#
# $Id: Makefile,v 1.15 2009-09-25 17:49:08 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
# we used to do this, but now "util" sucks in alot more.
# util should be split into parts.
# $(call sdk_set_link_libs,affyutil m)

# dont process the other files
sdk_cpp_ignore:=$(wildcard *.cpp)
$(call sdk_define_check_exe,test-util,$(wildcard *.cpp))
# if you just want to do one test....
#$(call sdk_define_check_exe,test-util-fspath,FsTest.cpp)
$(call sdk_define_check_exe,test-util-fspath,FsTest.cpp FsPathTest.cpp)
#
include ${sdk_makefile_post}

#
_valgrind:
	valgrind -v --tool=memcheck --leak-check=full --show-reachable=no \
	  ${sdk_output_bin}/test-util

_test_fspath: ${sdk_output_bin}/test-util-fspath
	${sdk_output_bin}/test-util-fspath

# Local Variables:
# mode: makefile
# End:
