#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.4 2009-09-25 17:48:54 mspald Exp $
#

#
sdk_root:=../../../..
include ${sdk_root}/Makefile.defs
include ../../../Makefile.inc
#
#
testfiles=$(wildcard *.cpp)
sdk_cpp_ignore=$(wildcard *.cpp)
#
$(call sdk_define_check_exe,test-calvin-converter-util,$(testfiles))
#
include ${sdk_makefile_post}
