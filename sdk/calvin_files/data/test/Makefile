#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.9 2009-04-16 16:32:00 awilli Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
include ../../Makefile.inc
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
# TODO: Handle missing data
ignore_files=\
	DataSetTest_LargeData.cpp \
	DataSetTest_RemapData.cpp \
	MemMapFileTest.cpp
#
$(call sdk_define_check_exe,test-calvin-data,$(filter-out $(ignore_files), $(wildcard *.cpp)))
#
include ${sdk_makefile_post}
