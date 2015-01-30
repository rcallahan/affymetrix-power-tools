#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.7 2009-04-16 16:32:11 awilli Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
include ../../Makefile.inc
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
ignore_files=\
	MMDataSetDataWriterTest.cpp \
	MMDATFileWriterTest.cpp
$(call sdk_define_check_exe,test-calvin-writers,$(filter-out $(ignore_files), $(wildcard *.cpp)))
#
include ${sdk_makefile_post}

