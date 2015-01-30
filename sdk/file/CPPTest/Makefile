#
# cvs:affy/sdk/util/Makefile2 ---
#
# $Id: Makefile,v 1.20 2009-09-25 17:49:01 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
# FIXME: Find data for excluded tests
ignore_files=\
	BPMAPFileDataTest.cpp \
	BPMAPFileWriterTest.cpp \
	CHPFileDataTest.cpp \
	DttArrayFileReaderTest.cpp \
	GRDFileDataTest.cpp \
	MDLFileDataTest.cpp \
	TSVFileDataTest.cpp
$(call sdk_define_check_exe,test-file,$(filter-out $(ignore_files), $(wildcard *.cpp)))
#

#
include ${sdk_makefile_post}
