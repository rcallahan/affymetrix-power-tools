#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.9 2009-09-25 17:48:53 mspald Exp $
#

#
sdk_root:=../../../..
include ${sdk_root}/Makefile.defs
include ../../../Makefile.inc
#
# Not all the files have
#
sdk_cpp_ingore=$(wildcard *.cpp)

testfiles=$(wildcard *.cpp) \
 ../src/CHPConversionUtilities.cpp \
 ../src/CHPFileConvertToCalvin.cpp \
 ../src/CHPFileConvertToXDA.cpp \
 ../src/CHPFileConverter.cpp \
 ../src/CHPFileConverterErrorCode.cpp \
 ../src/CHPFileVersion.cpp \
 ../src/CmdLine.cpp

#
$(call sdk_define_check_exe,test-calvin-convert-chp,$(testfiles))
#
include ${sdk_makefile_post}
