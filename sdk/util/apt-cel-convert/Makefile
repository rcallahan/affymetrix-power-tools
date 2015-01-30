#
# cvs:affy/sdk/file/util/chp_to_txt/Makefile ---
#
# $Id: Makefile,v 1.4 2009-09-25 17:49:09 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
conv_files:=../../calvin_files/converters/cel/src/CELConversionUtilities.cpp ../../calvin_files/converters/cel/src/CELFileConvertToASCII.cpp ../../calvin_files/converters/cel/src/CELFileConvertToCalvin.cpp ../../calvin_files/converters/cel/src/CELFileConvertToXDA.cpp ../../calvin_files/converters/cel/src/CELFileConverter.cpp ../../calvin_files/converters/cel/src/CELFileConverterErrorCode.cpp ../../calvin_files/converters/cel/src/CELFileVersion.cpp

ctt_files:=$(wildcard *.cpp ${conv_files})


$(call sdk_define_exe,apt-cel-convert,${ctt_files})
$(call sdk_define_install_exe,apt-cel-convert)
#
include ${sdk_makefile_post}
