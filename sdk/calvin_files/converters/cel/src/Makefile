#
# cvs:affy/sdk/file/util/chp_to_txt/Makefile ---
#
# $Id: Makefile,v 1.5 2009-11-07 05:51:08 mspald Exp $
#

#
sdk_root:=../../../..
include ${sdk_root}/Makefile.defs
include ../../../Makefile.inc
#
ctt_files:=$(filter-out CELFileConversionOptions.cpp, $(wildcard *.cpp ../../../converters/utils/src/CmdLine.cpp))
$(call sdk_define_exe,calvin-converter-cel,${ctt_files})
#
include ${sdk_makefile_post}
