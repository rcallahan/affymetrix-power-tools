#
# cvs:affy/sdk/dabg/Makefile ---
#
# $Id: Makefile,v 1.14 2009-09-25 17:48:56 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
#
$(call sdk_define_exe,apt-cel-transformer,$(wildcard *.cpp) ../../calvin_files/converters/cel/src/CELFileVersion.cpp)
$(call sdk_define_install_exe,apt-cel-transformer)
#
include ${sdk_makefile_post}
