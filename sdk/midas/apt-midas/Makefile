#
# cvs:affy/sdk/midas/apt/Makefile ---
#
# $Id: Makefile,v 1.13 2009-09-25 17:49:05 mspald Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
#
# ignore these windows files
sdk_cpp_ignore:=MidasProcessingBox.cpp apt-midas-win.cpp
#
$(call sdk_define_exe,apt-midas,apt-midas.cpp)
$(call sdk_define_install_exe,apt-midas)
#
include ${sdk_makefile_post}

