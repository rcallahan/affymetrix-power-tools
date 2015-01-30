#
# cvs:affy/sdk/file/util/chp_to_txt/Makefile ---
#
# $Id: Makefile,v 1.3 2009-09-25 17:49:09 mspald Exp $
#

ctt_files:=$(wildcard *.cpp)
#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#
$(call sdk_define_exe,apt-chp-to-txt,${ctt_files})
$(call sdk_define_install_exe,apt-chp-to-txt)
#
include ${sdk_makefile_post}

