#
# cvs:affy/sdk/dabg/Makefile ---
#
# $Id: Makefile,v 1.23 2009-09-25 17:48:57 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
sdk_cpp:=$(wildcard *.cpp)
#
#
$(call sdk_define_exe,apt-chp-to-table,apt-chp-to-table.cpp)
$(call sdk_define_install_exe,apt-chp-to-table)
$(call sdk_define_exe,apt-table-to-impute,apt-table-to-impute.cpp)
$(call sdk_define_install_exe,apt-table-to-impute)
#
include ${sdk_makefile_post}
