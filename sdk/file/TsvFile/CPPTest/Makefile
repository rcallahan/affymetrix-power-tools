#
# cvs:affy/sdk/util/Makefile2 ---
#
# $Id: Makefile,v 1.7 2009-09-25 17:49:02 mspald Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil m)
#
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
$(call sdk_define_check_exe,test-tsv-file,$(wildcard *.cpp) ../TsvJoin.cpp ../DumpPgf.cpp)
#
include ${sdk_makefile_post}
#
