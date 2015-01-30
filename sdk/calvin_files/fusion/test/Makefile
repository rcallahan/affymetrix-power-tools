#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.12 2009-04-16 16:32:04 awilli Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
include ../../Makefile.inc
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
#$(call sdk_define_check_exe,test-calvin-fusion,$(filter-out $(ignore_files), $(wildcard *.cpp)))
$(call sdk_define_check_exe,test-calvin-fusion,$(wildcard *.cpp))
#
include ${sdk_makefile_post}
