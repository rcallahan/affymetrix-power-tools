#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.1 2009-04-16 14:40:18 awilli Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=-I${sdk_root}/mas5-stat/src -I${sdk_root}/mas5-stat/workflow
sdk_cpp_ignore=$(wildcard *.cpp)
#
$(call sdk_define_check_exe,test-mas5-stat,$(wildcard *.cpp))
#
include ${sdk_makefile_post}
