#
# cvs:affy/sdk/util/Makefile2 ---
#
# $Id: Makefile,v 1.6 2007-08-11 00:06:52 awilli Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=-I..
#
sdk_cpp_ingore=$(wildcard *.cpp)
#
$(call sdk_define_check_exe,test-portability,$(wildcard *.cpp))
#
include ${sdk_makefile_post}

# Local Variables:
# mode: makefile
# End:
