#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.22 2009-04-22 21:13:54 awilli Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affysdk)
$(call sdk_set_link_libs,affysdk affyutil m)
#
$(call sdk_define_exe,test-affy-random-sample,affy_random_sample_test.cpp)
#
sdk_subdirs:=CPPTest
include ${sdk_makefile_post}
