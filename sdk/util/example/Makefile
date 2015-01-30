#
# cvs:affy/sdk/util/example/Makefile ---
#
# $Id: Makefile,v 1.6 2007-08-31 22:58:17 awilli Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
#$(call sdk_set_lib_name,affyutil)
#$(call sdk_set_link_libs,affyutil m)

# Append this to a variable which is used when running ${CXX}
$(call sdk_define_exe,cc-chp-compare,cc-chp-compare.cpp)
$(call sdk_define_exe,xda-chp-compare,xda-chp-compare.cpp)
$(call sdk_define_exe,cel-compare,cel-compare.cpp)
$(call sdk_define_exe,err-example,ErrExample.cpp)
$(call sdk_define_exe,mem-info,MemInfo.cpp)
$(call sdk_define_exe,socket-listener-demo,socket-listener-demo.cpp)
$(call sdk_define_exe,engine-runner,engine-runner.cpp)

include ${sdk_makefile_post}
