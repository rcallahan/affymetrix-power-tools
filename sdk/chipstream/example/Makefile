#
# cvs:affy/sdk/chipstream/example/Makefile ---
#
# $Id: Makefile,v 1.23 2009-04-26 20:14:42 harley Exp $
#

# before include
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore+=file-stats.cpp
#$(call sdk_define_exe,file-stats,file-stats.cpp)
$(call sdk_define_exe,simpleRmaExample,simpleRmaExample.cpp)
$(call sdk_define_exe,multiChannelExample,multiChannelExample.cpp)
#sdk_cpp_ignore+=simpleRmaExample.cpp
#
include ${sdk_makefile_post}
