#
# cvs:affy/sdk/util/apt-file5-equivalent/Makefile ---
#
# $Id: Makefile,v 1.1 2009-10-27 21:39:09 rsatin Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=-I../../external/xerces/src
#
$(call sdk_set_lib_name,affyutil)
# $(call sdk_set_link_libs,affyutil asqlite3 m)
#
$(call sdk_define_exe,apt-file5-equivalent,apt-file5-equivalent.cpp)
$(call sdk_define_install_exe,apt-file5-equivalent)

include ${sdk_makefile_post}
