#
# cvs:affy/sdk/util/apt-calvin-equivalent/Makefile ---
#
# $Id: Makefile,v 1.1 2009-10-27 21:42:52 rsatin Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=-I../../external/xerces/src
#
$(call sdk_set_lib_name,affyutil)
#$(call sdk_set_link_libs,affyutil m xerces)
#
$(call sdk_define_exe,apt-calvin-equivalent,apt-calvin-equivalent.cpp)
$(call sdk_define_install_exe,apt-calvin-equivalent)

include ${sdk_makefile_post}
