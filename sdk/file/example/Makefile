#
# cvs:affy/sdk/file/example/Makefile ---
#
# $Id: Makefile,v 1.17 2009-09-25 17:49:02 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
$(call sdk_set_link_libs,affyfile affyutil z m)
#
$(call sdk_define_exe,file-example,main.cpp ReadAffxFiles.cpp WriteAffxFiles.cpp)
#
include ${sdk_makefile_post}
