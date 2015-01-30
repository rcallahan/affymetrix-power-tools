#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.38 2009-04-16 16:32:13 awilli Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
$(call sdk_set_link_libs,affyfile affyutil z m)

sdk_cpp_ignore:=\
	test-fileio.cpp
#
include ${sdk_makefile_post}
