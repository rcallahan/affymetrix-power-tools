# This Makefile is for use within the APT package
#
# affy/sdk/calvin_files/Makefile.inc ---
#
# $Id: Makefile.inc,v 1.7 2009-09-25 17:48:52 mspald Exp $
#
# This makefile is included by the sub makefiles
# so we don not have to repeat all this in each calvin makefile.
#
$(call sdk_set_lib_name,affycalvin)
$(call sdk_set_link_libs,affysdk affycalvin affyfile affyutil xerces-c m)

# Not all the files have the path in the "#include"
# MYBRID CHANGE
#sdk_cpp_cflags+=-I${sdk_root}/portability

# Local Variables:
# mode: makefile
# End:
