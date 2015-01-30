#
# cvs:affy/sdk/translation/CPPTest/Makefile ---
#
# $Id: Makefile,v 1.10 2009-09-25 17:49:08 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore=$(wildcard *.cpp)
# Add the output directory to the list.
_check_run+=_mk_output_dir
_mk_output_dir: | output
output:
	mkdir ${@}
#
$(call sdk_define_check_exe,test-apt-dmet-translation ,$(wildcard *.cpp))
# say what we link against
#$(call sdk_set_link_libs, affysdk affycalvin affyfile affyutil pcrecpp pcre  m)

#
include ${sdk_makefile_post}

# Local Variables:
# mode: makefile
# End:
