#
# cvs:affy/sdk/file/cdfexport/regression/Makefile ---
#
# $Id: Makefile,v 1.9 2009-09-24 00:10:12 alwon Exp $
#

t_pgms:=cdf-export-quick-test
a_pgms:=${t_pgms} apt-cdf-export

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil m)
$(call sdk_define_exe,apt-rt-cdf-export,apt-rt-cdf-export.cpp)
#
sdk_exe+=${t_pgms}
_integration_run+=cp_pgms integration_1
_clean_run+=clean_integration
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./cdf-export-quick-test

clean_integration:
	-rm -r test-generated
