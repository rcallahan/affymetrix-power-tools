#
# cvs:affy/sdk/util/matrix-diff/regression/Makefile ---
#
# $Id: Makefile,v 1.7 2009-09-24 00:10:11 alwon Exp $
#

t_pgms:=md5-check-test
a_pgms:=${t_pgms} apt-md5-check

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
#$(call sdk_set_lib_name,affyutil)
#$(call sdk_set_link_libs,affyutil m)
$(call sdk_define_exe,md5-check-test,md5-check-test.cpp)
#
sdk_exe+=${t_pgms}
_integration_run+=cp_pgms integration_1
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./md5-check-test

