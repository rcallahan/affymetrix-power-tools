#
# cvs:affy/sdk/file/util/chp_to_txt/regression/Makefile ---
#
# $Id: Makefile,v 1.4 2009-09-24 00:10:11 alwon Exp $
#

t_pgms:=chp-to-txt-test
i_pgms:=chp-to-txt-quick-test
a_pgms:=${t_pgms} ${i_pgms} apt-chp-to-txt

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil m)
$(call sdk_define_exe,apt-rt-chp-to-txt,apt-rt-chp-to-txt.cpp)
#
sdk_exe+=${t_pgms} ${i_pgms}
_integration_run+=cp_pgms integration_1
_clean_run+=clean_integration
_regression_run+=cp_pgms regression_1
_clean_run+=clean_regression
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./chp-to-txt-quick-test

regression_1:
	./chp-to-txt-test

clean_integration:
	-rm -r test-out

clean_regression:
	-rm -r test-out
