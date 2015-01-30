#
# cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.6 2009-09-24 00:10:11 alwon Exp $
#

t_pgms:=cel-convert-test
i_pgms:=cel-convert-quick-test
a_pgms:=${t_pgms} ${i_pgms} apt-cel-convert

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms} ${i_pgms}
_integration_run+=cp_pgms integration_1
_regression_run+=cp_pgms regression_1
_clean_run+=clean_integration clean_regression
$(call sdk_define_exe,apt-rt-cel-convert,apt-rt-cel-convert.cpp)
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./cel-convert-quick-test

regression_1:
	./cel-convert-test

clean_integration:
	-rm -r test-generated

clean_regression:
	-rm -r test-generated
