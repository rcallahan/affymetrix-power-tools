#
# # cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.1 2009-10-01 00:41:09 wshort Exp $
#

t_pgms:=test-regression-copynumber-wave #test-regression-copynumber-wave-quick
a_pgms:=${t_pgms} apt-copynumber-wave

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
#_integration_run+=cp_pgms integration_1
_clean_run+=clean_regression
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

#integration_1:
#	./test-regression-copynumber-wave-quick

regression_1:
	./test-regression-copynumber-wave

clean_regression:
	-rm -r test-generated

