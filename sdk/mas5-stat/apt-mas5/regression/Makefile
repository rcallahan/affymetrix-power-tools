#
# cvs:affy/sdk/chipstream/apt-mas5/regression/Makefile ---
#
# $Id: Makefile,v 1.27 2009-09-24 00:10:11 alwon Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
act_exe:=apt-mas5
ctt_exe:=mas5-test
_regression_run+=regression_1
_clean_run+=clean_regression
sdk_exe+=${ctt_exe}
#
include ${sdk_makefile_post}

regression_1: ${sdk_output_bin}/${apg_exe} cp_pgms
	./${ctt_exe}

cp_pgms:
	${CP} ${sdk_output_bin}/${act_exe} ${sdk_output_bin}/${ctt_exe} .

clean_regression:
	-rm -r qc-test-generated
