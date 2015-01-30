#
# cvs:affy/sdk/chipstream/apt-cel-transformer/regression/Makefile ---
#
# $Id: Makefile,v 1.27 2009-09-24 00:10:11 alwon Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
act_exe:=apt-cel-transformer
ctt_exe:=cel-transformer-test
cti_exe:=cel-transformer-quick-test
_integration_run+=integration_1
_regression_run+=regression_1
_clean_run+=clean_integration clean_regression
$(call sdk_define_exe,apt-rt-cel-transformer,apt-rt-cel-transformer.cpp)
sdk_exe+=${ctt_exe} ${cti_exe}
#
include ${sdk_makefile_post}

integration_1: ${sdk_output_bin}/${apg_exe} cp_pgms
	./${cti_exe}

regression_1: ${sdk_output_bin}/${apg_exe} cp_pgms
	./${ctt_exe}

cp_pgms:
	${CP} ${sdk_output_bin}/${act_exe} ${sdk_output_bin}/${ctt_exe} ${sdk_output_bin}/${cti_exe} .

clean_integration:
	-rm -r qc-test-generated

clean_regression:
	-rm -r qc-test-generated
