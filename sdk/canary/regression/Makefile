#
# cvs:affy/sdk/chipstream/apt-probeset-genotype/regression/Makefile ---
#
# $Id: Makefile,v 1.6 2009-10-20 17:34:07 csugne Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
apg_exe:=apt-canary
pgt_exe:=canary-test
pgi_exe:=canary-quick-test
_integration_run+=integration_1 
_regression_run+=regression_1 
_clean_run+=clean_integration clean_regression

sdk_exe+=${pgt_exe} ${pgi_exe}
$(call sdk_define_exe,apt-rt-canary,apt-rt-canary.cpp)
#
include ${sdk_makefile_post}

_integration_exe: ${sdk_output_bin}/${apg_exe}

cp_pgms: ${sdk_output_bin}/${pgt_exe} ${sdk_output_bin}/${apg_exe} ${sdk_output_bin}/${pgi_exe}
	${CP} ${sdk_output_bin}/${apg_exe} ${sdk_output_bin}/${pgt_exe} ${sdk_output_bin}/${pgi_exe} .

integration_1: cp_pgms
	./${pgi_exe}

regression_1: cp_pgms
	./${pgt_exe}

clean_integration:
	-rm -r test-generated

clean_regression:
	-rm -r test-generated

