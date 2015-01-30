#
# cvs:affy/sdk/chipstream/apt-probeset-genotype/regression/Makefile ---
#
# $Id: Makefile,v 1.5 2009-09-30 18:39:52 csugne Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
apg_exe:=apt-dmet-genotype
pgt_exe:=dmet-genotype-test
_integration_run+=integration_1 
_clean_run+=clean_integration
sdk_exe+=${pgt_exe}
$(call sdk_define_exe,apt-rt-dmet-genotype,apt-rt-dmet-genotype.cpp)
#
include ${sdk_makefile_post}

_integration_exe: ${sdk_output_bin}/${apg_exe}

cp_pgms: ${sdk_output_bin}/${pgt_exe} ${sdk_output_bin}/${apg_exe}
	${CP} ${sdk_output_bin}/${apg_exe} ${sdk_output_bin}/${pgt_exe} .

integration_1: cp_pgms
	./${pgt_exe}

clean_integration:
	-rm -r test-generated

