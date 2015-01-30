#
# cvs:affy/sdk/chipstream/apt-probeset-genotype/regression/Makefile ---
#
# $Id: Makefile,v 1.22 2009-10-14 17:50:38 mgilch Exp $
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
apg_exe:=apt-probeset-genotype-stages 
# pgt_exe:=probeset-genotype-test
a_pgms:=${apg_exe} probeset-genotype-stages-quick-test
_regression_run+=regression_1 
_integration_run+=cp_pgms integration_1 
_valgrind_run+=cp_pgms valgrind_1 
_clean_run+=clean_regression
sdk_exe+=probeset-genotype-stages-quick-test ${pgt_exe} 
#$(call sdk_define_exe,apt-rt-probeset-genotype,apt-rt-probeset-genotype.cpp)
#
include ${sdk_makefile_post}

_regression_exe: ${sdk_output_bin}/${apg_exe}

cp_pgms: 
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1: cp_pgms
	./probeset-genotype-stages-quick-test

regression_1: cp_pgms

clean_regression:
	-rm -r test-generated

