#
# affy/sdk/chipstream/apt-geno-qc/regression/Makefile ---
#

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
t_pgms:=geno-qc-test
i_pgms:=geno-qc-quick-test
a_pgms:=${t_pgms} ${i_pgms} apt-geno-qc
_integration_run+=cp_pgms integration_1
_regression_run+=cp_pgms regression_1
_clean_run+=clean_regression
sdk_exe+=${t_pgms} ${i_pgms}
$(call sdk_define_exe,apt-rt-geno-qc,apt-rt-geno-qc.cpp)
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./geno-qc-quick-test

regression_1:
	./geno-qc-test

clean_integration:
	-rm -r geno-qc-snp5-quick-test.gqc  geno-qc-snp6-quick-test.gqc  geno-qc-sty-quick-test.gqc

clean_regression:
	-rm -r geno-qc-snp5-test.gqc  geno-qc-snp6-test.gqc  geno-qc-sty-test.gqc
