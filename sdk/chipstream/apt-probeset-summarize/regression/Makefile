#
# cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.45 2009-10-14 23:31:45 mgilch Exp $
#

t_pgms:=probeset-summarize-quick-test probeset-summarize-quick-test.valgrind probeset-summarize-test summary-vis-test 
a_pgms:=${t_pgms} apt-probeset-summarize apt-summary-vis

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
_integration_run+=cp_pgms integration_1
_valgrind_run+=cp_pgms valgrind_1
_clean_run+=clean_regression
$(call sdk_define_exe,apt-rt-probeset-summarize,apt-rt-probeset-summarize.cpp)
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

valgrind_1:
	./probeset-summarize-quick-test.valgrind

integration_1:
	./probeset-summarize-quick-test ${valgrind_opts} ${database_opts}

regression_1:
	./probeset-summarize-test ${valgrind_opts} ${database_opts}
	./summary-vis-test

clean_regression:
	-rm -r summary-vis-test-generated test-generated

_run_gdb:
	-rm -rf test-generated/qt-doHGT
	-mkdir test-generated
	gdb --args ${sdk_output_bin}/apt-probeset-summarize -a plier-gcbg-sketch -a rma-sketch -a quant-norm,pm-only,med-polish,pca-select -a quant-norm,pm-only,med-polish,spect-select -b ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.bgp -c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf -p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf --qc-probesets ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.qcc -m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.mps -o test-generated/qt-doHGT --use-disk=false --cel-files  ../../../regression-data/data/idata/p-sum/HuGene-1_0-st-v1/cel-files.txt --feat-effects --feat-details --write-sketch --cc-chp-output

