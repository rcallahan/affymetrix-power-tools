#
# # cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.9 2009-10-14 19:28:17 mgilch Exp $
#

t_pgms:=test-regression-copynumber-cyto test-regression-copynumber-cyto-quick test-regression-copynumber-cyto-quick.valgrind 
a_pgms:=${t_pgms} apt-copynumber-cyto

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
_integration_run+=cp_pgms integration_1
_valgrind_run+=cp_pgms valgrind_1
_clean_run+=clean_regression
$(call sdk_define_exe,apt-rt-copynumber-cyto,apt-rt-copynumber-cyto.cpp)
$(call sdk_define_exe,createTsv,createTsv.cpp)
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

valgrind_1:
	./test-regression-copynumber-cyto-quick.valgrind

integration_1:
	./test-regression-copynumber-cyto-quick   ${valgrind_opts} ${database_opts}

regression_1:
	./test-regression-copynumber-cyto  ${valgrind_opts} ${database_opts}

clean_regression:
	-rm -r test-generated

