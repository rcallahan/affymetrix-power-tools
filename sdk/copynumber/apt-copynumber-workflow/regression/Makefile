#
# # cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.11 2009-10-15 16:37:25 mgilch Exp $
#

t_pgms:=test-regression-copynumber-workflow test-regression-copynumber-workflow-quick test-regression-copynumber-workflow-quick.valgrind

a_pgms:=${t_pgms} apt-copynumber-workflow

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
_integration_run+=cp_pgms integration_1
_valgrind_run+=cp_pgms valgrind_1
_clean_run+=clean_regression
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

valgrind_1:
	./test-regression-copynumber-workflow-quick.valgrind

integration_1:
	./test-regression-copynumber-workflow-quick  ${valgrind_opts} ${database_opts}
regression_1:
	./test-regression-copynumber-workflow  ${valgrind_opts} ${database_opts}

clean_regression:
	-rm -r test-generated


