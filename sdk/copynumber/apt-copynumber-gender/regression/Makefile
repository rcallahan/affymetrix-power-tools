#
# # cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.9 2009-10-14 19:28:17 mgilch Exp $
#

t_pgms:=test-regression-gender  
a_pgms:=${t_pgms} apt-copynumber-gender

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
_clean_run+=clean_regression
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

regression_1:
	./test-regression-gender  ${valgrind_opts} ${database_opts}

clean_regression:
	-rm -r test-generated

