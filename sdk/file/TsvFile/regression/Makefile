#
# cvs:affy/sdk/chipstream/apt-probeset-summarize/regression/Makefile ---
#
# $Id: Makefile,v 1.10 2009-09-30 18:39:52 csugne Exp $
#

t_pgms:=dump-pgf-test tsv-join-test
a_pgms:=${t_pgms} apt-dump-pgf apt-tsv-join

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil m)
$(call sdk_define_exe,apt-rt-tsv-join,apt-rt-tsv-join.cpp)
#
sdk_exe+=${t_pgms}
_integration_run+=cp_pgms integration_1
_clean_run+=clean_integration
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./tsv-join-test
	./dump-pgf-test

clean_integration:
	-rm -r tsv-join-test-generated dump-pgf-test-generated
