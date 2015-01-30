#
# cvs:affy/sdk/util/apt-calvin-equivalent/regression/Makefile ---
#
# $Id: Makefile,v 1.1 2009-10-27 21:42:10 rsatin Exp $
#

t_pgms:=calvin-equivalent-test
a_pgms:=${t_pgms} apt-calvin-equivalent

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyutil)
$(call sdk_set_link_libs,affyutil xerces-c asqlite3 m)
#
sdk_exe+=${t_pgms}
_integration_run+=cp_pgms integration_1
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

integration_1:
	./calvin-equivalent-test

