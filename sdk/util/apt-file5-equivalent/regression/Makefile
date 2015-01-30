#
# cvs:affy/sdk/util/apt-file5-equivalent/regression/Makefile ---
#
# $Id: Makefile,v 1.1 2009-10-27 21:37:22 rsatin Exp $
#

t_pgms:=file5-equivalent-test
a_pgms:=${t_pgms} apt-file5-equivalent

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
	./file5-equivalent-test

