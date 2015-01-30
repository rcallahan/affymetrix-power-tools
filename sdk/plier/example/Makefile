#
# cvs:affy/sdk/plier/Makefile ---
#
# $Id: Makefile,v 1.12 2009-09-25 17:49:06 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_define_exe,plier-example,test.cpp ../plier_impl.cxx)
_check_run+=_plier_example_run
#
include ${sdk_makefile_post}

_plier_example_run: ${sdk_output_bin}/plier-example
	${sdk_output_bin}/plier-example inputs.txt > outputs.txt
	diff -w inputs.reference outputs.txt

