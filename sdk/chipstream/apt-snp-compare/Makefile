#
# cvs:affy/sdk/chipstream/apt-snp-compare/Makefile ---
#
# $Id: Makefile,v 1.9 2009-09-25 17:48:57 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
#
$(call sdk_define_exe,apt-snp-compare,apt-snp-compare.cpp)
$(call sdk_define_install_exe,apt-snp-compare)
#
_check_run+=_test_snp_compare
#
include ${sdk_makefile_post}

_test_snp_compare: ${sdk_output_bin}/apt-snp-compare
	${sdk_output_bin}/apt-snp-compare test/calls.gold.txt test/calls.generated.txt 2> test/output.generated.txt
	diff -w test/output.generated.txt test/output.expected.txt

