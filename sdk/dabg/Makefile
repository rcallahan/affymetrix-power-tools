#
# cvs:affy/sdk/dabg/Makefile ---
#
# $Id: Makefile,v 1.8 2006-05-05 21:38:56 awilli Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
#
sdk_exe:=test-dabg
_check_run+=_test_dabg_run
#
include ${sdk_makefile_post}

# Hang this program on "_check_run"

_test_dabg_run: ${sdk_output_bin}/test-dabg
	${sdk_output_bin}/test-dabg
