#
# affy/sdk/bboard/bboard/Makefile ---
#
# $Id: Makefile,v 1.4 2009-10-29 22:28:58 harley Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
include ${sdk_root}/bboard/Makefile.inc
#
$(call sdk_define_exe,test-bboard,test-bboard.cpp)
#
include ${sdk_makefile_post}


_test: _test_bboard

test_bboard_exe:=${sdk_output_bin}/test-bboard

_test_bboard: ${test_bboard_exe}
	${test_bboard_exe}

_test_bboard_gdb: ${test_bboard_exe}
	${gdb} ${test_bboard_exe}

_test_bboard_valgrind: ${test_bboard_exe}
	${valgrind} ${test_bboard_exe}
