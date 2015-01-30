#
# cvs:affy/sdk/canary/Makefile ---
#
# $Id: Makefile,v 1.9 2009-04-22 21:13:53 awilli Exp $
#

sdk_root:=..
include ${sdk_root}/Makefile.defs

$(call sdk_define_exe,apt-canary,apt-canary.cpp)
$(call sdk_define_install_exe,apt-canary)

include ${sdk_makefile_post}

#
apt_canary:=${sdk_output_bin}/apt-canary

#
_run-apt-canary: ${apt_canary}
	${apt_canary}

_test_help: ${apt_canary}
	${apt_canary} --help

