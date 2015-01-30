#
# cvs:affy/sdk/file/apt-cdf-export/Makefile ---
#
# $Id: Makefile,v 1.8 2009-09-25 17:49:02 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil z m)
#
#
$(call sdk_define_exe,apt-cdf-export,apt-cdf-export.cpp)
$(call sdk_define_install_exe,apt-cdf-export)
#
include ${sdk_makefile_post}



#
apt_cdf_export:=${sdk_output_bin}/apt-cdf-export

_test: _test_1 _test_2

_test_1:
	${apt_cdf_export} --help
	${apt_cdf_export} --version

test_2_cdf:=../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf
_test_2:
	${apt_cdf_export} --file-md5 ${test_2_cdf}
	md5sum ${test_2_cdf}
