#
# affy/sdk/chipstream/apt-cdf-to-spf/Makefile ---
#
# $Id: Makefile,v 1.3 2009-05-01 00:43:15 harley Exp $
#

# before include
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_define_exe,apt-cdf-to-spf,apt-cdf-to-spf.cpp)
$(call sdk_define_install_exe,apt-cdf-to-spf)
#
include ${sdk_makefile_post}

apt_cdf_to_spf:=../../output/amd64-pc-linux/bin/apt-cdf-to-spf

_test: _test_1

_test_1: ${apt_cdf_to_spf}
	${apt_cdf_to_spf} \
	  --cdf-file ${sdk_root}/regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf \
	  --spf-file Mapping250K_Sty.spf
