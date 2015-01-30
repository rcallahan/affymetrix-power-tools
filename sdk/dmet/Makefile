#
# cvs:affy/sdk/chipstream/Makefile ---
#
# $Id: Makefile,v 1.9 2009-05-27 20:35:53 awilli Exp $
#

# before include
sdk_root:=..
include ${sdk_root}/Makefile.defs
$(call sdk_define_exe,apt-dmet-copynumber,apt-dmet-copynumber.cpp)
$(call sdk_define_install_exe,apt-dmet-copynumber)
$(call sdk_define_exe,apt-dmet-genotype,apt-dmet-genotype.cpp)
$(call sdk_define_install_exe,apt-dmet-genotype)
#
sdk_subdirs=regression
include ${sdk_makefile_post}
