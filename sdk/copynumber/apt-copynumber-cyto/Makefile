#
# affy/sdk/chipstream/apt-probeset-genotype/Makefile ---
#
# $Id: Makefile,v 1.2 2009-04-22 21:13:54 awilli Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs

# dont include "mosaicism" yet.
sdk_subdirs:=regression

# classic
$(call sdk_define_exe,apt-copynumber-cyto,apt-copynumber-cyto.cpp)
$(call sdk_define_install_exe,apt-copynumber-cyto)

include ${sdk_makefile_post}
