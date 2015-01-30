#
# cvs:affy/sdk/dabg/Makefile ---
#
# $Id: Makefile,v 1.23 2009-09-25 17:48:57 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
sdk_cpp:=$(wildcard *.cpp)
#
$(call sdk_define_exe,apt-summary-vis,apt-summary-vis.cpp SummaryVis.cpp)
$(call sdk_define_install_exe,apt-summary-vis)
$(call sdk_define_exe,apt-probeset-summarize,apt-probeset-summarize.cpp)
$(call sdk_define_install_exe,apt-probeset-summarize)
#
include ${sdk_makefile_post}
