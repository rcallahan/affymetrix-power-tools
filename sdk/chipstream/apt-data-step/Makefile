#
# cvs:affy/sdk/util/apt-data-step/Makefile ---
#
# $Id: Makefile,v 1.3 2009-10-19 20:58:57 csugne Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs
sdk_cpp:=$(wildcard *.cpp)
#
#
#
$(call sdk_define_exe,apt-data-step, apt-data-step.cpp)
$(call sdk_define_exe,cel-to-data-store, cel-to-data-store.cpp)
$(call sdk_define_exe,apt-probeset-summarize-stages,apt-probeset-summarize-stages.cpp)
$(call sdk_define_exe,apt-probeset-genotype-stages,apt-probeset-genotype-stages.cpp)
#$(call sdk_define_install_exe,apt-data-subset)

include ${sdk_makefile_post}
