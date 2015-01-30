#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.3 2009-09-25 17:48:53 mspald Exp $
#

#
sdk_root:=../../../..
include ${sdk_root}/Makefile.defs
include ../../../Makefile.inc
#
# Not all the files have
#
testfiles=$(wildcard *.cpp)
#
$(call sdk_define_exe,calvin-compare-cel,$(testfiles))
#
include ${sdk_makefile_post}
