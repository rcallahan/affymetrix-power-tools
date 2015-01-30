#
# sdk/chipstream/CPPTest/Makefile ---
#
# $Id: Makefile,v 1.15 2009-09-25 17:48:56 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
ignore_files=\
  PmAdjusterTest.cpp \
  SparseMartTest.cpp
#
$(call sdk_define_check_exe,test-chipstream,$(filter-out $(ignore_files), $(wildcard *.cpp)))
#
include ${sdk_makefile_post}

# Local Variables:
# mode: makefile
# End:
