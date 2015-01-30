#
# cvs:affy/sdk/util/Makefile2 ---
#
# $Id: Makefile,v 1.7 2009-09-25 17:49:05 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore=$(wildcard *.cpp)
# Put these in order in hopes they will be run in this order.
to_test_cpp:=Satuint16Test.cpp VectormapTest.cpp NormalizationTest.cpp
$(call sdk_define_check_exe,test-normalization,${to_test_cpp})
#
include ${sdk_makefile_post}

# Local Variables:
# mode: makefile
# End:
