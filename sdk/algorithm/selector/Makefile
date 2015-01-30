#
# affy/sdk/algorithm/selector/Makefile ---
#
# $Id: Makefile,v 1.5 2009-09-25 17:48:50 mspald Exp $
#

sdk_root:=../..
include ${sdk_root}/Makefile.defs

# to ignore a scratch file do this
sdk_cpp_ignore:=tmp.cpp
#                  
#                  
$(call sdk_define_exe,test-ToySelect,ToySelect.cpp)
#        
include ${sdk_makefile_post}

_check: ${sdk_output_bin}/test-ToySelect
	${sdk_output_bin}/test-ToySelect >& test-ToySelect.out
	diff -w test-ToySelect.out expected/test-ToySelect.out && echo "ok."
