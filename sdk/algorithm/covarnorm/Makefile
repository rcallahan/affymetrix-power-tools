sdk_root:=../..
include ${sdk_root}/Makefile.defs
# to ignore scratch file do this
sdk_cpp_ignore:= tmp.cpp
#                  
#                  
$(call sdk_define_exe,test-toy-norm,toy.norm.cpp)
#        
include ${sdk_makefile_post}

test_out_tmp:=./test-toy-norm.tmp
_check: ${sdk_output_bin}/test-toy-norm
	${sdk_output_bin}/test-toy-norm >& ${test_out_tmp}
	diff -w -u ${test_out_tmp} expected/test-toy-norm.out && echo "ok."
