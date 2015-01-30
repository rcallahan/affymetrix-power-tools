#
# affy/sdk/file5/Makefile ---
#
# $Id: Makefile,v 1.26 2009-10-30 19:52:38 harley Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=${hdf5_cflags}
#
$(call sdk_set_lib_name,affyfile5)
$(call sdk_set_link_libs,affyfile affyfile5 affyutil xerces-c hdf5 z m)

# File5_Object.cpp File5_Vector.cpp
#sdk_cpp_ignore:=File5_File.cpp File5_Group.cpp File5_Matrix.cpp \
#   File5_Tsv.cpp File5_util.cpp test-file5.cpp
#sdk_cpp_ignore:=File5_Tsv.cpp
#sdk_subdirs:=

# classic
$(call sdk_define_exe,test-file5,test-file5.cpp)

_clean: _clean_core

_clean_core:
	-rm core.*

# Add this before 'include ${sdk_makefile_post}'
_check_run+=_test

#
include ${sdk_makefile_post}


###
valgrind_exe:=/nfs/ss11/harley/disk/amd64-pc-linux/bin/valgrind
valgrind:=${valgrind_exe} -v --tool=memcheck --leak-check=full --show-reachable=yes --num-callers=20 
#valgrind:=

###
_test: _test_1 _test_2

test_file5_exe:=${sdk_output_bin}/test-file5

#
_test_1: ${test_file5_exe}
	-rm -f test-file5-*.file5
	${test_file5_exe} test-vector
#	${GDB} ${test_file5_exe} test-vector
#	for file5 in test-file5-*.file5* ; do ${sdk_output_bin}/h5dump $$file5; done

_test_1_gdb: ${test_file5_exe}
	gdb --args ${test_file5_exe} test-vector

_test_1_valgrind: ${test_file5_exe}
	${valgrind} ${test_file5_exe} test-vector

###
_test_2: ${test_file5_exe}
	-rm test-file5-usermeta.file5
	${test_file5_exe}
	${sdk_output_bin}/h5dump test-file5-usermeta.file5

_test_tsv_delete: ${sdk_output_bin}/test-file5
	-rm test-file5-tsv-delete.file5
	${sdk_output_bin}/test-file5 --test-tsv
	${sdk_output_bin}/h5dump -n test-file5-tsv-delete.file5

###
_test_big_vectors: ${sdk_output_bin}/test-file5
	${sdk_output_bin}/test-file5 --test-big-vectors


###
ray_exe:=${sdk_output_bin}/ray.file5.test

${ray_exe}: ray.file5.test.cpp
	g++ -O0 -g -Wall -Wno-sign-compare -Wno-non-virtual-dtor -I.. -I../output/amd64-pc-linux/include -I../output/share/include -I../../external/newmat -I/nfs/linux/pkg/lang/cppunit/1.10.2/include -o ${@} ray.file5.test.cpp -L../output/amd64-pc-linux/lib -laffysdk -laffycalvin -laffyfile5 -laffyfile -laffyutil -lasqlite3 -lhdf5 -lnewmat -lxerces-c -lpcrecpp -lpcreposix -lpcre -lm -lz

#gdb --args 
_ray: ${ray_exe}
	${ray_exe}
	h5dump /tmp/raytest.file5
