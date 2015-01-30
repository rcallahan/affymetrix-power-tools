#
# affy/sdk/file5/example/Makefile ---
#
# $Id: Makefile,v 1.5 2009-04-22 21:13:54 awilli Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_cflags+=${hdf5_cflags}
#
$(call sdk_set_lib_name,affyfile5)
$(call sdk_set_link_libs,affyfile affyfile5 affyutil xerces-c hdf5 z m asqlite3)
# classic
$(call sdk_define_exe,example-file5-1,example-file5-1.cpp)
$(call sdk_define_exe,simple-file5-example,simple-file5-example.cpp)
$(call sdk_define_exe,apt-file5-to-sqlite, apt-file5-to-sqlite.cpp)

include ${sdk_makefile_post}

EF_EXE:=${sdk_output_bin}/example-file5-1

_test: ${EF_EXE} _test_1

_test_1: ${EF_EXE}
	@-rm test-1.a5
	${EF_EXE} --create-file test-1.a5
	# this will be ok
	${EF_EXE} --open-file test-1.a5
	# this will fail (already exists)
	-${EF_EXE} --create-file test-1.a5

_test_2:
	@-rm test-2.a5
	${EF_EXE} -write-tsv test-2.a5 --num-rows 1000
	${EF_EXE} -read-tsv test-2.a5 --num-rows 10
	${EF_EXE} -read-tsv-vec test-2.a5 
