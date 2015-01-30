#
# cvs:affy/sdk/util/Makefile2 ---
#
# $Id: Makefile,v 1.31 2009-09-25 17:49:08 mspald Exp $
#

# sdk_root has a relative path to the top of the tree.
# Once set, suck in our variables and definitions.
sdk_root:=..
include ${sdk_root}/Makefile.defs

#
$(call sdk_set_lib_name,affyutil)
#$(call sdk_set_link_libs,affyutil m)

# Dont get rid of RowFile and TableFile just yet.
# sdk_cpp_lib:=$(filter-out RowFile.cpp TableFile.cpp,$(wildcard *.cpp))

# Append this to a variable which is used when running ${CXX}
# These calls define a series of rules used to build executables.
# (files which arent used in the calls will be fed to a default set
# of actions defined in "Makefile2.post".)
$(call sdk_define_exe,dump-guid,dump-guid.cpp)
$(call sdk_define_exe,util-meminfo,util-meminfo.cpp)
$(call sdk_define_exe,apt-engine-wrapper,apt-engine-wrapper.cpp)
$(call sdk_define_install_exe,apt-engine-wrapper)
$(call sdk_define_exe,apt-check-matrix,apt-check-matrix.cpp)
$(call sdk_define_exe,apt-check-mixedfile,apt-check-mixedfile.cpp)
$(call sdk_define_exe,apt-check-cel,apt-check-cel.cpp)
$(call sdk_define_exe,apt-check-chp,apt-check-chp.cpp)
$(call sdk_define_exe,apt-check-calvinchp,apt-check-calvinchp.cpp)
$(call sdk_define_exe,apt-check-textfile,apt-check-textfile.cpp)
#
$(call sdk_define_exe,test-tmpfilefactory,test-tmpfilefactory.cpp)

# While we could use ${sdk_root} to name the file,
# this allows one of the prior makefiles or calls
# to change the "post" makefile by setting the variable.
include ${sdk_makefile_post}


gTEST_MD5SUM_FILES:=test-md5sum.cpp md5sum.cpp md5.cpp
TEST_MD5SUM_EXE:=./test-md5sum-${CPUCOMSYS}
#${TEST_MD5SUM_EXE}: md5sum.h ${TEST_MD5SUM_FILES}
${TEST_MD5SUM_EXE}: md5sum.h ${TEST_MD5SUM_FILES}
	g++ -g -O0 -o $@ ${TEST_MD5SUM_FILES}


_test_md5sum: ${TEST_MD5SUM_EXE}
	${TEST_MD5SUM_EXE}

_test_meminfo: ${sdk_output_bin}/util-meminfo
	${sdk_output_bin}/util-meminfo
ifeq ($(findstring -linux,${CPUCOMSYS}),-linux)
	${sdk_output_bin}/util-meminfo ./CPPTest/input/proc-meminfo-24.txt
	${sdk_output_bin}/util-meminfo ./CPPTest/input/proc-meminfo-26.txt
endif


###

tff_exe:=${sdk_output_bin}/test-tmpfilefactory

define tff_file_clean
  -rm -rf tff-test*
endef

define tff_file_test
  @if ls tff-test* >& /dev/null ; then ${1} ; else ${2} ; fi
endef 

define tff_file_test_assert_none
  $(call tff_file_test,echo "ERR: Shouldnt have any test files."; exit 1, echo "CHECK: none ok.")
endef 

define tff_file_test_assert_some
  $(call tff_file_test,echo "CHECK: some ok",echo "ERR: Should have some test files."; exit 1)
endef 

# do all three tests
_test_tmpfilefactory: \
	_test_tmpfilefactory_debug \
	_test_tmpfilefactory_1 \
	_test_tmpfilefactory_2

# debugging should leave the files behind.
_test_tmpfilefactory_debug: ${ttf_exe}
	$(call tff_file_clean)
	$(call tff_file_test_assert_none)
#	${tff_exe} -v 3 -d 3
	${tff_exe} -v 0 -d 3
	$(call tff_file_test_assert_some)
	$(call tff_file_clean)
	@echo "ok."

# shouldnt leave files behind
_test_tmpfilefactory_1: ${ttf_exe}
	$(call tff_file_clean)
	$(call tff_file_test_assert_none)
	${tff_exe}
	$(call tff_file_test_assert_none)
	@echo "ok."

# shouldnt leave files behind but should have ouput.
_test_tmpfilefactory_2: ${ttf_exe}
	$(call tff_file_clean)
	$(call tff_file_test_assert_none)
	${tff_exe} -v 10
	$(call tff_file_test_assert_none)
	@echo "ok."
