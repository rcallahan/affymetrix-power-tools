#
# cvs:affy/sdk/file/TsvFile/Makefile ---
#
# $Id: Makefile,v 1.56 2009-09-16 19:11:13 harley Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile)
#$(call sdk_set_link_libs,affyfile affyutil m)
#
sdk_cpp_ignore:=tsv-util.cpp
#
$(call sdk_define_exe,test-tsv,tsv-test.cpp)
$(call sdk_define_exe,tsv-example,tsv-example.cpp)
$(call sdk_define_exe,tsv-util,tsv-util.cpp)
$(call sdk_define_exe,apt-tsv-join,apt-tsv-join.cpp TsvJoin.cpp)
$(call sdk_define_install_exe,apt-tsv-join)
$(call sdk_define_exe,apt-dump-pgf,apt-dump-pgf.cpp DumpPgf.cpp)
$(call sdk_define_install_exe,apt-dump-pgf)
#
_check_run+=_tsv_test


#
include ${sdk_makefile_post}

###

valgrind_exe:=/nfs/ss11/harley/disk/${CPUCOMSYS}/bin/valgrind
valgrind:=${valgrind_exe} --tool=memcheck -v --leak-check=full --show-reachable=yes

###

tsv_util:=${sdk_output_bin}/tsv-util

_tsv_test: ${sdk_output_bin}/test-tsv
	${sdk_output_bin}/test-tsv
	@diff -w data-test-4.tsv data-test-4.out && echo "diff ok."

_tsv_test_valgrind: ${sdk_output_bin}/test-tsv
	${valgrind} ${sdk_output_bin}/test-tsv

_tsv_time: ${sdk_output_bin}/tsv-example
	time ${sdk_output_bin}/tsv-example /var/tmp/HuEx-1_0-st-v2.pgf

_tsv_gprof: ${sdk_output_bin}/tsv-example
	${sdk_output_bin}/tsv-example /var/tmp/HuEx-1_0-st-v2.pgf
	gprof ${sdk_output_bin}/tsv-example > gprof.out

_tsv_example: ${sdk_output_bin}/tsv-example
	${sdk_output_bin}/tsv-example ./data-test-1.pgf

_tsv_util_test: ${tsv_util} \
	_tsv_util_test_help \
	_tsv_util_test_headers \
	_tsv_util_test_linecount \
	_tsv_util_test_paste

_tsv_util_test_help:
	${tsv_util} --help

_tsv_util_test_headers:
	${tsv_util} --headers data-test-1.pgf

_tsv_util_test_linecount:
	${tsv_util} --linecount data-test-1.pgf

_tsv_util_test_paste_1: ${tsv_util}
	${tsv_util} --paste paste.tsv data-test-1.pgf data-test-2.csv data-test-3.tsv
_tsv_util_test_paste_2: ${tsv_util}
	${tsv_util} --seg-lines 1 --paste paste.tsv data-test-1.pgf data-test-2.csv data-test-3.tsv
_tsv_util_test_paste_3: ${tsv_util}
	-rm paste.tsv
	${tsv_util} --paste paste.tsv --check-col probeset_id input/input-0?.tsv

_tsv_util_test_paste: \
	_tsv_util_test_paste_ok1 \
	_tsv_util_test_paste_ok2 \
	_tsv_util_test_paste_bad1 \
	_tsv_util_test_paste_bad2 \
	_tsv_util_test_paste_bad3

tsv_util_paste_test_exe:=${tsv_util} --paste paste.tsv --check-col probeset_id input/input-0?.tsv

##
_tsv_util_test_paste_ok1: ${tsv_util}
	-rm paste.tsv
	${tsv_util} --paste paste.tsv input/input-0?.tsv || exit 1
	[ -f paste.tsv ] || ( echo "fail ${@}" && exit 0)
_tsv_util_test_paste_ok2: ${tsv_util}
	-rm paste.tsv
	${tsv_util_paste_test_exe} || exit 1
	[ -f paste.tsv ] || ( echo "fail ${@}" && exit 0)
_tsv_util_test_paste_bad1: ${tsv_util}
	-rm paste.tsv
	-${tsv_util_paste_test_exe} input/input-bad-col.tsv
	-[ ! -f paste.tsv ] && echo "fail ${@}" && exit 1
_tsv_util_test_paste_bad2: ${tsv_util}
	-rm paste.tsv
	-${tsv_util_paste_test_exe} input/input-bad-short.tsv
	-[ ! -f paste.tsv ] && echo "fail ${@}" && exit 1
_tsv_util_test_paste_bad3: ${tsv_util}
	-rm paste.tsv
	-${tsv_util_paste_test_exe} input/input-bad-val.tsv
	-[ ! -f paste.tsv ] && echo "fail ${@}" && exit 1

_tsv_util_test_paste_valgrind:
	${valgrind} ${tsv_util} --paste paste.tsv data-test-1.pgf data-test-2.csv data-test-3.tsv

###

tsv_util_paste:=time ${tsv_util} -v 2 --key-col 0 --paste

_tsv_util_test_paste_passes: \
	_tsv_util_test_paste_20 \
	_tsv_util_test_paste_21 _tsv_util_test_paste_22 _tsv_util_test_paste_2122_diff

test-paste-inputs:
	mkdir ${@}
	for i in $$(seq 100) ; do \
	  ii=$$(printf "%04d" $$i) ; \
	  cp input/input-01.tsv ${@}/input-$${ii}.tsv ; \
	done

_tsv_util_test_paste_20: | ${tsv_util} test-paste-inputs
	-rm paste-20.tsv
	${tsv_util_paste} paste-20.tsv --max-paste-fh 2 \
	  input/input-0?.tsv input/input-0?.tsv || exit 1
	[ -f paste-20.tsv ] || ( echo "fail ${@}" && exit 0)

_tsv_util_test_paste_21: | ${tsv_util} test-paste-inputs
	-rm paste-21.tsv
	${tsv_util_paste} paste-21.tsv --max-paste-fh 5 \
	  test-paste-inputs/input-*.tsv || exit 1
	[ -f paste-21.tsv ] || ( echo "fail ${@}" && exit 0)

_tsv_util_test_paste_22: | ${tsv_util} test-paste-inputs
	-rm paste-22.tsv
	${tsv_util_paste} paste-22.tsv \
	 --paste paste-22.tsv \
	  test-paste-inputs/input-*.tsv || exit 1
	[ -f paste-22.tsv ] || ( echo "fail ${@}" && exit 0)

_tsv_util_test_paste_2122_diff:
	diff -u paste-21.tsv paste-22.tsv

###

# 'data-test-5.pgf' is missing; using 'data-test-5.tsv' for now.
define tsv_diff_test
	${tsv_util} ${1} ; \
	rv=$$? ; \
	echo '$$rv=' $$rv ; \
	if [ "$$rv" == "${2}" ] ; then exit 0 ; else exit 1 ; fi
endef

_tsv_util_diff: ${tsv_util} _tsv_util_diff_1 _tsv_util_diff_2 _tsv_util_diff_3

_tsv_util_diff_1: ${tsv_util}
	$(call tsv_diff_test,--diff-print-same=false --diff data-test-1.pgf data-test-1.pgf,0)
	$(call tsv_diff_test,--diff-print-max=10     --diff data-test-1.pgf data-test-5.tsv,0)

_tsv_util_diff_2: ${tsv_util}
	$(call tsv_diff_test,--diff-print-linenums=0 --diff-print-format=2 --diff data-test-1.pgf data-test-5.tsv,0)
	$(call tsv_diff_test,--diff-print-linenums=0 --diff-print-format=2 --diff data-test-1.pgf data-test-5.tsv,0)

_tsv_util_diff_3: ${tsv_util}
	$(call tsv_diff_test,--diff-headers data-test-5.tsv data-test-6.tsv,0)
	$(call tsv_diff_test,--diff-headers data-test-6.tsv data-test-5.tsv,0)

#
_tsv_util_dup: _tsv_util_dup_empty _tsv_util_dup_5 _tsv_util_dup_6

_tsv_util_dup_empty:
	-${tsv_util} --print-duplicate-headerkeys input/input-empty.tsv
	@if [ $$(${tsv_util} --print-duplicate-headerkeys input/input-empty.tsv | wc -l) != "0" ] ; \
	then \
	  echo "Shouldnt have duplicate headers." ; exit 1; \
	fi

_tsv_util_dup_5:
	-${tsv_util} --print-duplicate-headerkeys data-test-5.tsv
	@if [ $$(${tsv_util} --print-duplicate-headerkeys data-test-5.tsv | wc -l) != "0" ] ; \
	then \
	  echo "Shouldnt have duplicate headers." ; exit 1; \
	fi
_tsv_util_dup_6:
	-${tsv_util} --print-duplicate-headerkeys data-test-6.tsv
	@if [ $$(${tsv_util} --print-duplicate-headerkeys data-test-6.tsv | wc -l) != "5" ] ; \
	then \
	  echo "Should have duplicate headers." ; exit 1; \
	fi

### for the lack of a better place to build this in the apt tree
### this is parked here.
### to read about it: 
### http://cvs/viewvc/analysis/selection/qc/compare/apt-merge-shortmers/Makefile?view=markup
###
_apt_merge_shortmers: apt-merge-shortmers

apt-merge-shortmers:
	cvs co -d ${@} analysis/selection/qc/compare/apt-merge-shortmers 
