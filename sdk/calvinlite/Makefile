#
# CalvinLite/Makefile ---
#
# $Id: Makefile,v 1.2 2009-10-29 22:28:59 harley Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
include ${sdk_root}/bboard/Makefile.inc
#
$(call sdk_set_lib_name,affyfile)
$(call sdk_define_exe,apt-calvinlite-test,apt-calvinlite-test.cpp)
$(call sdk_define_exe,apt-calvinlite-util,apt-calvinlite-util.cpp)
#
include ${sdk_makefile_post}


#
apt_calvinlite_test:=${sdk_output_bin}/apt-calvinlite-test
apt_calvinlite_util:=${sdk_output_bin}/apt-calvinlite-util

_test: _test_1 _test_2 _test_3 _test_3_1 _test_4 _test_5

_test_1: ${apt_calvinlite_test}
#	-rm foo.cychp
	${apt_calvinlite_test}
#	hexdump -C foo.cychp
#	head -40 Fullerton.cychp.hex

#test_2_cel:=Fullerton.cychp
test_2_cel:=../regression-data/data/idata/cel/GenomeWideSNP_6/NA11993_GW6_C.CEL
_test_2: _test_2_1 _test_2_2

_test_2_1:
#	-rm Fullerton.cychp-*.tsv
#	${apt_calvinlite_util} --print-params ${test_2_cel}
	cp ${test_2_cel} test-2.cel
	${apt_calvinlite_util} --tsv-export-gdh test-2.cel
#	${apt_calvinlite_util} --tsv-export-all ${test_2_cel}

_test_2_2:
# prints a bunch of headers
	${apt_calvinlite_util} --print-params ${test_2_cel}
# just prints the one.
	${apt_calvinlite_util} --print-params ${test_2_cel} affymetrix-cel-cols

test_3_cel:=../regression-data/data/idata/cel/GenomeWideSNP_6/NA11993_GW6_C.CEL
_test_3: ${apt_calvinlite_util}
	cp ${test_2_cel} test-3.cel
	${apt_calvinlite_util} --print-params   test-3.cel
	${apt_calvinlite_util} --dump-segs      test-3.cel
	${apt_calvinlite_util} --tsv-export-all test-3.cel

_test_3_1: ${apt_calvinlite_util}
	${apt_calvinlite_util} --is-calvin THIS_FILE_DOES_NOT_EXIST.cel ${test_3_cel}

_test_4: _test_4_1 _test_4_2

_test_4_1: ${apt_calvinlite_util}
	cp ${test_3_cel} test_4_1.cel
	${apt_calvinlite_util} --print-params test_4_1.cel
	${apt_calvinlite_util} --change-param test_4_1.cel GenomeWideSNP_6 GenomeWideSNP_666
	${apt_calvinlite_util} --print-params test_4_1.cel.changed
	${apt_calvinlite_util} --change-param-substr test_4_1.cel Wide ReallyReallyReallyFat
	${apt_calvinlite_util} --print-params test_4_1.cel.changed
	${apt_calvinlite_util} --change-param-substr test_4_1.cel Wide ""
	${apt_calvinlite_util} --print-params test_4_1.cel.changed

_test_4_2: ${apt_calvinlite_util}
	cp NA07029_AxiomGWASHuSNP1_20090812_MCKI_T01_B11_v1.AxiomGT1.chp test_4_2.cel
#	cp ${test_3_cel} test_4.cel
# should print all
	${apt_calvinlite_util} --print-params test_4_2.cel
# should print just the date
	${apt_calvinlite_util} --print-params test_4_2.cel date
# should print date, affymetrix-array-type
	${apt_calvinlite_util} --print-params test_4_2.cel date affymetrix-array-type

#test_5_cel:=${test_3_cel}
test_5_cel:=NA07029_AxiomGWASHuSNP1_20090812_MCKI_T01_B11_v1.AxiomGT1.chp
_test_5: ${apt_calvinlite_util} ${test_3_cel}
	-rm -f TEST5*
	cp ${test_5_cel} TEST5.cel
	${apt_calvinlite_util} --tsv-export-all TEST5.cel
	${apt_calvinlite_util} --tsv-import TEST5NEW.cel TEST5.cel-*.tsv
	${apt_calvinlite_util} --tsv-export-all TEST5NEW.cel.out
# we do expect diffs in the filenames.
	for part in TEST5.cel-*.tsv ; do \
	  p=TEST5NEW.cel.out-$${part#TEST5.cel-} ; \
	  echo "diff -a -u $$part $$p" ; \
	  diff -a -u $$part $$p ; \
	  done ; true
	chmod +x TEST5NEW.cel.out-gdh.tsv
	./TEST5NEW.cel.out-gdh.tsv

##########

#valgrind_exe:=/nfs/ss11/harley/disk/${CPUCOMSYS}/bin/valgrind
valgrind_exe:=valgrind
valgrind:=${valgrind_exe} --tool=memcheck -v --leak-check=full --show-reachable=yes

v:=
v:=${valgrind}

_valgrind_test_1: ${apt_calvinlite_test} ${apt_calvinlite_util}
	${v} ${apt_calvinlite_test}
	${apt_calvinlite_util} --tsv-export test-writearray-1.calvin
#	${apt_calvinlite_util} --debug --dump-segs test_writecychpfile_2.cychp
#	${apt_calvinlite_util} --debug --dump-segs test_writecychpfile_3.cychp
#	hexdump -C test_writecychpfile_2.cychp > test_writecychpfile_2.cychp.hex
#	${apt_calvinlite_util} --tsv-export test_writecychpfile_2.cychp
#	ls -l ./test_writecychpfile_2*
##	${apt_calvinlite_util} --tsv-export test_writecychpfile_3.cychp
#	cat ./test_writecychpfile_3.cychp-gdh.tsv
#	gdb --args ${apt_calvinlite_test}
#	${apt_calvinlite_util} --tsv-export foo.cychp
