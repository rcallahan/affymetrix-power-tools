#
# affy/sdk/file5/apt-file5-util/Makefile ---
#
# $Id: Makefile,v 1.24 2009-09-10 00:26:17 harley Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
$(call sdk_set_lib_name,affyfile5)
$(call sdk_set_link_libs,affyfile affyfile5 affyutil xerces-c hdf5 z m)
# classic
$(call sdk_define_exe,apt-file5-util,apt-file5-util.cpp)
$(call sdk_define_install_exe,apt-file5-util)

#
include ${sdk_makefile_post}

#
apt_file5_util:=${sdk_output_bin}/apt-file5-util

#
h5dump:=${sdk_output_bin}/h5dump

#
valgrind_exe:=/nfs/ss11/harley/disk/amd64-pc-linux/bin/valgrind
valgrind_mem:=${valgrind_exe} --tool=memcheck --leak-check=full --show-reachable=yes -v
valgrind_mem_gdb:=${valgrind_mem} --db-attach=yes
valgrind_massif:=${valgrind_exe} --tool=massif --depth=100 -v

# GDB:=gdb --args

###
test_gen_file:=test_gen_file.tsv

#
TEST_TSV_PGF1:=data-test-1.pgf
${TEST_TSV_PGF1}:
	cp ../../file/TsvFile/data-test-1.pgf .

#
TEST_TSV_SNP6:=GenomeWideSNP_6.spf
${TEST_TSV_SNP6}:
	cp ../../regression-data/data/lib/GenomeWideSNP_6/${TEST_TSV_SNP6} .

TEST_TSV_SNP6_SMALL:=GenomeWideSNP_6.spf.small
${TEST_TSV_SNP6_SMALL}: ${TEST_TSV_SNP6}
	head -10000 ${TEST_TSV_SNP6} > ${@}

# which file?
TEST_TSV:=${TEST_TSV_PGF1}
#TEST_TSV:=${TEST_TSV_SNP6_SMALL}

###
_test: ${TEST_TSV} _test_1 _test_2

_test_1: ${apt_file5_util}
	${apt_file5_util} --help

_test_2: ${apt_file5_util} _test_2_1 _test_2_2

_test_2_1:
	-rm ${TEST_TSV}5
	${GDB} ${apt_file5_util} --to-tsv5 ${TEST_TSV} --internal-name ${TEST_TSV} -o ${TEST_TSV}5

_test_2_2:
	${GDB} ${apt_file5_util} --to-tsv ${TEST_TSV}5 --internal-name ${TEST_TSV} -o ${TEST_TSV}.tsv-out

###
_test_3: _test_3_1 _test_3_2

_test_3_1: ${apt_file5_util} ${TEST_TSV}
	${apt_file5_util} --to-tsv5 ${TEST_TSV} --internal-name ${TEST_TSV} -o ${TEST_TSV}5
	${h5dump} -A ${TEST_TSV}5 | grep -e DATASET -e GROUP | grep "_user_meta"

_test_3_2: ${apt_file5_util} ${TEST_TSV}
	${apt_file5_util} --no-metadata --to-tsv5 ${TEST_TSV} --internal-name ${TEST_TSV} -o ${TEST_TSV}5
#	${h5dump} -A ${TEST_TSV}5 | grep -e DATASET -e GROUP 
	${h5dump} -A ${TEST_TSV}5 | grep -e DATASET -e GROUP | grep -v "_user_meta"

###
_test_5: _test_5_1 _test_5_2 _test_5_check

_test_5_1:
	-rm ${data_test}.pgf5 ${data_test}-out.pgf5
	${GDB} ${apt_file5_util} --to-tsv5 ${data_test}.pgf
	${h5dump} ${data_test}.pgf5

_test_5_2:
#	cp ${data_test}.pgf5-ref ${data_test}-out.pgf5
	cp ${data_test}.pgf5 ${data_test}-out.pgf5
	${GDB} ${apt_file5_util} --to-tsv ${data_test}-out.pgf5

_test_5_1_valgrind:
	-rm ${data_test}.pgf5 ${data_test}-out.pgf5
	 ${valgrind_mem} \
	  ${apt_file5_util} --to-tsv5 ${data_test}.pgf

_test_5_check:
	diff ${data_test}.pgf ${data_test}-out.pgf

###
_test_gen_file: ${test_gen_file}

${test_gen_file}: ${apt_file5_util} Makefile
	time ${apt_file5_util} --gen-test-file ${test_gen_file} --gen-test-file-size=10000



# time to write:
# (10+10) x 100000 => 30 sec
# (10+0)  x 100000 => 13
# (0+10)  x 100000 => 17

# space for strings:
#
# (10+0) x 1000000 => tsv format  11000100_B (11MB)  (time to make= 2.35s)
#                  => hdf5 format 49709184_B (48MB)  (time to convert= 20s)

# == tested on ostia:

# chunk size: 1000
# Same file but with fixed length strings (20B) 20889040_B
# fixed length ( 20B) with compression           3579462_B
# fixed length (256B) with compression          12598913_B  (50sec)

# all string data is written, but h5dump truncs at the first null.

#
# CS=10,000 fixed length (256B) with compression  => too long to run!
# CS=100 =>  13974649_B (35sec)
# CS=100, A-Z, no compress => 257298272_B (19sec)
# CS=100, A-Z, compress =>      2586576_B (21sec)
# CS=100, strval, compress =>   4404904_B (17sec)
# CS=100, cleared & strval, 1KB, no compress  => 1025298272_B (31sec)
# (above) compress                            =>    6334984_B (52sec)
# (above) 512B                                =>    5462982_B (44sec)
#         256                                 =>    4430640_B (18sec)
#         128                                 =>    4188112_B (15sec)
#          64                                 =>    3950544_B (17sec)


# CS=100, cleared & strval, 256B, no compress =>  257298272_B (17sec)
# (above) compress                            =>    4430640_B (20sec)
#       , uncleared & strval, 256B,           =>   14090092_B (38sec)


#test_gen_file:=test_gen_file-5000.tsv
_test_6: ${test_gen_file}
#${valgrind_mem}
	${apt_file5_util} --to-tsv5 ${test_gen_file} --internal-name ${test_gen_file}
#	time ${apt_file5_util} --to-tsv5 ${test_gen_file}
#	h5dump ${test_gen_file}5

# gdb --args
CNREF_A5:=../../regression-data/data/copy-number/cppunit/input/CNReference.a5
_test_7: ${apt_file5_util}
	${apt_file5_util} --to-tsv ${CNREF_A5} --internal-name CN5/CN5.calls -o CN5-calls.tsv
	${apt_file5_util} --to-tsv ${CNREF_A5} --internal-name CN5/CN5.snp-posteriors -o CN5-snp-posteriors.tsv

###

# This test is for the user headers and appending/replacing datasets.

test_8:=test-8
test_8_tsv:=${test_8}.tsv
test_8_a5:=${test_8}.a5
data_test_1_pgf:=../../file/TsvFile/data-test-1.pgf

_test_8: _test_8_clean ${test_8_tsv} _test_8_1 _test_8_2 _test_8_diff _test_8_3 _test_8_4

_test_8_clean:
	@-rm ${test_8_a5} ${test_8_tsv} ${test_8_tsv}.out

${test_8_tsv}: ${apt_file5_util} Makefile
	${valgrind_mem} \
	${apt_file5_util} --gen-test-file ${test_8_tsv} --gen-test-file-levels=3 --gen-test-file-size=100 \
	>& ${@}.valgrind

_test_8_1: ${apt_file5_util}
	${valgrind_mem} \
	${apt_file5_util} --to-tsv5 ${test_8_tsv} --internal-name ${test_8} -o ${test_8_a5} \
	>& ${@}.valgrind

_test_8_2: ${apt_file5_util}
	${valgrind_mem} \
	${apt_file5_util} --to-tsv ${test_8_a5} --internal-name ${test_8} -o ${test_8_tsv}.out \
	>& ${@}.valgrind

_test_8_diff:
	diff -u ${test_8_tsv} ${test_8_tsv}.out >& ${@}.diff

# append to add
_test_8_3: ${apt_file5_util}
	${apt_file5_util} --to-tsv5 ${test_8_tsv} --internal-name append -o ${test_8_a5} --append
#	${h5dump} -n ${test_8_a5}

# append to replace
_test_8_4: ${apt_file5_util}
	${apt_file5_util} --to-tsv5 ${data_test_1_pgf} --internal-name append -o ${test_8_a5} --append
#	${h5dump} -n ${test_8_a5}



###

DMET_Plus_r3_genomic_ref_a5:=/nfs/dmet/usr/awilli/dmet-build-r3-lib-files/DMET_Plus.r3.genomic.ref.a5
test9:=test9
test9_tsv:=${test9}.tsv
test9_file5:=${test9}.tsv5

_test9: ${apt_file5_util} _test_9_1

${test9_tsv}:
	${apt_file5_util} --gen-test-file ${@} --gen-test-str-cols 2 --gen-test-int-cols 0 --gen-test-double-cols 10

${test9_file5}: ${test9_tsv}
	${apt_file5_util} --to-tsv5 ${test9_tsv} --internal-name ${test9}

_test_9_1: ${test9_file5}
	${apt_file5_util} -i ${test9_file5} --internal-name ${test9} --find-col test-s-000 --find-val str-000300
	${apt_file5_util} -i ${test9_file5} --internal-name ${test9} --find-col test-s-000 --find-val notfound

_test_9_dump: ${apt_file5_util}
	${apt_file5_util} \
	  --to-tsv \
	  --internal-name ProbesetGenotypeEngine/dmet.snp-posteriors \
	  -o test-9.tsv \
	  -p 3 \
	  ${DMET_Plus_r3_genomic_ref_a5} \
# now reverse the process
	${apt_file5_util} --to-tsv5 -o test-9.a5 test-9.tsv

_test_9_2:
	${apt_file5_util} \
	-i /nfs/dmet/usr/awilli/dmet-build-r3-lib-files/DMET_Plus.r3.genomic.ref.a5 \
	--internal-name ProbesetGenotypeEngine/dmet.snp-posteriors \
	--find-col id \
	--find-val AM_10004-0-1

###
_test_10: _test_10_1 _test_10_2 _test_10_3

_test_10_1: ${apt_file5_util}
	${apt_file5_util} --line-count --file-name ${CNREF_A5} --internal-name CN5/CN5.snp-posteriors
_test_10_2: ${apt_file5_util}
	${apt_file5_util} --line-count --file-name ${CNREF_A5}-xxx --internal-name CN5/CN5.snp-posteriors
_test_10_3: ${apt_file5_util}
	${apt_file5_util} --line-count --file-name ${CNREF_A5} --internal-name CN5/CN5.snp-posteriors-xxx
