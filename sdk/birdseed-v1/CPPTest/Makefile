#
# cvs:affy/sdk/dabg/Makefile ---
#
# $Id: Makefile,v 1.1 2007-09-11 20:10:54 awilli Exp $
#

# add these test to the check target
## TODO: This is busted. The later two fail because of small numeric differences
##       and the first fails (seg fault) on solaris. 
_check_run+=\
 _check_convert \
 _check_birdseed_nogender_tsv \
 _check_birdseed_nogender

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs

#
#sdk_cpp_cflags+=-I${sdk_root}/broadutil -I${sdk_root}/birdseed
#
# This doesnt doe what we want... it runs without args.
#$(call sdk_define_check_exe,birdseed,$(wildcard *.cpp))
#
$(call sdk_define_exe,birdseed-v1,birdseed.cpp)
# test and converter program
$(call sdk_define_exe,birdseed-priors-util-v1,birdseed-priors-util.cpp)
#
include ${sdk_makefile_post}

##########

GDB_CMD:=gdb --args
I:=./input
O:=./output
E:=./expected

#
${O}:
	mkdir $@

#
_check_convert: ${sdk_output_bin}/birdseed-v1 | ${O}
	@-rm -f ${O}/priors-convert.{bin,tsv,txt}
	${sdk_output_bin}/birdseed-v1 --priors-text ${I}/priors.txt --write-priors-binary ${O}/priors-convert.bin
	${sdk_output_bin}/birdseed-v1 --priors-text ${I}/priors.txt --write-priors-tsv    ${O}/priors-convert.tsv 
	${sdk_output_bin}/birdseed-v1 --priors-text ${I}/priors.txt --write-priors-text   ${O}/priors-convert.txt 

# cribbed from birdseedtest.py
_check_birdseed_nogender: ${sdk_output_bin}/birdseed-v1 | output
	@-rm -f ${O}/nogender.calls ${O}/nogender.confidences 
	${sdk_output_bin}/birdseed-v1 -v 1 \
	  --priors-text input/priors.txt \
	 -c 1270.366651 --chrX-snps ${I}/BI_SNP.chrx \
	 ${I}/nogender.intensities \
	 ${O}/nogender.calls ${O}/nogender.confidences
	${sdk_output_bin}/apt-matrix-diff -e 0.0001 -l 1 -c 1 expected/nogender.calls output/nogender.calls
	${sdk_output_bin}/apt-matrix-diff -e 0.0001 -l 1 -c 1 expected/nogender.confidences output/nogender.confidences


# VALGRIND:=valgrind  --tool=memcheck --db-attach=yes --leak-check=full

# cribbed from birdseedtest.py
_check_birdseed_nogender_tsv: ${sdk_output_bin}/birdseed-v1 | output
	@-rm ${O}/nogender.confidences.tsv ${O}/nogender.confidences.tsv ${O}/priors.tsv
	../../birdseed-dev/apt-birdseed-priors-to-tsv -p ${I}/priors.txt -o ${O}/priors.tsv
	${VALGRIND} ${sdk_output_bin}/birdseed-v1 -v 1 \
	  --priors-tsv output/priors.tsv \
	 -c 1270.366651 --chrX-snps ${I}/BI_SNP.chrx \
	 ${I}/nogender.intensities \
	 ${O}/nogender.calls.tsv ${O}/nogender.confidences.tsv
	${sdk_output_bin}/apt-matrix-diff -e 0.0001 -l 1 -c 1 expected/nogender.calls output/nogender.calls.tsv
	${sdk_output_bin}/apt-matrix-diff -e 0.0001 -l 1 -c 1 expected/nogender.confidences output/nogender.confidences.tsv

_check_apt_convert:
	@-rm ${O}/priors.tsv
	../../birdseed-dev/apt-birdseed-priors-to-tsv -ct chocolate -p ${I}/priors.txt -o ${O}/priors.tsv


_jhg: ${sdk_output_bin}/birdseed-priors-util-v1
