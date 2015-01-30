#
# affy/sdk/chipstream/apt-convert-priors/Makefile ---
#
# $Id: Makefile,v 1.4 2009-04-22 21:13:53 awilli Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
# sdk_cpp_ignore:=apt-convert-priors.cpp
$(call sdk_define_exe,convert-priors,apt-convert-priors.cpp)
#
include ${sdk_makefile_post}

### tests
convert_priors_exe:=${sdk_output_bin}/convert-priors

# so we dont have to wait
_quick_build:
	cd .. ; make _build_lib 
	make

#
_test: _test_1 _test_2

_test_1: ${convert_priors_exe}
	${convert_priors_exe} --format-ver=1 test-posteriors.txt
	-diff -u test-posteriors.txt test-posteriors.txt.v1.txt | less

_test_2: ${convert_priors_exe}
	${convert_priors_exe} --format-ver=2 test-posteriors.txt

snp_a5:=../../regression-data/data/copy-number/GenomeWideSNP_6/doCopynumberWorkflowCopynumberRun/CN5.snp-posteriors.a5

_test_3: ${convert_priors_exe}
	${convert_priors_exe} --format-ver=2 -o test-3.tsv --tsv-name CN5.snp-posteriors ${snp_a5}
