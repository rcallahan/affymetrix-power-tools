#
# cvs:affy/sdk/chipstream/Makefile ---
#
# $Id: Makefile,v 1.4 2009-09-25 17:48:51 mspald Exp $
#

# before include
sdk_root:=..
include ${sdk_root}/Makefile.defs
#
#
include ${sdk_makefile_post}


### converter

P_DIR:= /nfs/ss11/harley/Exact/nfs-python/lib/birdseed/GenomeWideSNP_6/20070327

_convert:
	./apt-birdseed-priors-to-tsv -p ${P_DIR}/SNP6.0.priors.txt -o ${P_DIR}/SNP6.0.priors.tsv

_convert_test:
	./apt-birdseed-priors-to-tsv -p CPPTest/input/priors.txt -o priors.tsv && cat priors.tsv
