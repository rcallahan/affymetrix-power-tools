#
# /nfs/zoldo1/harley/Apt/apt1-mosaicism-1/sdk/copynumber/apt-copynumber-cyto/mosaicism/Makefile ---
#
# $Id$
#

_default:

sdk_root:=../../..
include ${sdk_root}/Makefile.defs

##########

apt_copynumber_cyto:=${sdk_output_bin}/apt-copynumber-cyto

##########

# Bens report about the error in APT:
# http://infowiki/index.php/User:Ben_Bolstad:Reports:Cyto2:201002-Mosaicism_QuantifyingDifferences_APT_vs_Prototype

##########

# compare the apt to ref with one case.

data-lib data-cels:
	mkdir ${@}

###
cf+=20NA13330_80NA07126_A05_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.CEL
cf+=20NA13330_80NA07126_A05_MosaicismStudy_51pt0_CytoScanHD_MS_20110412.CEL
cf+=20NA13330_80NA07126_A05_MosaicismStudy_51pt0_CytoScanHD_WW_20110411.CEL
cf+=20NA13330_80NA07126_B05_MosaicismStudy_50pt0_CytoScanHD_MS_20110412.CEL
cf+=20NA13330_80NA07126_B05_MosaicismStudy_50pt0_CytoScanHD_WW_20110411.CEL
cf+=20NA13330_80NA07126_B05_MosaicismStudy_50pt5_CytoScanHD_CC_20110411.CEL
cf+=20NA13330_80NA07126_C05_MosaicismStudy_50pt5_CytoScanHD_MS_20110412.CEL
cf+=20NA13330_80NA07126_C05_MosaicismStudy_50pt5_CytoScanHD_WW_20110411.CEL
cf+=20NA13330_80NA07126_C05_MosaicismStudy_51pt0_CytoScanHD_CC_20110411.CEL
cf+=30NA13330_70NA07126_A04_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.CEL
cf+=30NA13330_70NA07126_A04_MosaicismStudy_51pt0_CytoScanHD_MS_20110412.CEL
cf+=30NA13330_70NA07126_A04_MosaicismStudy_51pt0_CytoScanHD_WW_20110411.CEL
cf+=30NA13330_70NA07126_B04_MosaicismStudy_50pt0_CytoScanHD_MS_20110412.CEL
cf+=30NA13330_70NA07126_B04_MosaicismStudy_50pt0_CytoScanHD_WW_20110411.CEL
cf+=30NA13330_70NA07126_B04_MosaicismStudy_50pt5_CytoScanHD_CC_20110411.CEL
cf+=30NA13330_70NA07126_C04_MosaicismStudy_50pt5_CytoScanHD_MS_20110412.CEL
cf+=30NA13330_70NA07126_C04_MosaicismStudy_50pt5_CytoScanHD_WW_20110411.CEL
cf+=30NA13330_70NA07126_C04_MosaicismStudy_51pt0_CytoScanHD_CC_20110411.CEL

data-cels/%: | data-cels
	cp /nfs/cyto/cyto3G/CytoScanHD/Mosaicism/${@F} ${@}

_data_cel_files: ${cf:%=data-cels/%}

###
data_lib:=./data-lib

dlf+=${data_lib}/CytoScanHD_Array.cdf 
dlf+=${data_lib}/CytoScanHD_Array.chrXprobes 
dlf+=${data_lib}/CytoScanHD_Array.chrYprobes 
dlf+=${data_lib}/CytoScanHD_Array.na32.1.annot.db 
dlf+=${data_lib}/CytoScanHD_Array.r1.qca 
dlf+=${data_lib}/CytoScanHD_Array.r1.qcc 
dlf+=${data_lib}/CytoScanHD_Array.snplist.txt 

${dlf}: | ${data_lib}
	cp /nfs/cyto/cyto3G/CytoScanHD/Pipeline/libs/${@F} ${@}.tmp
	mv ${@}.tmp ${@}
data_lib_files+=${dlf}

${data_lib}/CytoScanHD_Array.na32.1.v3.REF_MODEL:
	cp /nfs/cyto/cyto3G/CytoScanHD/Pipeline/References/${@F} ${@}.tmp
	mv ${@}.tmp ${@}
data_lib_files+=${data_lib}/CytoScanHD_Array.na32.1.v3.REF_MODEL

_data_lib_files: ${data_lib} ${data_lib_files}

##########

test_1_cels:=test-1-three.cels
#test_1_cels:=test-1.cels

#valgrind_cmd=valgrind --tool=memcheck -v --show-reachable=yes --leak-check=full
#wrapper_cmd:=${valgrind_cmd}
wrapper_cmd:=time
#wrapper_cmd:=gdb --args

_test_1_apt:
	-rm -rf out-1-apt
	${wrapper_cmd} \
	${sdk_output_bin}/apt-copynumber-cyto \
	  -v 4 \
	  --cyto2 false \
	  --doDualNormalization true \
	  --keep-intermediate-data false \
	  --run-geno-qc true \
	  --qca-file ${data_lib}/CytoScanHD_Array.r1.qca \
	  --qcc-file ${data_lib}/CytoScanHD_Array.r1.qcc \
	  --snp-qc-use-contrast true \
	  --snp-qc-snp-list ${data_lib}/CytoScanHD_Array.snplist.txt \
	  --force true \
	  --adapter-type-normalization false \
	  --text-output false \
	  --cnchp-output false \
	  --cychp-output true \
	  --set-analysis-name cyhd \
	  --cdf-file ${data_lib}/CytoScanHD_Array.cdf \
	  --chrX-probes ${data_lib}/CytoScanHD_Array.chrXprobes \
	  --chrY-probes ${data_lib}/CytoScanHD_Array.chrYprobes \
	  --annotation-file ${data_lib}/CytoScanHD_Array.na32.1.annot.db \
	  --reference-input ${data_lib}/CytoScanHD_Array.na32.1.v3.REF_MODEL \
	  --out-dir ./out-1-apt-outdir \
	  --male-gender-ratio-cutoff 1.5 \
	  --female-gender-ratio-cutoff 0.9 \
	  --xx-cutoff 0.61 \
	  --xx-cutoff-high 0.95 \
	  --y-cutoff 0.58 \
	  --cn-calibrate-parameters 'calibrated-log2ratios.alpha-cn-calibrate=0.564278.alpha-X-cn-calibrate=0.619453.alpha-Y-cn-calibrate=0.494620.beta-cn-calibrate=1.beta-X-cn-calibrate=1.beta-Y-cn-calibrate=1' \
	  --local-gc-background-intensity-adjustment-method none \
	  --image-correction-intensity-adjustment-method none \
	  --wave-correction-log2ratio-adjustment-method wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=6.wave-smooth=true \
	  --log2ratio-adjustment-method log2ratio-adjustment-method-high-pass-filter.use=true \
	  --analysis genotype \
	  --analysis log2-ratio.gc-correction=false.median-autosome-median-normalization=true.median-smooth-marker-count=5 \
	  --analysis kernel-smooth.sigma_span=50 \
	  --analysis "cn-cyto2.hmmCN_state='0,1,2,3,4'.hmmCN_mu='-2,-0.5,0,0.33,0.54'.hmmCN_sigma='0.3,0.3,0.3,0.3,0.3'.hmmCN_state-X='0,1,2,3,4'.hmmCN_mu-X='-2,-0.55,0,0.36,0.57'.hmmCN_sigma-X='0.3,0.3,0.3,0.3,0.3'.hmmCN_state-Y='0,1,2,3,4'.hmmCN_mu-Y='-1.4,-0.5,0,0.3,0.51'.hmmCN_sigma-Y='0.3,0.3,0.3,0.3,0.3'.diagonal-weight-Y=0.995.mapd-weight-Y=0.2.min-segment-size-Y=5.hmm-confidence-weight-Y=0.6.diagonal-weight=0.995.mapd-weight=0.2.shrink=true.min-segment-size=5.hmm-confidence-weight=0.6.diagonal-weight-X=0.995.mapd-weight-X=0.2.min-segment-size-X=5.hmm-confidence-weight-X=0.6" \
	  --analysis cn-cyto2-gender.cutoff=0.5 \
	  --analysis cn-segment \
	  --analysis "mosaicism.gains-boundries='0.08764945,0.15380349,0.21465931,0.27100300'.losses-boundries='-0.08293345,-0.17551812,-0.28048196,-0.40165383'.marker-bandwidth=6000.confidence-window=251" \
	  --cel-files ${test_1_cels} >& out-1-apt.output
#	mv -i ?0NA*.tsv out-1-apt-outdir
#|& tee 

# --analysis allelic-difference-CytoScan.outlier-trim=3.0.step=20.window=100.point-count=128.bandwidth=0.25.cutoff=0.05.clean-threshold=0.35.symmetry=true
# --analysis loh-segment
# --analysis cn-neutral-loh
# --analysis lohCytoScan.lohCS_errorrate=0.05.lohCS_beta=0.001.lohCS_alpha=0.01.lohCS_separation=1000000.lohCS_nMinMarkers=10.lohCS_NoCallThreshold=0.05.lohCS_minGenomicSpan=1000000

_test_1_apt_tsv_export:
	cd ./out-1-apt-outdir && \
	  ../../../../output/x86_64-intel-linux/bin/apt-calvinlite-util --tsv-export-all *.cychp
# less ./out-1-apt-outdir/*.cyhd.cychp-g03-s01.tsv

_test_1_gdb_core:
	gdb ${sdk_output_bin}/apt-copynumber-cyto core.*

_test_1_cmp:
	emacs -nw \
	  out-1-R-outdir/20NA13330_80NA07126_A05_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.cyhd.cychp-01--03-clean2 \
	  out-1-apt-outdir/20NA13330_80NA07126_A05_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.CEL-chr1.segclean-2.tsv

##########

_test_1_R:
	-mkdir -p out-1-R-outdir/
#	R CMD BATCH GetPrototypeSegments.R - |& tee reference.Rout
	(time R --vanilla --verbose < GetPrototypeSegments.R) >& out-1-R-outdir/out-1-R.Rout
# |& tee test-1-R-outdir/test-1-R.Rout

_test_1_R_diff:
	for seg in test-1-R-outdir/*.protoMosaicSeg ; \
	do \
	  seg_base=$$(basename $${seg}) ; \
	  echo "### $${seg_base}" ; \
	  diff -u $${seg} /nfs/cyto/usr/bbolst/Cyto3G/CytoScanHD/DataAnalysis/Mosaic2.7M/Prototype/$${seg_base} ; \
	done
