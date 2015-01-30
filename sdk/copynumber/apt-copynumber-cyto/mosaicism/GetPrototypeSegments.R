#
# sdk/copynumber/apt-copynumber-cyto/mosaicism/GetPrototypeSegments.R ---
#

library(bitops,lib.loc="/nfs/cyto/usr/bbolst/R-Packs")
library(caTools,lib.loc="/nfs/cyto/usr/bbolst/R-Packs")

source("MosaicAlgorithmCode4.R")

library(affyio,lib.loc="/nfs/cyto/usr/bbolst/R-Packs")

Read.CYCHP <- function(filename) {
  return(.Call("Read_Generic_R_List",filename))
}

## 
out_dir <- "out-1-R-outdir"

## these are exact boundaries as provided to APT
gains.boundaries <- c(0.08764945, 0.15380349, 0.21465931, 0.27100300)
losses.boundaries <- c(-0.08293345, -0.17551812, -0.28048196, -0.40165383)

# should read from "test-1.cels"
#filenames <- dir("../",full=TRUE,pattern="cychp$")
#filenames <- dir("test-1-apt/",full=TRUE,pattern="cychp$")
filenames <- dir("/nfs/cyto/usr/bbolst/Cyto3G/CytoScanHD/DataAnalysis/Mosaic2.7M",full=TRUE,pattern="cychp$")
#filenames <- dir("./test-1-apt-outdir",full=TRUE,pattern="cychp$")

print(filenames)

#
filenames <- c(
"/nfs/cyto/usr/bbolst/Cyto3G/CytoScanHD/DataAnalysis/Mosaic2.7M/20NA13330_80NA07126_A05_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.cyhd.cychp"
#"data-cels/20NA13330_80NA07126_A05_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.CEL",
#"data-cels/20NA13330_80NA07126_A05_MosaicismStudy_51pt0_CytoScanHD_MS_20110412.CEL",
#"data-cels/20NA13330_80NA07126_A05_MosaicismStudy_51pt0_CytoScanHD_WW_20110411.CEL",
#"data-cels/20NA13330_80NA07126_B05_MosaicismStudy_50pt0_CytoScanHD_MS_20110412.CEL",
#"data-cels/20NA13330_80NA07126_B05_MosaicismStudy_50pt0_CytoScanHD_WW_20110411.CEL",
#"data-cels/20NA13330_80NA07126_B05_MosaicismStudy_50pt5_CytoScanHD_CC_20110411.CEL",
#"data-cels/20NA13330_80NA07126_C05_MosaicismStudy_50pt5_CytoScanHD_MS_20110412.CEL",
#"data-cels/20NA13330_80NA07126_C05_MosaicismStudy_50pt5_CytoScanHD_WW_20110411.CEL",
#"data-cels/20NA13330_80NA07126_C05_MosaicismStudy_51pt0_CytoScanHD_CC_20110411.CEL",
#"data-cels/30NA13330_70NA07126_A04_MosaicismStudy_50pt0_CytoScanHD_CC_20110411.CEL",
#"data-cels/30NA13330_70NA07126_A04_MosaicismStudy_51pt0_CytoScanHD_MS_20110412.CEL",
#"data-cels/30NA13330_70NA07126_A04_MosaicismStudy_51pt0_CytoScanHD_WW_20110411.CEL",
#"data-cels/30NA13330_70NA07126_B04_MosaicismStudy_50pt0_CytoScanHD_MS_20110412.CEL",
#"data-cels/30NA13330_70NA07126_B04_MosaicismStudy_50pt0_CytoScanHD_WW_20110411.CEL",
#"data-cels/30NA13330_70NA07126_B04_MosaicismStudy_50pt5_CytoScanHD_CC_20110411.CEL",
#"data-cels/30NA13330_70NA07126_C04_MosaicismStudy_50pt5_CytoScanHD_MS_20110412.CEL",
#"data-cels/30NA13330_70NA07126_C04_MosaicismStudy_50pt5_CytoScanHD_WW_20110411.CEL",
#"data-cels/30NA13330_70NA07126_C04_MosaicismStudy_51pt0_CytoScanHD_CC_20110411.CEL"
)

print(filenames)

##########

for (fn in filenames) {
  cat("### filename: ",fn,"\n")

  debug_current_filename <<- fn
  basename <-  tail(strsplit(fn,"/")[[1]],n=1)

  temp <- Read.CYCHP(fn)

  temp <- data.frame(temp$DataGroup$ProbeSets$Datasets$CopyNumber$DataColumns)
  temp <- temp[!is.na(temp$Log2Ratio),]

  All.seg <- NULL
  for (chr in c(1:22,24,25)){
    cur.chr <- temp[temp[,2]==chr,]

    debug_current_chr <<- formatC(chr,format="d",width=2,flag="0")
    cat("### debug_current_chr",debug_current_chr,"\n")
    debug_current_basename <<- paste(out_dir,"/",basename,"-",debug_current_chr,sep="")
    cat("### debug_current_basename",debug_current_basename,"\n")

    cur.seg <- mosaic.segmentation.algorithm(cur.chr[,4], gains.boundaries, losses.boundaries,bandwidth=6000)

    cur.seg <- cbind(cur.seg,0,0)

    for (i in 1:nrow(cur.seg)) {
      cur.seg[i,5] <-  cur.seg[i,2] -  cur.seg[i,1] + 1
      cur.seg[i,6] <-trunc(mean(diff(cur.chr[cur.seg[i,1]:cur.seg[i,2],3])))

      cur.seg[i,1] <- cur.chr[ cur.seg[i,1],3]
      cur.seg[i,2] <- cur.chr[ cur.seg[i,2],3]
    }
    All.seg <- c(All.seg,list(cbind(chr,cur.seg)))
  }

  All.seg <- do.call("rbind",All.seg)

  colnames(All.seg) <- c("Chromosome","StartPosition","StopPosition","Mixture","Confidence","MarkerCount","MeanMarkerDistance")

  write.table(All.seg,file=paste(out_dir,"/",basename,".protoMosaicSeg",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
}
