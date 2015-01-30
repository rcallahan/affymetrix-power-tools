
set.seed(1)
tail.type <- list(
  "less" = 0,
  "greater" = 1,
  "two.sided" = 2,
)
n <- 100

# Test data for Wilcoxon signed rank test
outFile <- "signrank_test.tab"
sampleSize <- sample(1:70,n,replace=TRUE)
for(i in 1:n) {
  dataset <- round(rnorm(sampleSize[i]),3)
  thisTail <- sample(unlist(tail.type),1)
  pval <- wilcox.test(dataset,alternative=names(thisTail))$p.value
  outRow <- matrix(c(pval,thisTail,sampleSize[i],dataset),nrow=1)
  if(i==1) {
    write.table(outRow,file=outFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  } else {
    write.table(outRow,file=outFile,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }
}

# Test data for Wilcoxon rank sum test
outFile <- "ranksum_test.tab"
sampleSize1 <- sample(1:70,n,replace=TRUE)
sampleSize2 <- sample(1:70,n,replace=TRUE)
for(i in 1:n) {
  dataset1 <- round(rnorm(sampleSize1[i]),3)
  dataset2 <- round(rnorm(sampleSize2[i]),3)
  thisTail <- sample(unlist(tail.type),1)
  pval <- wilcox.test(dataset1,dataset2,alternative=names(thisTail))$p.value
  outRow <- matrix(c(pval,thisTail,sampleSize1[i],sampleSize2[i],dataset1,dataset2),nrow=1)
  if(i==1) {
    write.table(outRow,file=outFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  } else {
    write.table(outRow,file=outFile,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }
}
