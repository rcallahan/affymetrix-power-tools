kitA1 <- read.table("kitA_cond1_signal_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitO2 <- read.table("kitO_cond2_signal_size.txt",header=T,sep="\t",stringsAsFactors=F)
new2 <- read.table("new_cond2_signal_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitA1sub <- kitA1[which(kitA1$probeset_id %in% gs_ps_sub),]
kitO2sub <- kitO2[which(kitO2$probeset_id %in% gs_ps_sub),]
new2sub <- new2[which(new2$probeset_id %in% gs_ps_sub),]
kitA1subT <- t(kitA1sub)
kitO2subT <- t(kitO2sub)
new2subT <- t(new2sub)
new2subT <- as.data.frame(t(apply(new2subT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitO2subT <- as.data.frame(t(apply(kitO2subT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitA1subT <- as.data.frame(t(apply(kitA1subT[-1, ],1,as.numeric)),stringsAsFactors=F)
names(new2subT) <- new2sub$probeset_id
names(kitO2subT) <- kitO2sub$probeset_id
names(kitA1subT) <- kitA1sub$probeset_id
kitA1subT <- data.frame(matrix("A",nrow(kitA1subT),1),kitA1subT,stringsAsFactors=F)
names(kitA1subT)[1] <- "Kit"
kitO2subT <- data.frame(matrix("O",nrow(kitO2subT),1),kitO2subT,stringsAsFactors=F)
names(kitO2subT)[1] <- "Kit"
new2subT <- data.frame(matrix("O",nrow(new2subT),1),new2subT,stringsAsFactors=F)
names(new2subT)[1] <- "Kit"

allData <- rbind(kitA1subT,kitO2subT,new2subT)
mySvm <- svm(as.factor(Kit) ~ .,data=allData)

plot(myPCA300$x[1:24,1],myPCA300$x[1:24,2],xlim=c(-190,-178),ylim=c(-30,15),col=1,xlab="PC1",ylab="PC2")
points(myPCA300$x[25:66,1],myPCA300$x[25:66,2],col=2)

kitA1CEUsub <- kitA1CEU[which(kitA1CEU$probeset_id %in% top300),]
kitO2CEUsub <- kitO2CEU[which(kitO2CEU$probeset_id %in% top300),]

gs_chr_pos_300 <- paste(annot[which(annot$probeset_id %in% top300),]$chr,"_",annot[which(annot$probeset_id %in% top300),]$snp_pos,sep="")
ceu_in_gs_ind <- which(ceu_chr_pos %in% gs_chr_pos_300 )
ceuTop300 <- ceuAnnot[ceu_in_gs_ind,]$probeset_id

kitA1CEU300T <- t(kitA1CEU300)
kitO2CEU300T <- t(kitO2CEU300)

kitA1CEU300T <- as.data.frame(t(apply(kitA1CEU300T[-1, ],1,as.numeric)),stringsAsFactors=F)
kitO2CEU300T <- as.data.frame(t(apply(kitO2CEU300T[-1, ],1,as.numeric)),stringsAsFactors=F)
names(kitO2CEU300T) <- kitO2CEU300$probeset_id
names(kitA1CEU300T) <- kitA1CEU300$probeset_id

plot(ceuPC1[1:24],ceuPC2[1:24],xlim=c(-182,-170),ylim=c(-30,15),col=1,xlab="PC1",ylab="PC2")
points(ceuPC1[25:48],ceuPC2[25:48],col=2)


############################################################################################
dups <- which(duplicated(ceuAnnot$probeset_id))
ceuAnnot <- ceuAnnot[-dups,]

gs_chr_pos <- paste(annot$chr,"_",annot$snp_pos,sep="")
ceu_chr_pos <- paste(ceuAnnot$chr_b36,"_",ceuAnnot$pos_b36,sep="")

gs_ps_57K <- kitA1sub$probeset_id
annot <- data.frame(gs_chr_pos,annot,stringsAsFactors=F)
annot57K <- annot[which(annot$probeset_id %in% gs_ps_57K),]
ceuAnnot <- data.frame(ceu_chr_pos,ceuAnnot,stringsAsFactors=F)
ceuGSAnnot <- merge(annot57K,ceuAnnot,by.x="gs_chr_pos",by.y="ceu_chr_pos")
ceu_to_gs_map <- data.frame(GS_probeset_id = ceuGSAnnot$probeset_id.x,CEU_probeset_id = ceuGSAnnot$probeset_id.y,stringsAsFactors=F)

kitA1CEUsub <- kitA1CEU[which(kitA1CEU$probeset_id %in% ceu_asi_gs_map[,1]),]
kitO2CEUsub <- kitO2CEU[which(kitO2CEU$probeset_id %in% ceu_asi_gs_map[,1]),]
kitA1CEUsubT <- t(kitA1CEUsub)
kitO2CEUsubT <- t(kitO2CEUsub)
kitA1CEUsubT <- as.data.frame(t(apply(kitA1CEUsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitO2CEUsubT <- as.data.frame(t(apply(kitO2CEUsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
names(kitO2CEUsubT) <- kitO2CEUsub$probeset_id
names(kitA1CEUsubT) <- kitA1CEUsub$probeset_id

allData48 <- allData[1:48,]

firstPCA <- list()
top300list <- list()
secondPCA300 <- list()
projectionPC1 <- list()
projectionPC2 <- list()
for (i in 1:5){
	jpeg(paste("PCA_validation",i,".jpg",sep=""),width=1200,height=400)
	par(mfrow=c(1,3))
	rand24 <- sort(sample(1:24,19))
	rand48 <- sort(sample(25:48,19))

	randSample <- c(rand24,rand48)
	allCEUDataRand <- allCEUData[randSample,]

	myRandPCA <- prcomp(allCEUDataRand,center=F)
	firstPCA[[i]] <- myRandPCA
	plot(myRandPCA$x[1:19,1],myRandPCA$x[1:19,2],xlim=c(-2432,-2420),ylim=c(-150,100),col=1,xlab="PC1",ylab="PC2",main="CEU 57000+")
	points(myRandPCA$x[20:38,1],myRandPCA$x[20:38,2],col=2)

	myPC <- data.frame(row.names(myRandPCA$rotation),abs(myRandPCA$rotation[,2]),stringsAsFactors=F)
	myPC <- myPC[order(myPC[,2],decreasing=T),]
	testTop300 <- myPC[1:300,]
	top300list[[i]] <- testTop300
	
	allCEUDataRand300 <- allCEUDataRand[,which(names(allCEUDataRand) %in% testTop300[,1])]
	myRandPCA300 <- prcomp(allCEUDataRand300,center=F)
	secondPCA300[[i]] <- myRandPCA300

	plot(myRandPCA300$x[1:19,1],myRandPCA300$x[1:19,2],xlim=c(-190,-170),ylim=c(-30,20),col=1,xlab="PC1",ylab="PC2",main="CEU Top 300")
	points(myRandPCA300$x[20:38,1],myRandPCA300$x[20:38,2],col=2)

	gs_ps_sub_300 <- ceu_to_gs_map$GS_probeset_id[which(ceu_to_gs_map$CEU_probeset_id %in% testTop300[,1])]
	allData48Rand300 <- allData48[-randSample,which(names(allData48) %in% gs_ps_sub_300)]

	PC1 <- as.matrix(allData48Rand300) %*% as.matrix(myRandPCA300$rotation[,1])
	PC2 <- as.matrix(allData48Rand300) %*% as.matrix(myRandPCA300$rotation[,2])
	projectionPC1[[i]] <- PC1
	projectionPC2[[i]] <- PC2

	plot(PC1[1:5],PC2[1:5],xlim=c(-190,-170),ylim=c(-30,20),col=1,xlab="PC1",ylab="PC2",main="CEU Projection onto GS")
	points(PC1[6:10],PC2[6:10],col=2)
	dev.off()
}

map_all_prods <- gsCeuAsiUCSFHanAnnot[,c(1,3,21)]
dups <- which(duplicated(map_all_prods[,1]))
map_all_prods <- map_all_prods[-dups,]
allCEUData <- allCEUData[,(which(names(allCEUData) %in% map_all_prods$CEU_probeset_id))]
myPCA <- prcomp(allCEUData)
myPC <- data.frame(row.names(myPCA$rotation),abs(myPCA$rotation[,1]),stringsAsFactors=F)
myPC <- myPC[order(myPC[,2],decreasing=T),]
testTop300 <- myPC[1:300,]
allCEUData300 <- allCEUData[,which(names(allCEUData) %in% testTop300[,1])]
myPCA300 <- prcomp(allCEUData300)
means <- apply(allCEUData300,2,mean)
psPC1means <- data.frame(CEU_probeset_id = names(myPCA300$rotation[,1]),PC1 = myPCA300$rotation[,1],means,stringsAsFactors=F)
psPC1means <- merge(map_all_prods,psPC1means,by="CEU_probeset_id")

simCEU <- c()
for (i in 1:48){
	myVector <- as.vector(allCEUData[i,])
	cat(i,"\n")
	for (j in 1:100){
		thisSampleInd <- sample(1:ncol(allCEUData),300)
		sampleVector <- myVector[thisSampleInd]
		names(sampleVector) <- c(1:300)
		simCEU <- rbind(simCEU,sampleVector)
	}
}

simCEU2 <- c()
for (i in 1:4800){
	if ((i %% 100) == 0){
		cat(i,"\n")
	}
	rowVec <- c()
	for (j in 1:300){
		thisSampleInd <- sample(1:48,1)
		mySample <- allCEUData300[thisSampleInd,j]
		rowVec <- cbind(rowVec,mySample)
	}
	simCEU2 <- rbind(simCEU2,rowVec)
}

jpeg("2_simulations.jpg",width=1200,height=600)
par(mfrow=c(1,2))
simCEUC <- t(apply(simCEU,1,function(x) x-psPC1means$means))
simCEUC <- as.data.frame(simCEUC,stringsAsFactors=F)
names(simCEUC) <- names(simCEU)
PC1 <- as.matrix(simCEUC[,order(names(simCEUC))]) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),1])
PC2 <- as.matrix(simCEUC[,order(names(simCEUC))]) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),2])
plot(PC1[1:2400],PC2[1:2400],xlim=c(-30,30),ylim=c(-10,10),col=1,xlab="PC1",ylab="PC2",main="4800 Simulated 1")
points(PC1[2401:4800],PC2[2401:4800],col=2)
points(myPCA300$x[,1],myPCA300$x[,2])
PC1nul1 <- PC1

simCEU2C <- t(apply(simCEU2,1,function(x) x-psPC1means$means))
simCEUC2 <- as.data.frame(simCEU2C,stringsAsFactors=F)
names(simCEU2C) <- names(simCEU2)
PC1 <- as.matrix(simCEU2C) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),1])
PC2 <- as.matrix(simCEU2C) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),2])
plot(PC1[1:4800],PC2[1:4800],xlim=c(-30,30),ylim=c(-10,10),col=1,xlab="PC1",ylab="PC2",main="4800 Simulated 2")
points(myPCA300$x[,1],myPCA300$x[,2])
PC1nul2 <- PC1

dev.off()

newCond1GS <- read.table("new_cond1_size.txt",header=T,sep="\t",stringsAsFactors=F)
newCond2GS <- read.table("new_cond2_size.txt",header=T,sep="\t",stringsAsFactors=F)
newCond3GS <- read.table("new_cond3_size.txt",header=T,sep="\t",stringsAsFactors=F)
newCond4GS <- read.table("new_cond4_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitOCond1GS <- read.table("kitO_cond1_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitOCond2GS <- read.table("kitO_cond2_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitACond1GS <- read.table("kitA_cond1_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitACond2GS <- read.table("kitA_cond2_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitACond3GS <- read.table("kitA_cond3_size.txt",header=T,sep="\t",stringsAsFactors=F)
kitACond4GS <- read.table("kitA_cond4_size.txt",header=T,sep="\t",stringsAsFactors=F)

psPC1means <- psPc1means[order(psPc1means$GS_probeset_id),]

newCond1GSsub <- newCond1GS[which(newCond1GS$probeset_id %in% psPC1means[,2]),]
newCond2GSsub <- newCond2GS[which(newCond2GS$probeset_id %in% psPC1means[,2]),]
newCond3GSsub <- newCond3GS[which(newCond3GS$probeset_id %in% psPC1means[,2]),]
newCond4GSsub <- newCond4GS[which(newCond4GS$probeset_id %in% psPC1means[,2]),]
kitOCond1GSsub <- kitOCond1GS[which(kitOCond1GS$probeset_id %in% psPC1means[,2]),]
kitOCond2GSsub <- kitOCond2GS[which(kitOCond2GS$probeset_id %in% psPC1means[,2]),]
kitACond1GSsub <- kitACond1GS[which(kitACond1GS$probeset_id %in% psPC1means[,2]),]
kitACond2GSsub <- kitACond2GS[which(kitACond2GS$probeset_id %in% psPC1means[,2]),]
kitACond3GSsub <- kitACond3GS[which(kitACond3GS$probeset_id %in% psPC1means[,2]),]
kitACond4GSsub <- kitACond4GS[which(kitACond4GS$probeset_id %in% psPC1means[,2]),]

newCond1GSsubT <- t(newCond1GSsub)
newCond2GSsubT <- t(newCond2GSsub)
newCond3GSsubT <- t(newCond3GSsub) 
newCond4GSsubT <- t(newCond4GSsub) 
kitOCond1GSsubT <- t(kitOCond1GSsub) 
kitOCond2GSsubT <- t(kitOCond2GSsub) 
kitACond1GSsubT <- t(kitACond1GSsub) 
kitACond2GSsubT <- t(kitACond2GSsub) 
kitACond3GSsubT <- t(kitACond3GSsub) 
kitACond4GSsubT <- t(kitACond4GSsub) 

newCond1GSsubT <- as.data.frame(t(apply(newCond1GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
newCond2GSsubT <- as.data.frame(t(apply(newCond2GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
newCond3GSsubT <- as.data.frame(t(apply(newCond3GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
newCond4GSsubT <- as.data.frame(t(apply(newCond4GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitOCond1GSsubT <- as.data.frame(t(apply(kitOCond1GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitOCond2GSsubT <- as.data.frame(t(apply(kitOCond2GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitACond1GSsubT <- as.data.frame(t(apply(kitACond1GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitACond2GSsubT <- as.data.frame(t(apply(kitACond2GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitACond3GSsubT <- as.data.frame(t(apply(kitACond3GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)
kitACond4GSsubT <- as.data.frame(t(apply(kitACond4GSsubT[-1, ],1,as.numeric)),stringsAsFactors=F)

names(newCond1GSsubT) <- newCond1GSsub$probeset_id
names(newCond2GSsubT) <- newCond2GSsub$probeset_id
names(newCond3GSsubT) <- newCond3GSsub$probeset_id
names(newCond4GSsubT) <- newCond4GSsub$probeset_id
names(kitOCond1GSsubT) <- kitOCond1GSsub$probeset_id
names(kitOCond2GSsubT) <- kitOCond2GSsub$probeset_id
names(kitACond1GSsubT) <- kitACond1GSsub$probeset_id
names(kitACond2GSsubT) <- kitACond2GSsub$probeset_id
names(kitACond3GSsubT) <- kitACond3GSsub$probeset_id
names(kitACond4GSsubT) <- kitACond4GSsub$probeset_id

newCond1GSsubT <-  newCond1GSsubT[,order(names(newCond1GSsubT))]  
newCond2GSsubT <-  newCond2GSsubT[,order(names(newCond2GSsubT))]  
newCond3GSsubT <-  newCond3GSsubT[,order(names(newCond3GSsubT))]  
newCond4GSsubT <-  newCond4GSsubT[,order(names(newCond4GSsubT))]  
kitOCond1GSsubT <- kitOCond1GSsubT[,order(names(kitOCond1GSsubT))]
kitOCond2GSsubT <- kitOCond2GSsubT[,order(names(kitOCond2GSsubT))]
kitACond1GSsubT <- kitACond1GSsubT[,order(names(kitACond1GSsubT))]
kitACond2GSsubT <- kitACond2GSsubT[,order(names(kitACond2GSsubT))]
kitACond3GSsubT <- kitACond3GSsubT[,order(names(kitACond3GSsubT))]
kitACond4GSsubT <- kitACond4GSsubT[,order(names(kitACond4GSsubT))]

newCond1GSsubTC <-  t(apply(newCond1GSsubT,1,function(x) x- psPC1means$means))
newCond2GSsubTC <-  t(apply(newCond2GSsubT,1,function(x) x- psPC1means$means)) 
newCond3GSsubTC <-  t(apply(newCond3GSsubT,1,function(x) x- psPC1means$means))
newCond4GSsubTC <-  t(apply(newCond4GSsubT,1,function(x) x- psPC1means$means)) 
kitOCond1GSsubTC <- t(apply(kitOCond1GSsubT,1,function(x) x- psPC1means$means))
kitOCond2GSsubTC <- t(apply(kitOCond2GSsubT,1,function(x) x- psPC1means$means))
kitACond1GSsubTC <- t(apply(kitACond1GSsubT,1,function(x) x- psPC1means$means))
kitACond2GSsubTC <- t(apply(kitACond2GSsubT,1,function(x) x- psPC1means$means))
kitACond3GSsubTC <- t(apply(kitACond3GSsubT,1,function(x) x- psPC1means$means))
kitACond4GSsubTC <- t(apply(kitACond4GSsubT,1,function(x) x- psPC1means$means))

newCond1GSsubTproj <- as.matrix(newCond1GSsubTC) %*% as.matrix(psPC1means$PC1)
newCond2GSsubTproj <- as.matrix(newCond2GSsubTC) %*% as.matrix(psPC1means$PC1)
newCond3GSsubTproj <- as.matrix(newCond3GSsubTC) %*% as.matrix(psPC1means$PC1)
newCond4GSsubTproj <- as.matrix(newCond4GSsubTC) %*% as.matrix(psPC1means$PC1)
kitOCond1GSsubTproj <-as.matrix(kitOCond1GSsubTC)%*% as.matrix(psPC1means$PC1)
kitOCond2GSsubTproj <-as.matrix(kitOCond2GSsubTC)%*% as.matrix(psPC1means$PC1)
kitACond1GSsubTproj <-as.matrix(kitACond1GSsubTC)%*% as.matrix(psPC1means$PC1)
kitACond2GSsubTproj <-as.matrix(kitACond2GSsubTC)%*% as.matrix(psPC1means$PC1)
kitACond3GSsubTproj <-as.matrix(kitACond3GSsubTC)%*% as.matrix(psPC1means$PC1)
kitACond4GSsubTproj <-as.matrix(kitACond4GSsubTC)%*% as.matrix(psPC1means$PC1)


h1 <- density(newCond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h2 <- density(newCond2GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h3 <- density(newCond3GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h4 <- density(newCond4GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h5 <- density(kitOCond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h6 <- density(kitOCond2GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h7 <- density(kitACond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h8 <- density(kitACond2GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h9 <- density(kitACond3GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h10 <- density(kitACond4GSsubTproj)#,breaks=seq(-15, 15, 0.1))

plot(h1,xlim=c(-15,15),col=1,main="Projection on PC1",xlab="PC1")
lines(h2,col=2)
lines(h3,col=3)
lines(h4,col=4)
lines(h5,col=5)
lines(h6,col=6)
lines(h7,col=7)
#lines(h8,col=8)
#lines(h9,col=9)
#lines(h10,col=10)
h11 <- density(PC1nul1)
h12 <- density(PC1nul2)
lines(h11,col=11)
lines(h12,col=12)
legend("topleft", legend=c("newCond1","newCond2","newCond3","newCond4","kitOCond","kitOCond2","kitACond1","Null-1","Null-2"),pch=15,col=c(1,2,3,4,5,6,7,11,12))

,"kitACond2","kitACond3","kitACond4"
,8,9,10

addZeros <- list()
addZeros[[1]] <-  newCond1GSsubT   
addZeros[[2]] <-  newCond2GSsubT   
addZeros[[3]] <-  newCond3GSsubT   
addZeros[[4]] <-  newCond4GSsubT   
addZeros[[5]]  <- kitOCond1GSsubT  
addZeros[[6]]  <- kitOCond2GSsubT 
addZeros[[7]]  <- kitACond1GSsubT  
#addZeros[[8]]  <- kitACond2GSsubT  
#addZeros[[9]]  <- kitACond3GSsubT  
#addZeros[[10]] <- kitACond4GSsubT
myProj <- list()  
for (i in 1:7){
	myRows <- nrow(addZeros[[i]])
	for (j in 1:myRows){
		mySample <- sample(1:300,30)
		addZeros[[i]][j,mySample] <- 0
	}
	print(summary(names(addZeros[[i]])==psPC1means$GS_probeset_id))
	addZeros[[i]] <-  t(apply(addZeros[[i]],1,function(x) x- psPC1means$means))
	myProj[[i]] <- as.matrix(addZeros[[i]]) %*% as.matrix(psPC1means$PC1)
	h <- density(myProj[[i]])
	if (i == 1){
		plot(h,xlim=c(-20,20),ylim=c(0,0.6),col=i,main="Projection on PC with Zeros",xlab="PC1")
	} else {
		lines(h,col=i)
	}	
}
h11 <- density(PC1nul1)
h12 <- density(PC1nul2)
lines(h11,col=11)
lines(h12,col=12)
legend("topleft", legend=c("newCond1","newCond2","newCond3","newCond4","kitOCond","kitOCond2","kitACond1","Null-1","Null-2"),pch=15,col=c(1,2,3,4,5,6,7,11,12))

,"kitACond2","kitACond3","kitACond4"
,8,9,10
addZeros <- list()
addZeros[[1]] <-  newCond1GSsubT   
addZeros[[2]] <-  newCond2GSsubT   
addZeros[[3]] <-  newCond3GSsubT   
addZeros[[4]] <-  newCond4GSsubT   
addZeros[[5]]  <- kitOCond1GSsubT  
addZeros[[6]]  <- kitOCond2GSsubT  
addZeros[[7]]  <- kitACond1GSsubT  
#addZeros[[8]]  <- kitACond2GSsubT  
#addZeros[[9]]  <- kitACond3GSsubT  
#addZeros[[10]] <- kitACond4GSsubT 
myProj <- list()  
for (i in 1:7){
	myRows <- nrow(addZeros[[i]])
	for (j in 1:myRows){
		mySample <- rnorm(300)
		addZeros[[i]][j,] <- addZeros[[i]][j,] + mySample
	}
	addZeros[[i]] <-  t(apply(addZeros[[i]],1,function(x) x- psPC1means$means))
	myProj[[i]] <- as.matrix(addZeros[[i]]) %*% as.matrix(psPC1means$PC1)
	h <- density(myProj[[i]])
	if (i == 1){
		plot(h,xlim=c(-15,15),ylim=c(0,0.6),col=i,main="Projection on PC with Random Noise",xlab="PC1")
	} else {
		lines(h,col=i)
	}	
}
h11 <- density(PC1nul1)
h12 <- density(PC1nul2)
lines(h11,col=11)
lines(h12,col=12)
legend("topleft", legend=c("newCond1","newCond2","newCond3","newCond4","kitOCond","kitOCond2","kitACond1","Null-1","Null-2"),pch=15,col=c(1,2,3,4,5,6,7,11,12))

"kitACond2","kitACond3","kitACond4",
,8,9,10

firstPCA <- list()
top300list <- list()
secondPCA300 <- list()
projectionPC2KA <- list()
projectionPC2KO <- list()
for (i in 1:5){
	rand24 <- sort(sample(1:24,19))
	rand48 <- sort(sample(25:48,19))

	randSample <- c(rand24,rand48)
	allCEUDataRand <- allCEUData[randSample,]

	myRandPCA <- prcomp(allCEUDataRand)
	firstPCA[[i]] <- myRandPCA
	
	myPC <- data.frame(row.names(myRandPCA$rotation),abs(myRandPCA$rotation[,1]),stringsAsFactors=F)
	myPC <- myPC[order(myPC[,2],decreasing=T),]
	testTop300 <- myPC[1:300,]
	top300list[[i]] <- testTop300
	
	allCEUDataRand300 <- allCEUDataRand[,which(names(allCEUDataRand) %in% testTop300[,1])]
	means <- apply(allCEUDataRand300,2,mean)
	myRandPCA300 <- prcomp(allCEUDataRand300)
	secondPCA300[[i]] <- myRandPCA300


	allCEUData10Rand300 <- allCEUData[-randSample,which(names(allCEUData) %in% testTop300[,1])]
	allCEUData10Rand300 <- allCEUData10Rand300[,order(names(allCEUData10Rand300))]
	myRotation <- myRandPCA300$rotation[order(row.names(myRandPCA300$rotation)),1]
	allCEUData10Rand300 <-  t(apply(allCEUData10Rand300,1,function(x) x- means))
	PC2 <- as.matrix(allCEUData10Rand300) %*% as.matrix(myRotation)
	projectionPC2KA[[i]] <- PC2[1:5]
	projectionPC2KO[[i]] <- PC2[6:10]
	h1 <- density(projectionPC2KA[[i]])
	h2 <- density(projectionPC2KO[[i]])
	if (i == 1){
		plot(h1,xlim=c(-20,20),ylim=c(0,1),col=1,main="CEU Projection on PC1",xlab="PC1")
		lines(h2,col=2)
	
	} else {
		myColor <- sample(1:100,1)
		lines(h1,col=1)
		lines(h2,col=2)
	}
	
}

bestPC2 <- myPCA300$rotation[,2]
bestPC2 <- data.frame(names(bestPC2),bestPC2,stringsAsFactors=F)
names(bestPC2) <- c("CEU_probeset_id","Value")
bestPC2 <- merge(bestPC2,ceu_to_gs_map,by.x="CEU_probeset_id",by.y="CEU_probeset_id")
bestPC2 <- bestPC2[,c(1,3,2)]
bestPC2 <- bestPC2[order(bestPC2$GS_probeset_id),]

newCond1GSsub <- newCond1GS[which(newCond1GS$probeset_id %in% bestPC2$GS_probeset_id),]
newCond2GSsub <- newCond2GS[which(newCond2GS$probeset_id %in% bestPC2$GS_probeset_id),]
newCond3GSsub <- newCond3GS[which(newCond3GS$probeset_id %in% bestPC2$GS_probeset_id),]
newCond4GSsub <- newCond4GS[which(newCond4GS$probeset_id %in% bestPC2$GS_probeset_id),]
kitOCond1GSsub <- kitOCond1GS[which(kitOCond1GS$probeset_id %in% bestPC2$GS_probeset_id),]
kitOCond2GSsub <- kitOCond2GS[which(kitOCond2GS$probeset_id %in% bestPC2$GS_probeset_id),]
kitACond1GSsub <- kitACond1GS[which(kitACond1GS$probeset_id %in% bestPC2$GS_probeset_id),]

newCond1GSsubTproj <- as.matrix(newCond1GSsubT) %*% as.matrix(bestPC2[,3])
newCond2GSsubTproj <- as.matrix(newCond2GSsubT) %*% as.matrix(bestPC2[,3])
newCond3GSsubTproj <- as.matrix(newCond3GSsubT) %*% as.matrix(bestPC2[,3])
newCond4GSsubTproj <- as.matrix(newCond4GSsubT) %*% as.matrix(bestPC2[,3])
kitOCond1GSsubTproj <-as.matrix(kitOCond1GSsubT)%*% as.matrix(bestPC2[,3])
kitOCond2GSsubTproj <-as.matrix(kitOCond2GSsubT)%*% as.matrix(bestPC2[,3])
kitACond1GSsubTproj <-as.matrix(kitACond1GSsubT)%*% as.matrix(bestPC2[,3])

h1 <- density(newCond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h2 <- density(newCond2GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h3 <- density(newCond3GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h4 <- density(newCond4GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h5 <- density(kitOCond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h6 <- density(kitOCond2GSsubTproj)#,breaks=seq(-15, 15, 0.1))
h7 <- density(kitACond1GSsubTproj)#,breaks=seq(-15, 15, 0.1))

plot(h1,xlim=c(-15,15),col=1,main="GS Projection on CEU PC2",xlab="PC2")
lines(h2,col=2)
lines(h3,col=3)
lines(h4,col=4)
lines(h5,col=5)
lines(h6,col=6)
lines(h7,col=7)
lines(h11,col=11)
lines(h12,col=12)
legend("topleft", legend=c("newCond1","newCond2","newCond3","newCond4","kitOCond","kitOCond2","kitACond1","Null-1","Null-2"),pch=15,col=c(1,2,3,4,5,6,7,11,12))

classify <- list()
classify[[1]] <-  newCond1GSsubTproj   
classify[[2]] <-  newCond2GSsubTproj   
classify[[3]] <-  newCond3GSsubTproj   
classify[[4]] <-  newCond4GSsubTproj   
classify[[5]]  <- kitOCond1GSsubTproj  
classify[[6]]  <- kitOCond2GSsubTproj  
classify[[7]]  <- kitACond1GSsubTproj

classMatrix <- matrix(0,3,3)
#classify <- myProj
thres <- 4.5
for (i in 1:7){
	myClassData <- classify[[i]]
	if(i != 7){
		for (j in 1:length(myClassData)){
			if (myClassData[j] < -thres){
				classMatrix[2,1] <- classMatrix[2,1] + 1
			}
			if (myClassData[j] > thres){
				classMatrix[2,2] <- classMatrix[2,2] + 1
			}
			if (myClassData[j] >= -thres & myClassData[j] <= thres){
				classMatrix[2,3] <- classMatrix[2,3] + 1
			}
		}
	} else {
		for (j in 1:length(myClassData)){
			if (myClassData[j] < -thres){
				classMatrix[1,1] <- classMatrix[1,1] + 1
			}
			if (myClassData[j] > thres){
				classMatrix[1,2] <- classMatrix[1,2] + 1
			}
			if (myClassData[j] >= -thres & myClassData[j] <= thres){
				classMatrix[1,3] <- classMatrix[1,3] + 1
			}
		}
	}
}
classMatrix


for (i in 1:5){
	h1 <- density(projectionPC2KA[[i]])
	h2 <- density(projectionPC2KO[[i]])
	if (i == 1){
		plot(h1,xlim=c(-20,20),ylim=c(0,1),col=1,main="Cross Validation: CEU Projection on PC1",xlab="PC2")
		lines(h2,col=2)
	
	} else {
		myColor <- sample(1:100,1)
		lines(h1,col=1)
		lines(h2,col=2)
	}
	
}
lines(h11,col=11)
lines(h12,col=12)
legend("topright", legend=c("KitA","KitO","Null-1","Null-2"),pch=15,col=c(1,2,11,12))

