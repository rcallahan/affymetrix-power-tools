allFiles <- dir()[grep("log2",dir())]

fileHeaders <- unique(sapply(strsplit(allFiles, "\\."), "[", 1))

kitaInds <- grep("KitA",fileHeaders)
kitoInds <- grep("KitO",fileHeaders)

fileHeaders <- fileHeaders[c(kitaInds,kitoInds)]

allSizes <- c()
for (myHeader in fileHeaders){
	myFile <- paste(myHeader,".log2",sep="")
	cat("Working on",myFile,"\n")
	myData <- read.table(myFile,header=T,sep="\t",stringsAsFactors=F)
	myData <- myData[!is.na(myData$log2.A.B.),]
	myMed <- tapply(myData$log2.A.B.,myData$probeset_id,median)
	allSizes <- rbind(allSizes,myMed)

}

allSizes <- data.frame(allSizes, row.names=fileHeaders, stringsAsFactors=F)
write.table(allSizes,"KitO_CEU_signal_size.txt",row.names=F,quote=F,sep="\t")

map_all_prods <- read.table("map_all_prods.txt",header=T,sep="\t",stringsAsFactors=F)
psInd <- which(names(allSizes) %in% map_all_prods$GS_probeset_id)
allGSData <- allSizes[,psInd]
kitoCond2Ind <- grep("KitO_Cond2",row.names(allGSData))
allGSData48 <- allGSData[c(1:24,kitoCond2Ind),]

ceuInd <- which(names(allSizes) %in% map_all_prods$CEU_probeset_id)
allCEUData <- allSizes[,ceuInd]

myPCA <- prcomp(allCEUData)
myPC <- data.frame(row.names(myPCA$rotation),abs(myPCA$rotation[,1]),stringsAsFactors=F)
myPC <- myPC[order(myPC[,2],decreasing=T),]
testTop300 <- myPC[1:300,]
allCEUData300 <- allCEUData[,which(names(allCEUData) %in% testTop300[,1])]
myPCA300 <- prcomp(allCEUData300)
means <- apply(allCEUData300,2,mean)
psPC1means <- data.frame(CEU_probeset_id = names(myPCA300$rotation[,1]),PC1 = myPCA300$rotation[,1],means,stringsAsFactors=F)
psPC1means <- merge(map_all_prods,psPC1means,by="CEU_probeset_id")

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
		plot(h1,xlim=c(-15,15),ylim=c(0,1),col=1,main="CEU Projection on PC1",xlab="PC1")
		lines(h2,col=2)
	
	} else {
		myColor <- sample(1:100,1)
		lines(h1,col=1)
		lines(h2,col=2)
	}
	
}

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
plot(PC1[1:2400],PC2[1:2400],xlim=c(-15,15),ylim=c(-10,10),col=1,xlab="PC1",ylab="PC2",main="4800 Simulated 1")
points(PC1[2401:4800],PC2[2401:4800],col=2)
points(myPCA300$x[,1],myPCA300$x[,2])
PC1nul1 <- PC1

simCEU2C <- t(apply(simCEU2,1,function(x) x-psPC1means$means))
simCEUC2 <- as.data.frame(simCEU2C,stringsAsFactors=F)
names(simCEU2C) <- names(simCEU2)
PC1 <- as.matrix(simCEU2C) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),1])
PC2 <- as.matrix(simCEU2C) %*% as.matrix(myPCA300$rotation[order(row.names(myPCA300$rotation)),2])
plot(PC1[1:4800],PC2[1:4800],xlim=c(-15,15),ylim=c(-10,10),col=1,xlab="PC1",ylab="PC2",main="4800 Simulated 2")
points(myPCA300$x[,1],myPCA300$x[,2])
PC1nul2 <- PC1

dev.off()

allGSDataSub <- allGSData[,which(names(allGSData) %in% psPC1means[,2])]
allGSDataSubCent <-  t(apply(allGSDataSub,1,function(x) x- psPC1means$means))
allGSDataSubCentProj <- as.matrix(allGSDataSubCent) %*% as.matrix(psPC1means$PC1)

plot(density(allGSDataSubCentProj[1:24]),xlim=c(-15,15),main="CEU Projection onto GS",xlab="PC1")
lines(density(allGSDataSubCentProj[25:156]),col=2)
lines(density(PC1nul1),col=3)
lines(density(PC1nul2),col=4)

myRows <- nrow(allGSDataSubCent)
addZeros <- allGSDataSub
for (j in 1:myRows){
	mySample <- sample(1:300,30)
	addZeros[j,mySample] <- 0
}
addZerosCent <-  t(apply(addZeros,1,function(x) x- psPC1means$means))
addZerosCentProj <- as.matrix(addZerosCent) %*% as.matrix(psPC1means$PC1)
plot(density(addZerosCentProj[1:24]),xlim=c(-15,15),main="CEU Projection onto GS",xlab="PC1")
lines(density(addZerosCentProj[25:156]),col=2)
lines(density(PC1nul1),col=3)
lines(density(PC1nul2),col=4)

myRows <- nrow(allGSDataSubCent)
addZeros <- allGSDataSub
for (j in 1:myRows){
	mySample <- rnorm(300)
	addZeros[j,] <- addZeros[j,] + mySample
}
addZerosCent <-  t(apply(addZeros,1,function(x) x- psPC1means$means))
addZerosCentProj <- as.matrix(addZerosCent) %*% as.matrix(psPC1means$PC1)
plot(density(addZerosCentProj[1:24]),xlim=c(-15,15),main="CEU Projection onto GS",xlab="PC1")
lines(density(addZerosCentProj[25:156]),col=2)
lines(density(PC1nul1),col=3)
lines(density(PC1nul2),col=4)

