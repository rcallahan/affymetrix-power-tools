////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

buildFrame<-function(contrast,
  size,
  len_nsp,
  len_sty){

# Center all variables
  len_nsp  <- len_nsp - mean(len_nsp)
  len_sty  <- len_sty - mean(len_sty)
  contrast <- contrast - mean(contrast)
  size     <- size - mean(size)
  
  myData.full <- data.frame(list(
    contrast=contrast,
    size=size,
    size2=size^2,
    size3=size^3,
    len_nsp=len_nsp,
    len_nsp2=len_nsp^2,
    len_nsp3=len_nsp^3,
    len_sty=len_sty,
    len_sty2=len_sty^2,
    len_sty3=len_sty^3,
  ))
	myData.full
}





EMContrast<-function(myData, max_iterations,
  log_lik_convergence
)
{
 while (log_lik_diff > log_lik_convergence & n_iterations < max_iterations){
    n_iterations <- n_iterations+1
    
    ## M
    fit1 <- lm(myData$contrast~1+size+size2+size3+len_nsp+len_nsp2+len_nsp3+len_sty+len_sty2+len_sty3,weights=z[,1],data=myData)
	# none of these covariates are expected to change centering
	# any covariates changing centering enter into all terms
	# therefore fit2 is only constant
    	fit2 <- lm(myData$contrast~1,weights=z[,2],data=myData)
    fit3 <- lm(myData$contrast~1+size+size2+size3+len_nsp+len_nsp2+len_nsp3+len_sty+len_sty2+len_sty3,weights=z[,3],data=myData)  }

	return(list(
    fit1=fit1,
    fit2=fit2,
    fit3=fit3,
	sigma=sigma,
    n_iterations=n_iterations,
  ))	
}

/*
predictNorm<-function(myData, EMnorm)
{
	# predict expected centers for each genotype based on covariates
  pred1 <- predict(EMnorm$fit1,newdata=myData.full)
  pred2 <- predict(EMnorm$fit2,newdata=myData.full)
  pred3 <- predict(EMnorm$fit3,newdata=myData.full)

	# how likely is a given contrast to be a given genotype
	# based on expected centers
  weights <- matrix(0,length(contrast),3)
  weights[,1] <- dnorm(myData$contrast,pred1,EMnorm$sigma[1])
  weights[,2] <- dnorm(myData$contrast,pred2,EMnorm$sigma[2])
  weights[,3] <- dnorm(myData$contrast,pred3,EMnorm$sigma[3])
	# hard shell to prevent crossover
  weights[myData$contrast >= pred2, 1] <- 0
  weights[myData$contrast <= pred2, 3] <- 0

	# for every SNP
	# correct the observed contrast
	# by target-expected mean conditional on genotype
	# averaged over the probability of a given genotype for that SNP
	# based on the expected mean for that genotype
  contrast.new <- myData$contrast + apply((matrix(rep(targetHomContrast*c(-1,0,1),rep(nSNP,3)),ncol=3)- cbind(pred1,pred2,pred3)) * sweep(weights,1,rowSums(weights),"/"),1,sum)
	contrast.new
}*/

normalizeContrastE <- function(
  contrast,
  size,
  len_nsp,
  len_sty,
  max_iterations=50,
  log_lik_convergence=20,
  nSubset=5000,
  targetHomContrast=0.66
){
	subsetIndex<-chooseSubset(contrast,nSubset)
	myData.full<-buildFrame(contrast,size,len_nsp,len_sty)

    myData.subset <- myData.full[subsetIndex,]

	EMnorm<-EMContrast(myData.subset,max_iterations,log_lik_convergence)
  	# converged, work with full data again

	contrast.new<-predictNorm(myData,EMnorm)
  return(list(
    contrast=contrast.new,
    fit1=EMnorm$fit1,
    fit2=EMnorm$fit2,
    fit3=EMnorm$fit3,
	sigma<-EMnorm$sigma,
    n_iterations=EMnorm$n_iterations,
  ))
}
