###############################################################################
###############################################################################
#Definition of fonctions for diversity indices computation
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)
library(vegan)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#Function for Allelic Richness computation
###############################################################################

#data: a dataset at the 'genind' format from the package 'adegenet'
AllRich<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean,
                         quiet = TRUE)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=row.names(datapop@tab), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele 
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Ar for each population/loci combination
  #(see El Mousadik and Petit 1996 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                        nCm(effTot[i,j],echMin)))
    }
  }
  
  #Allelic richness in each population, for each LOCUS
  ArLOC<-(tapply(CoMat, conf, sum))
  rez<-list("Minimum Sampling Size"=echMin,"Allelic Richness Matrix"=ArLOC)
  return(rez)
}

BRAt<-BRAcc #name of the input file
#converting data to a genind format
BRADE<-df2genind(BRAt[,14:27],ncode=3,ind.names=as.character(BRAt$sample_ID), 
                 pop=BRAt$pop_ID,NA.char=c("0"),ploidy=1)
BRADE@other$xy<-BRAt[,4:5]
AllRich(BRADE)[[2]]
Ar<-apply(AllRich(BRADE)[[2]],1,mean)


###############################################################################
#Function for Private Allelic Richness computation
###############################################################################

#data: a dataset at the 'genind' format from the package 'adegenet'
PrivAllRich<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE,other.action=mean,quiet=TRUE)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=row.names(datapop@tab), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  #   effN<- effN[order(as.numeric(rownames(effN))),]
  #   colnames(effN)<-locNames(datapop)
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele 
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Pijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                        nCm(effTot[i,j],echMin)))
    }
  }
  #fourth, compute the product of Qijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat2<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat2[i,j]<-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                      nCm(effTot[i,j],echMin))
    }
  }
  #fifth, compute the product of Kijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat3<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  temp<-c()
  for (i in 1:(dim(datapop@tab)[1])) {
    temp<-as.matrix(CoMat2[-i,])
    ifelse(dim(temp)[2]==1,CoMat3[i,]<-apply(temp,1,prod),
           CoMat3[i,]<-apply(temp,2,prod))
  }
  CoMat4<-CoMat*CoMat3
  #Private Allelic richness in each population, for each LOCUS
  PrivArLOC<-(tapply(CoMat4, conf, sum))
  #   PrivArLOC<-PrivArLOC[order(as.numeric(dimnames(PrivArLOC)[[1]])),]
  #   colnames(PrivArLOC)<-locNames(datapop)
  ##determining mean Allelic Richness across site and loci
  #determining mean Allelic Richness across loci
  #Ar<-(apply(ArLOC,1,mean))
  rez<-list("Minimum Sampling Size"=echMin,
            "Private Allelic Richness Matrix"=PrivArLOC)
  return(rez)
}

PrivAllRich(BRADE)
PrivAr<-apply(PrivAllRich(BRADE)[[2]],1,mean)


###############################################################################
#Function for Heterozygosity computation
###############################################################################

#data: a dataset at the 'genind' format
HeterNei<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean,
                         quiet = TRUE)
  #Heterozygosity (Nei 1987) in each population, for each LOCUS
  HsLOC<-matrix(nrow=(dim(datapop@tab)[1]),
                ncol=(length(levels(datapop@loc.fac))), byrow=TRUE)
  for (i in (1:(dim(datapop@tab)[1]))) {
    dataLOC<-genind2loci(data[data$pop==levels(data$pop)[i]])
    ss<-summary(dataLOC)
    HsLOC[i,]<-sapply(ss, function(x) H(x$allele))
  }
  #determining mean Heterozygosity across loci
  Hs<-(apply(HsLOC,1,mean))
  attr(Hs,"names")<-row.names(datapop@tab)
  return(Hs)
}

HetNei<-HeterNei(BRADE)


###############################################################################
#END
###############################################################################