###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/PeachMyz_data")


###############################################################################
#loading and preparing the dataset
###############################################################################

#first, we load the genetic dataset
MyzPeach<-read.table("MyzPeach.dat",header=T,sep="\t")
#here is the structure of the datafile, for explanation of each columns, see 
#ReadMe.txt file in DRYAD repository
head(MyzPeach)
#a summary of the different variables
summary(MyzPeach)
colnames(MyzPeach)
#number of individuals in each sampled populations
table(MyzPeach$patch_ID)
#total number of individuals
sum(table(MyzPeach$patch_ID)) #312 individuals


#here you can select the sub dataset you want, for a first analysis, we take
#the entire dataset
JDD<-MyzPeach #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("red","green","blue","yellow","orchid")


###############################################################################
#DAPC on microsatellites only
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=6,ind.names=JDD$sample_ID, 
                    pop=JDD$patch_ID,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDmicro@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component 
#(PCs) to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=10)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=truenames(JDDade)$pop,legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)


###############################################################################
#DAPC on microsatellites and resistance genotypes
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDall<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                         "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                         "MP_4","MP_46","kdr","skdr","R81T","MACE")],
                  ncode=6,ind.names=JDD$sample_ID, 
                  pop=JDD$patch_ID,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDall@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDall
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component 
#(PCs) to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=7)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=truenames(JDDade)$pop,legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)


###############################################################################
#Identifying the best K for STRUCTURE run
###############################################################################

#Analyzes were performed using STRUCTURE2.3.4 software, with a model allowing 
#admixture and correlation of allele frequencies. Each run consisted of a 
#burn-in period of 10.000 iterations followed by 50.000 simulations. One 
#hundred repetitions of each run were performed for K ranging from 1 to 10. 
#Before importing the file, replace white space in the column header names with 
#underscore, replace "?1" by "alpha", and remove double white spaces or it
#will provoc importation problem or failure

resstr<-read.table(file="BRAoutput.out", header=T,sep=" ",blank.lines.skip=T)
resccstr<-read.table(file="BRAccoutput.out", header=T,sep=" ",
                     blank.lines.skip=T)
resccconsstr<-read.table(file="BRAccconsoutput.out", header=T,sep=" ",
                         blank.lines.skip=T)


#a function which compute delta K values, nb_K is the number of different K 
#considered, and nb_rep is the number of repetition of each K
chooseK<-function(str_out,nb_K,nb_rep) {
  datatable<-data.frame("K"=c(rep(1:nb_K,each=nb_rep)),"Ln(Pd)"=str_out[,4])
  Lprim<-c(rep("NA",nb_rep))
  for (i in ((nb_rep+1):(nb_K*nb_rep))) {
    Lprim<-c(Lprim,str_out[i,4]-str_out[i-nb_rep,4])
  }
  datatable<-data.frame(datatable,as.numeric(Lprim))
  Lsecond<-c(rep("NA",nb_rep))
  for (i in (((2*nb_rep)+1):(nb_K*nb_rep))) {
    Lsecond<-c(Lsecond,abs(datatable[i,3]-datatable[i-nb_rep,3]))
  }
  Lsecond<-c(Lsecond,rep("NA",nb_rep))
  datatable<-data.frame(datatable,as.numeric(Lsecond))
  reztable<-data.frame("K"=c(1:nb_K))
  meanL<-c()
  sdL<-c()
  for (i in (1:nb_K)) {
    meanL<-c(meanL,mean(datatable[datatable$K==i,2]))
    sdL<-c(sdL,sd(datatable[datatable$K==i,2]))
  }
  reztable<-data.frame(reztable,meanL,sdL)
  meanLprime<-c()
  sdLprime<-c()
  for (i in (1:nb_K)) {
    meanLprime<-c(meanLprime,mean(as.numeric(datatable[datatable$K==i,3])))
    sdLprime<-c(sdLprime,sd(datatable[datatable$K==i,3]))
  }
  reztable<-data.frame(reztable,meanLprime,sdLprime)
  meanLsecond<-c()
  sdLsecond<-c()
  for (i in (1:nb_K)) {
    meanLsecond<-c(meanLsecond,mean(as.numeric(datatable[datatable$K==i,4])))
    sdLsecond<-c(sdLsecond,sd(datatable[datatable$K==i,4]))
  }
  reztable<-data.frame(reztable,meanLsecond,sdLsecond)
  deltaK<-c()
  for (i in (1:nb_K)) {
    deltaK<-c(deltaK,reztable[reztable$K==i,6]/reztable[reztable$K==i,3])
  }
  reztable<-data.frame(reztable,deltaK)
  return(reztable)
}

deltastr<-chooseK(resstr,15,15)
deltaccstr<-chooseK(resccstr,15,15)
deltaccconsstr<-chooseK(resccconsstr,15,15)

#a function to plot variation of Delta K and Ln(P(X|K)) with K. 'datadeltak' 
#is the output file of 'chooseK' function, and nb_K is the number of different
#K considered
plotdeltaK<-function(datadeltaK,nb_K,titre){
  op<-par(pty="s")
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,cex=2.5,lwd=4,lty=1,
       col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,bty="n",xaxt="n",yaxt="n",
       ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Delta K",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}

plotdeltaK(deltastr,15,"Complete dataset (n=474)")
plotdeltaK(deltaccstr,15,"Clone Corrected dataset (n=338)")
plotdeltaK(deltaccconsstr,15,"Conservative Clone Corrected dataset (n=259)")

#you can obtain the same figure as in the manuscript by exporting the plot to 
#pdf format, with a width of 12 inches and an height of 11 inches

#You can also obtain a combined plot for the three different dataset
op<-par(mfrow=c(1,3))
plotdeltaK(deltastr,15,"Complete dataset (n=474)")
plotdeltaK(deltaccstr,15,"Clone Corrected dataset (n=338)")
plotdeltaK(deltaccconsstr,15,"Conservative Clone Corrected dataset (n=259)")
par(op) 
#then export with a width of 32 inches and an height of 10 inches




###############################################################################
#trash
###############################################################################

#scatter plot with the different K groups and then plotting the population
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=FALSE)
#oilseed_rape
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,2],col="black",
       ,cex=2,bg="black",pch=21)
#peach
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,2],col="black",
       ,cex=2,bg="black",pch=21)
#tobacco
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,2],col="black",
       ,cex=2,bg="black",pch=21)
#other_crop
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,2],col="black",
       ,cex=2,bg="black",pch=21)


scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=(as.numeric(as.factor(JDDade@other$host))+20),cex=2)

scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=21,bg=coloor[(as.numeric(as.factor(JDDade@other$host)))])


plot(JDDade@other$xy,cex=3,col=dapcJDDade$assign,pch=as.numeric(as.factor(JDDade@other$host)))



# BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")
# 
# image(alt,col=brewer.pal(9,"Greys"))
# stars(table(pop(JDDade),dapcJDDade$assign),draw.segment=TRUE,
#       locations=JDDade@other$xy,
#       #locations=cbind(jitter(BRADEpop@other$xy$longitude,200),
#       #                jitter(BRADEpop@other$xy$latitude,200)),
#       add=T,len=0.5)



###############################################################################
#END
###############################################################################