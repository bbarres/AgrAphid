###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#loading and preparing the dataset
###############################################################################

#first, we load the genetic dataset
MyzAgra<-read.table("AgrAph.dat",header=T,sep="\t")
#here is the structure of the datafile, for explanation of each columns, see 
#ReadMe.txt file in DRYAD repository
head(MyzAgra)
#a summary of the different variables
summary(MyzAgra)
colnames(MyzAgra)
#total number of individuals
sum(table(MyzAgra$patch_ID)) #457 individuals

#let's remove the repeated MLGs in the dataset. We can easily do that by 
#using the 'dup' column of the dataset. To be conservative we remove every 
#repeated MLGs as well as non affected MLGs
MyzAgraccons<-MyzAgra[MyzAgra$dup=="o",]
JDD<-MyzAgraccons #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("orange","green","blue","yellow","hotpink")


###############################################################################
#DAPC on microsatellites only
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=6,ind.names=JDD$sample_ID, 
                    pop=JDD$year,ploidy=2)
#include the coordinates of the samples
JDDmicro@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<-find.clusters(JDDade,max.n.clust=30)
#with 50 PCs, we lost nearly no information and after K=5, the decrease of 
#the BIC value is smaller, so we chose the maximum number of clusters to be 5 
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#We try to optimize the number of principal component (PCs) to retain to 
#perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=15)
temp<-optim.a.score(dapcJDDade)
#the optimal number of PCs fell between 5 and 9 (depending on the run), so we 
#chose the smallest number of PCs (5) in order to avoid overfitting of the 
#model. Then we do the actual DAPC anlysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=3,n.pca=5)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
scatter(dapcJDDade,xax=1, yax=3,col=coloor)
scatter(dapcJDDade,xax=2, yax=3,col=coloor)
#a more beautifull scatter plot
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        main="Axis 1 & 2",solid=0.3,pch=19,cex=3,scree.da=FALSE)
scatter(dapcJDDade,xax=2,yax=3,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=FALSE)

#Run the 'find.clusters' and DAPC analysis for K=2 to 5
clustJDDade2<-find.clusters(JDDade,n.pca=50,n.clust=2)
dapcJDDade2<-dapc(JDDade,clustJDDade2$grp,n.da=1,n.pca=2)
compoplot(dapcJDDade2,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade3<-find.clusters(JDDade,n.pca=50,n.clust=3)
dapcJDDade3<-dapc(JDDade,clustJDDade3$grp,n.da=2,n.pca=3)
compoplot(dapcJDDade3,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade4<-find.clusters(JDDade,n.pca=50,n.clust=4)
dapcJDDade4<-dapc(JDDade,clustJDDade4$grp,n.da=3,n.pca=4)
compoplot(dapcJDDade4,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade5<-find.clusters(JDDade,n.pca=50,n.clust=5)
dapcJDDade5<-dapc(JDDade,clustJDDade5$grp,n.da=3,n.pca=5)
compoplot(dapcJDDade5,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)

write.table(dapcJDDade2$posterior,file="AgrAphDAPCK2.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade3$posterior,file="AgrAphDAPCK3.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade4$posterior,file="AgrAphDAPCK4.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade5$posterior,file="AgrAphDAPCK5.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)

barplot(t(dapcJDDade5$posterior),col=coloor,beside=FALSE,border=NA,
        space=0,ylim=c(-0.1,1.1),las=1,axisnames=FALSE,axes=FALSE)
#drawing an external rectangle
rect((0-dim(JDDade$tab)[1]/600),
     0-1/500,
     dim(JDDade$tab)[1]+dim(JDDade$tab)[1]/600,
     1+1/500,
     lwd=3)
#deliminated the different populations
rect(c(0,1,35),rep(0,3),c(1,35,46),rep(1,3),lwd=2)
#add some legend
mtext("K=5",side=2,las=1,cex=1.5,adj=0,line=2)

#export as a pdf file 12 X 3 inches


###############################################################################
#DAPC on microsatellites and resistance genotypes
###############################################################################

#converting data to a genind format including the resistance genotypes
JDDall<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                         "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                         "MP_4","MP_46","KDR","sKDR","MACE","R81T")],
                  ncode=6,ind.names=JDD$sample_ID, 
                  pop=JDD$year,ploidy=2)
#include the coordinates of the samples
JDDall@other$xy<-JDD[,c("longitude","latitude")]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDall
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=30,max.n.clust=35) #chose 5 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component 
#(PCs) to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=15)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=3,n.pca=10)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=pop(JDDade),legend=FALSE,
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
#Before importing the file, replace white space in the column header names 
#with underscore, replace "?1" by "alpha", and remove double white spaces or 
#it will provoc importation problem or failure

resstr_cccons<-read.table(file="AgrAphout.str", header=T,sep=" ",
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

deltastr_cccons<-chooseK(resstr_cccons,10,100)

#a function to plot variation of Delta K and Ln(P(X|K)) with K. 
plotdeltaK<-function(datadeltaK,nb_K,titre){
  #'datadeltak': the output file of 'chooseK' function
  #'nb_K': the number of different K considered
  #'titre': the title of the plot you want to be displayed
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

#the same function using log(deltaK), just in order to see smaller variation 
#of deltaK

#a function to plot variation of Delta K and Ln(P(X|K)) with K. 
plotlogdeltaK<-function(datadeltaK,nb_K,titre){
  #'datadeltak': the output file of 'chooseK' function
  #'nb_K': the number of different K considered
  #'titre': the title of the plot you want to be displayed
  op<-par(pty="s")
  plot(log(datadeltaK[1:(nb_K-2),8]+1),type="b",pch=24,cex=2.5,lwd=4,lty=1,
       col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(log(datadeltaK[1:(nb_K-2),8]+1),type="b",pch=24,bty="n",xaxt="n",yaxt="n",
       ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Ln(Delta K+1)",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}

op<-par(mfrow=c(1,2))
plotdeltaK(deltastr_cccons,10,
           "Conservative clone correction dataset (n=173)")
plotlogdeltaK(deltastr_cccons,10,
              "Conservative clone correction dataset (n=173)")
par(op)
#you can obtain the same figure as in the manuscript by exporting the plot to 
#png format, with a width of 2400 X 1100 pixels


###############################################################################
#trash
###############################################################################

#scatter plot with the different K groups and then plotting the population
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=FALSE)
#oilseed_rape
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,2],
       col="black",cex=2,bg="black",pch=21)
#peach
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,2],
       col="black",cex=2,bg="black",pch=21)
#tobacco
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,2],
       col="black",cex=2,bg="black",pch=21)
#other_crop
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,2],
       col="black",cex=2,bg="black",pch=21)


scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],
       col=coloor[dapcJDDade$assign],cex=2,
       pch=(as.numeric(as.factor(JDDade@other$host))+20))

scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],
       col=coloor[dapcJDDade$assign],pch=21,
       bg=coloor[(as.numeric(as.factor(JDDade@other$host)))])


plot(JDDade@other$xy,cex=3,col=dapcJDDade$assign,
     pch=as.numeric(as.factor(JDDade@other$host)))



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