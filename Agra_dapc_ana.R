###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#Only Temporal data: loading and preparing the dataset
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
#we can also include the resistance genotypes as supplementary information
JDDmicro@other$KDR<-JDD[,"KDR"]
JDDmicro@other$sKDR<-JDD[,"sKDR"]
JDDmicro@other$MACE<-JDD[,"MACE"]
JDDmicro@other$R81T<-JDD[,"R81T"]

#now we format the dataset to analyse it with dapc from the adegenet package
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

#a more beautifull scatter plot, the colors are matching colors used in 
#the structure-like plot
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
scatter(dapcJDDade5,xax=2,yax=3,cstar=1,cell=0,clab=0,main="Axis 2 & 3",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)

#the same plot with resistance genotypes added
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
#adding KDR resistotype
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$KDR)))
legend("topright",levels(as.factor(JDDmicro@other$KDR)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$KDR)))),
       xpd=NA)
title("KDR resistotypes")

#adding sKDR resistotype
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$sKDR)))
legend("topright",levels(as.factor(JDDmicro@other$sKDR)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$sKDR)))),
       xpd=NA)
title("sKDR resistotypes")

#adding MACE resistotype
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$MACE)))
legend("topright",levels(as.factor(JDDmicro@other$MACE)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$MACE)))),
       xpd=NA)
title("MACE resistotypes")

#in case we need the q-matrix of the individuals for other purposes
write.table(dapcJDDade2$posterior,file="AgrAphDAPCK2.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade3$posterior,file="AgrAphDAPCK3.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade4$posterior,file="AgrAphDAPCK4.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade5$posterior,file="AgrAphDAPCK5.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)


###############################################################################
#Structure-like plot
###############################################################################

#some examples of the use of the function you can load from 
#'Agra_strplot_fun.R'

#first you need to gather the number of individuals in each populations
effpop<-as.numeric(table(JDDade$pop))
#the names of the different populations might be useful too
poptiquet<-levels(JDDade$pop)
#be careful to use the same dataset that has been used for the DAPC 
#computation
structplot(t(dapcJDDade5$posterior),rainbow(5),effpop,poptiquet)
structplot(t(dapcJDDade5$posterior),rainbow(5),effpop,poptiquet,
           colbord="grey70",leg_y="K=5",angl=30,distxax=0.01)
structplot(t(dapcJDDade5$posterior),coloor,effpop,poptiquet,
           leg_y="K=5",mef=c(1,0,1,1,1),cexpop=0.5,cexy=5)
structplot(t(dapcJDDade2$posterior),coloor,effpop,poptiquet,
           colbord=0,leg_y="K=2",mef=c(0,1,0,1,1))

#Now, we can easily plot several structure-like plot in the same figure
op<-par(mfrow=c(4,1),mar=c(0,4,0,0),oma=c(3,0,0,0))
structplot(t(dapcJDDade5$posterior)[c(4,1,2,5,3),],rainbow(5),effpop,poptiquet,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade4$posterior)[c(1,4,3,2),],rainbow(5),effpop,poptiquet,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade3$posterior)[c(1,3,2),],rainbow(5),effpop,poptiquet,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade2$posterior),rainbow(5),effpop,poptiquet,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord="grey70",
           distxax=0.08)
par(op)
#export to pdf 15 X 4 inches