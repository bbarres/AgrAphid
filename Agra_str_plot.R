###############################################################################
###############################################################################
#AgrAphid article's R code
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)

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
#we can also include the resistance genotypes as supplementary information
JDDmicro@other$KDR<-JDD[,"KDR"]
JDDmicro@other$sKDR<-JDD[,"sKDR"]
JDDmicro@other$MACE<-JDD[,"MACE"]
JDDmicro@other$R81T<-JDD[,"R81T"]

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

#alternatively, you can import a q-matrix output file from STRUCTURE software 
#and use the function in the same manner. Be careful howerer to respect the 
#order of the individuals and the order of their respective populations
strK2<-t(read.table("outK2.str",header=FALSE,sep="\t")[,c(-1)])
strK3<-t(read.table("outK3.str",header=FALSE,sep="\t")[,c(-1)])
strK4<-t(read.table("outK4.str",header=FALSE,sep="\t")[,c(-1)])
strK5<-t(read.table("outK5.str",header=FALSE,sep="\t")[,c(-1)])

#plot with different K values
op<-par(mfrow=c(4,1),mar=c(0,4,0,0),oma=c(3,0,0,0))
structplot(strK5,rainbow(5),effpop,poptiquet,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK4,rainbow(5),effpop,poptiquet,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK3,rainbow(5),effpop,poptiquet,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK2,rainbow(5),effpop,poptiquet,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord="grey70",
           distxax=0.08)
par(op)

#same plot with space between populations
coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2","darkorange")
op<-par(mfrow=c(4,1),mar=c(0,3,0,0),oma=c(8,0,0,0))
structplot(strK5,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK4,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK3,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(strK2,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord="grey70",
           distxax=0.15,angl=0,cexpop=1.5)
par(op)
#export to pdf 25 x 5 inches


###############################################################################
#plot a list of 100 STRUCTURE output files for each K
###############################################################################

#Usually, you run STRUCTURE several times for the same K values. After that, 
#you can reorganize the output file such as the labels of the different group 
#in the different run match (using CLUMPP for example). Here we import the 
#output file of CLUMPP and then we plot all the repetition in the same graph. 
#We first need to edit a little the output file in excel prior to the 
#importation: just keep the q matrix without any other information

#for K=2
K2_100runs<-read.table("AgrAccconsK2.perm_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")[,c(-1)]
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K2_100runs<-split(K2_100runs,rep(1:100,each=309))
K2_100runs[[1]]
coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  temp<-K2_100runs[[i]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=2",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=2",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=3
K3_100runs<-read.table("AgrAccconsK3.perm_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")[,c(-1)]
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K3_100runs<-split(K3_100runs,rep(1:100,each=309))
coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  temp<-K3_100runs[[i]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=3",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=3",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=4
K4_100runs<-read.table("AgrAccconsK4.perm_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")[,c(-1)]
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K4_100runs<-split(K4_100runs,rep(1:100,each=309))
coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  temp<-K4_100runs[[i]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=4",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=5
K5_100runs<-read.table("AgrAccconsK5.perm_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")[,c(-1)]
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K5_100runs<-split(K5_100runs,rep(1:100,each=309))
coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  temp<-K5_100runs[[i]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=5",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

###############################################################################
#plot of the clusterisation for different K value after CLUMPP averaging
###############################################################################

#the plot for the different K values
strK2<-t(read.table("AgrAccconsK2.outfile",header=FALSE,sep="\t")[,c(-1)])
strK3<-t(read.table("AgrAccconsK3.outfile",header=FALSE,sep="\t")[,c(-1)])
strK4<-t(read.table("AgrAccconsK4.outfile",header=FALSE,sep="\t")[,c(-1)])
strK5<-t(read.table("AgrAccconsK5.outfile",header=FALSE,sep="\t")[,c(-1)])

coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other\nCrops","Aerial Trap",
             "Multiple hosts")
op<-par(mfrow=c(4,1),mar=c(0,4,0,0),oma=c(5,0,0,0))
structplot(strK5[c(4,2,1,5,3),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK4[c(2,1,4,3),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK3[c(1,3,2),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK2,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.05,cexpop=1.5,angl=0)
par(op)
#export to pdf 21 X 7 inches


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
#or with the function we designed
poptiquet<-c("2000","2001","2002","2003","2004","2005","2006","2007")
structplot(t(dapcJDDade$posterior),coloor,summary(JDDade)$pop.eff,
           poptiquet,spacepop=2,leg_y="Assignement",cexy=1.2,mef=c(0,1,1,1,0),
           colbord="grey70",angl=0,distxax=0.001)
title(main="DAPC clusterisation",cex.main=1.5,outer=FALSE)
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
#burn-in period of 10.000 iterations followed by 100.000 simulations. One 
#hundred repetitions of each run were performed for K ranging from 1 to 10. 
#Before importing the file, replace white space in the column header names 
#with underscore, replace "?1" by "alpha", and remove double white spaces or 
#it will provoc importation problem or failure

#run the 'Agra_deltaKplot_fun.R' code before running this code

resstr_cccons<-read.table(file="AgrAphout.str", header=T,sep=" ",
                          blank.lines.skip=T)
deltastr_cccons<-chooseK(resstr_cccons,10,100)

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


###############################################################################
#END
###############################################################################