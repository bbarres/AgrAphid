##############################################################################/
##############################################################################/
#Plotting the results of STRUCTURE 100 runs
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#plot a list of 100 STRUCTURE output files for each K####
##############################################################################/

#Usually, you run STRUCTURE several times for the same K values. After that, 
#you can reorganize the output file such as the labels of the different group 
#in the different run match (using CLUMPP for example). Here we import the 
#output file of CLUMPAK (http://clumpak.tau.ac.il/index.html) and then we plot 
#all the repetitions in the same graph  
#We first need to edit a little the output file in excel prior to the 
#importation: just keep the q matrix without any other information

#for K=2, import the 100 run datafile (100 q matrix)
K2_100runs<-read.table("data/AgrAccconsK2.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K2_colorder<-read.table("data/AgrAccconsK2.ind_colorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K2_100runs<-split(K2_100runs,rep(1:100,each=309))
#reordering the columns so that the different repetition colorization fit
for (i in 1:100){
  K2_100runs[[i]]<-K2_100runs[[i]][as.numeric(K2_colorder[i,])]
}
#importing the order of the run so that the different repetition corresponding
#to the same clustering solution followed each other
K2_reporder<-read.table("data/AgrAccconsK2.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor<-c("chartreuse4","firebrick","khaki2","darkorange","royalblue4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0.2,0,0,0),oma=c(7,0,3,0))
for (i in 1:99){
  j<-as.numeric(K2_reporder[i])
  temp<-K2_100runs[[j]]
  structplot(t(temp),coloor[c(2,1)],effpop,poptiquet,spacepop=2,
             leg_y="K=2",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
j<-as.numeric(K2_reporder[100])
temp<-K2_100runs[[j]]
structplot(t(temp),coloor[c(2,1)],effpop,poptiquet,spacepop=2,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,
           distxax=1,angl=20,cexpop=2.5)
title(main="K=2",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 15 x 30

#for K=3, import the 100 run datafile (100 q matrix)
K3_100runs<-read.table("data/AgrAccconsK3.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K3_colorder<-read.table("data/AgrAccconsK3.ind_colorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K3_100runs<-split(K3_100runs,rep(1:100,each=309))
#reordering the columns so that the different repetition colorization fit
for (i in 1:100){
  K3_100runs[[i]]<-K3_100runs[[i]][as.numeric(K3_colorder[i,])]
}
#importing the order of the run so that the different repetition corresponding
#to the same clustering solution followed each other
K3_reporder<-read.table("data/AgrAccconsK3.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","royalblue4","khaki2","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0.2,0,0,0),oma=c(7,0,3,0))
for (i in 1:99){
  j<-as.numeric(K3_reporder[i])
  temp<-K3_100runs[[j]]
  structplot(t(temp),coloor[c(2,1,3)],effpop,poptiquet,spacepop=2,
             leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
j<-as.numeric(K3_reporder[100])
temp<-K3_100runs[[j]]
structplot(t(temp),coloor[c(2,1,3)],effpop,poptiquet,spacepop=2,
           leg_y="K=3",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,
           distxax=1,angl=20,cexpop=2.5)
title(main="K=3",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 15 x 30

#for K=4, import the 100 run datafile (100 q matrix)
K4_100runs<-read.table("data/AgrAccconsK4.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K4_colorder<-read.table("data/AgrAccconsK4.ind_colorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K4_100runs<-split(K4_100runs,rep(1:100,each=309))
#reordering the columns so that the different repetition colorization fit
for (i in 1:100){
  K4_100runs[[i]]<-K4_100runs[[i]][as.numeric(K4_colorder[i,])]
}
#importing the order of the run so that the different repetition corresponding
#to the same clustering solution followed each other
K4_reporder<-read.table("data/AgrAccconsK4.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("royalblue4","khaki2","firebrick","chartreuse4","darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0.2,0,0,0),oma=c(7,0,3,0))
for (i in 1:99){
  j<-as.numeric(K4_reporder[i])
  temp<-K4_100runs[[j]]
  structplot(t(temp),coloor[c(1,2,4,3)],effpop,poptiquet,spacepop=2,
             leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
j<-as.numeric(K4_reporder[100])
temp<-K4_100runs[[j]]
structplot(t(temp),coloor[c(1,2,4,3)],effpop,poptiquet,spacepop=2,
           leg_y="K=4",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,
           distxax=1,angl=20,cexpop=2.5)
title(main="K=4",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 15 x 30

#for K=5, import the 100 run datafile (100 q matrix)
K5_100runs<-read.table("data/AgrAccconsK5.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K5_colorder<-read.table("data/AgrAccconsK5.ind_colorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K5_100runs<-split(K5_100runs,rep(1:100,each=309))
#reordering the columns so that the different repetition colorization fit
for (i in 1:100){
  K5_100runs[[i]]<-K5_100runs[[i]][as.numeric(K5_colorder[i,])]
}
#importing the order of the run so that the different repetition corresponding
#to the same clustering solution followed each other
K5_reporder<-read.table("data/AgrAccconsK5.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor<-c("darkorange","khaki2","royalblue4","firebrick","chartreuse4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0.2,0,0,0),oma=c(7,0,3,0))
for (i in 1:99){
  j<-as.numeric(K5_reporder[i])
  temp<-K5_100runs[[j]]
  structplot(t(temp),coloor[c(1,2,3,5,4)],effpop,poptiquet,spacepop=2,
             leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
j<-as.numeric(K5_reporder[100])
temp<-K5_100runs[[j]]
structplot(t(temp),coloor[c(1,2,3,5,4)],effpop,poptiquet,spacepop=2,
           leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,
           distxax=1,angl=20,cexpop=2.5)
title(main="K=5",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 15 x 30

#for K=6, import the 100 run datafile (100 q matrix)
K6_100runs<-read.table("data/AgrAccconsK6.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K6_colorder<-read.table("data/AgrAccconsK6.ind_colorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K6_100runs<-split(K6_100runs,rep(1:100,each=309))
#reordering the columns so that the different repetition colorization fit
for (i in 1:100){
  K6_100runs[[i]]<-K6_100runs[[i]][as.numeric(K6_colorder[i,])]
}
#importing the order of the run so that the different repetition corresponding
#to the same clustering solution followed each other
K6_reporder<-read.table("data/AgrAccconsK6.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4",
            "grey50")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0.2,0,0,0),oma=c(7,0,3,0))
for (i in 1:99){
  j<-as.numeric(K6_reporder[i])
  temp<-K6_100runs[[j]]
  structplot(t(temp),coloor[c(5,2,3,1,4,6)],effpop,poptiquet,spacepop=2,
             leg_y="K=6",cexy=1.2,mef=c(0,1,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
j<-as.numeric(K6_reporder[100])
temp<-K6_100runs[[j]]
structplot(t(temp),coloor[c(5,2,3,1,4,6)],effpop,poptiquet,spacepop=2,
           leg_y="K=6",cexy=1.2,mef=c(0,1,1,1,0),colbord=NA,
           distxax=1,angl=20,cexpop=2.5)
title(main="K=6",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 15 x 30


##############################################################################/
#END
##############################################################################/