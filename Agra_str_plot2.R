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
#plot a list of 100 STRUCTURE output files for each K
###############################################################################

#Usually, you run STRUCTURE several times for the same K values. After that, 
#you can reorganize the output file such as the labels of the different group 
#in the different run match (using CLUMPP for example). Here we import the 
#output file of CLUMPAK (http://clumpak.tau.ac.il/index.html) and then we plot 
#all the repetitions in the same graph  
#We first need to edit a little the output file in excel prior to the 
#importation: just keep the q matrix without any other information

#for K=2, import the 100 run datafile (100 q matrix)
K2_100runs<-read.table("AgrAccconsK2.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K2_colorder<-read.table("AgrAccconsK2.ind_colorder",header=FALSE,
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
K2_reporder<-read.table("AgrAccconsK2.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  j<-as.numeric(K2_reporder[i])
  temp<-K2_100runs[[j]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=2",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=2",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=3, import the 100 run datafile (100 q matrix)
K3_100runs<-read.table("AgrAccconsK3.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K3_colorder<-read.table("AgrAccconsK3.ind_colorder",header=FALSE,
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
K3_reporder<-read.table("AgrAccconsK3.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  j<-as.numeric(K3_reporder[i])
  temp<-K3_100runs[[j]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=3",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=3",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=4, import the 100 run datafile (100 q matrix)
K4_100runs<-read.table("AgrAccconsK4.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K4_colorder<-read.table("AgrAccconsK4.ind_colorder",header=FALSE,
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
K4_reporder<-read.table("AgrAccconsK4.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  j<-as.numeric(K4_reporder[i])
  temp<-K4_100runs[[j]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=4",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=4",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=5, import the 100 run datafile (100 q matrix)
K5_100runs<-read.table("AgrAccconsK5.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K5_colorder<-read.table("AgrAccconsK5.ind_colorder",header=FALSE,
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
K5_reporder<-read.table("AgrAccconsK5.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  j<-as.numeric(K5_reporder[i])
  temp<-K5_100runs[[j]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=5",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=5",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12

#for K=6, import the 100 run datafile (100 q matrix)
K6_100runs<-read.table("AgrAccconsK6.ind_datafile",header=FALSE,
                       blank.lines.skip=TRUE,sep="\t")
#import the column order for the best CLUMPP "permutation"
K6_colorder<-read.table("AgrAccconsK6.ind_colorder",header=FALSE,
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
K6_reporder<-read.table("AgrAccconsK6.ind_linorder",header=FALSE,
                        blank.lines.skip=TRUE,sep="\t")+1

coloor <- c("chartreuse4","firebrick","khaki2","darkorange","royalblue4",
            "grey30")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
#now we can plot the 100 runs on the same figure
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  j<-as.numeric(K6_reporder[i])
  temp<-K6_100runs[[j]]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=6",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=6",cex.main=2.5,outer=TRUE)
par(op)
#export pdf 25 x 12


###############################################################################
#plot the different repetitions of structure runs grouped by solutions
###############################################################################

#for K=6, major clustering solutions
K6_maj<-read.table("ClumppIndFileK6maj",header=FALSE,
                   blank.lines.skip=TRUE,sep=" ")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K6_maj<-split(K6_maj[,6:11],rep(1:55,each=309))
coloor <- c("chartreuse4","royalblue4","firebrick","grey20","khaki2",
            "darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
K6majord<-read.table("K6majclustord.txt",header=FALSE,sep=" ")
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:55){
  temp<-K6_maj[[i]]
  temp<-temp[,as.numeric(K6majord[i,])]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=6",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=6",cex.main=2.5,outer=TRUE)
par(op)

#for K=6, all runs
K6_<-read.table("ClumppIndFileK6",header=FALSE,
                blank.lines.skip=TRUE,sep=" ")
#then we split the dataframe in as many repetition that has been made
#by the number of individuals (here 309)
K6_<-split(K6_[,6:11],rep(1:100,each=309))
coloor <- c("chartreuse4","royalblue4","firebrick","grey20","khaki2",
            "darkorange")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other Crops","Aerial Trap",
             "Multiple hosts")
K6ord<-read.table("K6clustord.txt",header=FALSE,sep=" ")
op<-par(mfrow=c(100,1),mar=c(0,0,0,0),oma=c(1,0,3,0))
for (i in 1:100){
  temp<-K6_[[i]]
  temp<-temp[,as.numeric(K6ord[i,])]
  structplot(t(temp),coloor,effpop,poptiquet,spacepop=2,
             leg_y="K=6",cexy=1.2,mef=c(0,0,0,0,0),colbord=NA,
             distxax=0.15,angl=0,cexpop=1.5)
}
title(main="K=6",cex.main=2.5,outer=TRUE)
par(op)



###############################################################################
#Identifying the best K for STRUCTURE run
###############################################################################

#Analyzes were performed using STRUCTURE2.3.4 software, with a model allowing 
#admixture and correlation of allele frequencies. Each run consisted of a 
#burn-in period of 100.000 iterations followed by 1.000.000 simulations. One 
#hundred repetitions of each run were performed for K ranging from 1 to 15. 
#Before importing the file, replace white space in the column header names 
#with underscore, replace "?1" by "alpha", and remove double white spaces or 
#it will provoc importation problem or failure

#run the 'Agra_deltaKplot_fun.R' code before running this code

resstr_cccons<-read.table(file="AgrAphout2.str", header=T,sep=" ",
                          blank.lines.skip=T)
deltastr_cccons<-chooseK(resstr_cccons,15,100)

plotdeltaK(deltastr_cccons,15,
           "Conservative clone correction dataset (n=309)")

#you can obtain the same figure as in the manuscript by exporting the plot to 
#png format, with a width of 2400 X 1100 pixels


###############################################################################
#Plotting the best K values for the complete dataset
###############################################################################

#for this study, K=3 and K=5 are the most interesting values
coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
poptiquet<-c("Peach","Oilseed\nrape","Tobacco","Other\nCrops","Aerial Trap",
             "Multiple hosts")
#let's order the dataset as we want
ordatAgra<-datAgra
levels(ordatAgra$host_corrected)<-c(5,2,4,1,6,3)
levels(ordatAgra$host)<-c(5,2,4,1,3)
ordatAgra<-ordatAgra[order(as.numeric(as.character(ordatAgra$host)),
                           ordatAgra$year,ordatAgra$sampling_date),]
ordatAgracc<-ordatAgra[ordatAgra$one_MLG==1,]
ordatAgracc<-ordatAgracc[order(as.numeric(as.character
                                          (ordatAgracc$host_corrected)),
                               ordatAgracc$year,ordatAgracc$sampling_date),]

#plot for the complete clone-corrected dataset
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
effpop<-table(ordatAgracc$host_corrected)[c(4,2,6,3,1,5)]
#K=3
structplot(t(ordatAgracc[,39:41]),coloor[c(1,2,3,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(ordatAgracc[,43:47]),coloor[c(1,3,2,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
par(op)

#export to pdf 20 x 5 inches

#plot for the complete dataset
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
#plot for K=5 with the complete dataset
effpop<-table(ordatAgra$host)[c(4,2,5,3,1)]
#K=3
structplot(t(ordatAgra[,39:41]),coloor[c(1,2,3,4,5)],effpop,poptiquet,
           spacepop=8,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(ordatAgra[,43:47]),coloor[c(1,3,2,4,5)],effpop,poptiquet[-6],
           spacepop=8,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
par(op)

#export to pdf 25 x 5 inches


###############################################################################
#Plotting the best K values for the temporal dataset
###############################################################################

tempordatAgra<-ordatAgra[ordatAgra$data_batch=="AgrAphid"
                         & !is.na(ordatAgra$data_batch),]
tempordatAgracc<-tempordatAgra[tempordatAgra$one_MLG_year==1,]

#plot the temporal dataset by year
poptiquet<-names(table(tempordatAgra$year))
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
effpop<-table(tempordatAgra$year)
#K=3
structplot(t(tempordatAgra[,39:41]),coloor[c(1,2,3,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(tempordatAgra[,43:47]),coloor[c(1,3,2,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.15,angl=45,cexpop=1.2)
par(op)

#export to pdf 20 x 5 inches

#plot the temporal dataset by semester
poptiquet<-names(table(tempordatAgra$semester))
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
effpop<-table(tempordatAgra$semester)
#K=3
structplot(t(tempordatAgra[-c(1:2),39:41]),coloor[c(1,2,3,4,5)],effpop,
           poptiquet,spacepop=2,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),
           colbord=NA,distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(tempordatAgra[-c(1:2),43:47]),coloor[c(1,3,2,4,5)],effpop,
           poptiquet,spacepop=2,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),
           colbord=NA,distxax=0.15,angl=45,cexpop=1.2)
par(op)

#export to pdf 20 x 5 inches

#plot the clone-corrected temporal dataset by year
poptiquet<-names(table(tempordatAgracc$year))
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
effpop<-table(tempordatAgracc$year)
#K=3
structplot(t(tempordatAgracc[,39:41]),coloor[c(1,2,3,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),colbord=NA,
           distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(tempordatAgracc[,43:47]),coloor[c(1,3,2,4,5)],effpop,poptiquet,
           spacepop=2,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.15,angl=45,cexpop=1.2)
par(op)

#export to pdf 20 x 5 inches

#plot the clone-corrected temporal dataset by semester
poptiquet<-names(table(tempordatAgracc$semester))
op<-par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(4,0,0,0))
effpop<-table(tempordatAgracc$semester)
#K=3
structplot(t(tempordatAgracc[-c(1),39:41]),coloor[c(1,2,3,4,5)],effpop,
           poptiquet,spacepop=2,leg_y="K=3",cexy=1.2,mef=c(0,1,1,0,1),
           colbord=NA,distxax=0.15,angl=0,cexpop=1.2)
#K=5
structplot(t(tempordatAgracc[-c(1),43:47]),coloor[c(1,3,2,4,5)],effpop,
           poptiquet,spacepop=2,leg_y="K=5",cexy=1.2,mef=c(0,1,1,1,1),
           colbord=NA,distxax=0.15,angl=45,cexpop=1.2)
par(op)

#export to pdf 20 x 5 inches


###############################################################################
#Structure-like plot
###############################################################################

#some examples of the use of the function you can load from 
#'Agra_strplot_fun.R'

#you can import a q-matrix output file from STRUCTURE software 
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
#plot of the clusterisation for different K value after CLUMPP averaging
###############################################################################

#the plot for the different K values
strK2<-t(read.table("AgrAccconsK2.outfile",header=FALSE,sep="\t")[,c(-1)])
strK3<-t(read.table("AgrAccconsK3.outfile",header=FALSE,sep="\t")[,c(-1)])
strK4<-t(read.table("AgrAccconsK4.outfile",header=FALSE,sep="\t")[,c(-1)])
strK5<-t(read.table("AgrAccconsK5.outfile",header=FALSE,sep="\t")[,c(-1)])
strK6<-t(read.table("AgrAccconsK6.outfile",header=FALSE,sep="\t")[,c(-1)])

coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange",
            "grey20")
effpop<-c(69,29,11,16,168,16)
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other\nCrops","Aerial Trap",
             "Multiple hosts")
op<-par(mfrow=c(5,1),mar=c(0,4,0,0),oma=c(5,0,0,0))
structplot(strK6[c(2,3,1,6,5,4),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=6",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK5[c(5,2,1,4,3),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK4[c(2,4,3,1),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK3[c(1,3,2),],coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK2,coloor,effpop,poptiquet,spacepop=2,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.05,cexpop=1.5,angl=0)
par(op)
#export to pdf 21 X 7 inches


###############################################################################
#END
###############################################################################