##############################################################################/
##############################################################################/
#Plotting the results of STRUCTURE post CLUMPP
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Identifying the best K for STRUCTURE run####
##############################################################################/

#Analyzes were performed using STRUCTURE2.3.4 software, with a model allowing 
#admixture and correlation of allele frequencies. Each run consisted of a 
#burn-in period of 100.000 iterations followed by 1.000.000 simulations. One 
#hundred repetitions of each run were performed for K ranging from 1 to 15. 
#Before importing the file, replace white space in the column header names 
#with underscore, replace "?1" by "alpha", and remove double white spaces or 
#it will provoke importation problem or failure

#run the 'Agra_deltaKplot_fun.R' code before running this code

resstr_cccons<-read.table(file="data/AgrAphout2.str",header=T,sep=" ",
                          blank.lines.skip=T)
deltastr_cccons<-chooseK(resstr_cccons,15,100)

op<-par(mar=c(5.1,5.1,4.1,6.1))
plotdeltaK(deltastr_cccons,15,
           "Clone corrected dataset (n=309)")
par(op)

#export to .pdf 8 x 7 inches


##############################################################################/
#plot of the clusterisation for different K value after CLUMPP averaging####
##############################################################################/

#the plot for the different K values
strK2<-t(read.table("data/AgrAccconsK2.outfile",header=FALSE,sep="\t")[,c(-1)])
strK3<-t(read.table("data/AgrAccconsK3.outfile",header=FALSE,sep="\t")[,c(-1)])
strK4<-t(read.table("data/AgrAccconsK4.outfile",header=FALSE,sep="\t")[,c(-1)])
strK5<-t(read.table("data/AgrAccconsK5.outfile",header=FALSE,sep="\t")[,c(-1)])
strK6<-t(read.table("data/AgrAccconsK6.outfile",header=FALSE,sep="\t")[,c(-1)])

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


##############################################################################/
#END
##############################################################################/

