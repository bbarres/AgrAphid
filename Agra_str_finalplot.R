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

#preparing the dataset
temp<-as.data.table(datAgraccHost)
#reformatting KDR genotypes
temp$KDRg<-temp$KDR
levels(temp$KDRg)<-c("K-RR","K-RS","K-SS")
temp$KDRg<-as.character(temp$KDRg)
temp$KDRg[is.na(temp$KDRg)]<-"K-miss"
temp<-spread(temp,KDRg,KDRg)
#reformatting sKDR genotypes
temp$sKDRg<-temp$sKDR
levels(temp$sKDRg)<-c("sK-RR","sK-RR","sK-RR","sK-RS","sK-RS",
                      "sK-SS","sK-RS","sK-RS","sK-RS")
temp$sKDRg<-as.character(temp$sKDRg)
temp$sKDRg[is.na(temp$sKDRg)]<-"sK-miss"
temp<-spread(temp,sKDRg,sKDRg)
#reformatting MACE genotypes
temp$MACEg<-temp$MACE
levels(temp$MACEg)<-c("M-SS","M-RS")
temp$MACEg<-as.character(temp$MACEg)
temp$MACEg[is.na(temp$MACEg)]<-"M-miss"
temp<-spread(temp,MACEg,MACEg)
#reformatting R81T genotypes
temp$R81Tg<-temp$R81T
levels(temp$R81Tg)<-c("Neo-RR","Neo-RS","Neo-SS")
temp$R81Tg<-as.character(temp$R81Tg)
temp$R81Tg[is.na(temp$R81Tg)]<-"Neo-miss"
temp<-spread(temp,R81Tg,R81Tg)

#turning the resistance data into a 0-1 matrix 
temp<-cbind(temp,ifelse(is.na(temp[,c(56:70)]),0,1))
temp<-temp[,-c(56:70)]

#reordering the levels of the host
temp$host<-factor(temp$host,levels(temp$host)[c(4,2,5,3,1)])
setorder(temp,host,year,sampling_date)
names(table(temp$host))

#preparing the structure plot
poptiquet<-c("Peach","Oilseed rape","Tobacco","Other\nCrops","Aerial Trap")
effpop<-as.numeric(table(temp$host))
strK3<-t(temp[,c("K3_Q1","K3_Q2","K3_Q3")])
strK4<-t(temp[,c("K4_Q1","K4_Q2","K4_Q3","K4_Q4")])
strK5<-t(temp[,c("K5_Q1","K5_Q2","K5_Q3","K5_Q4","K5_Q5")])
strKDR<-t(temp[,c("K-RR","K-RS","K-SS","K-miss")])
strsKDR<-t(temp[,c("sK-RR","sK-RS","sK-SS","sK-miss")])
strMACE<-t(temp[,c("M-RS","M-SS","M-miss")])
strNEO<-t(temp[,c("Neo-RR","Neo-RS","Neo-SS","Neo-miss")])

#the plot for the different K values
layout(matrix(c(1,1,1,
                2,2,2,
                3,3,3,
                4,5,6,7),13,1,byrow=TRUE))
op<-par(mar=c(0.1,1.1,0.1,0),oma=c(5.1,3,1,0))

coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
structplot(strK3,coloor,effpop,poptiquet,spacepop=4,
           leg_y="K=3",cexy=1,mef=c(0,1,1,0,1),colbord=NA)
structplot(strK4,coloor[c(1,3,2,4,5)],effpop,poptiquet,spacepop=4,
           leg_y="K=4",cexy=1,mef=c(0,1,1,0,1),colbord=NA)
structplot(strK5,coloor[c(1,3,2,4,5)],effpop,poptiquet,spacepop=4,
           leg_y="K=5",cexy=1,mef=c(0,1,1,0,1),colbord=NA)

coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],brewer.pal(9,"Greens")[5],"grey80")
structplot(strKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="KDR",cexy=1,mef=c(0,1,1,0,1),colbord=NA)
structplot(strsKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="sKDR",cexy=1,mef=c(0,1,1,0,1),colbord=NA)
structplot(strMACE,coloor[c(2,3,4)],effpop,poptiquet,spacepop=4,
           leg_y="MACE",cexy=1,mef=c(0,1,1,0,1),colbord=NA)
structplot(strNEO,coloor,effpop,poptiquet,spacepop=4,
           leg_y="R81T",cexy=1,mef=c(0,1,1,1,1),colbord=NA,
           cexpop=1.5,distxax=0.3)

par(op)

#export to .pdf 12 x 7 inches


##############################################################################/
#END
##############################################################################/