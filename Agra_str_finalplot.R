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
#because missing data will screw up the ordering we complete one missing
#data for year
temp$year[235]<-2014
setorder(temp,host,year,sampling_date,K4_Q1,-K4_Q4,-K4_Q2)
names(table(temp$host))

#creating the offset for the box per year
temp2<-as.data.frame(table(temp$year,temp$host))
temp2<-temp2[temp2$Freq!=0,]
temp2$cumu<-cumsum(temp2$Freq)
temp2$Var2b<-temp2$Var2[c(1,1:(length(temp2$Var2)-1))]
temp2$decal<-cumsum(ifelse(temp2$Var2==temp2$Var2b,0,4))

#preparing the structure plot
poptiquet<-c("Peach","Oilseed\nrape","Tobacco","Other\nCrops","Aerial Trap")
effpop<-as.numeric(table(temp$host))
strK3<-t(temp[,c("K3_Q1","K3_Q2","K3_Q3")])
strK4<-t(temp[,c("K4_Q1","K4_Q2","K4_Q3","K4_Q4")])
#strK5<-t(temp[,c("K5_Q1","K5_Q2","K5_Q3","K5_Q4","K5_Q5")])
strKDR<-t(temp[,c("K-RR","K-RS","K-SS","K-miss")])
strsKDR<-t(temp[,c("sK-RR","sK-RS","sK-SS","sK-miss")])
strMACE<-t(temp[,c("M-RS","M-SS","M-miss")])
strNEO<-t(temp[,c("Neo-RR","Neo-RS","Neo-SS","Neo-miss")])

#the plot for the different K values
layout(matrix(c(1,1,1,1,
                2,2,2,2,
                3,4,5,6),12,1,byrow=TRUE))
op<-par(mar=c(0.1,1.1,0.1,0),oma=c(4.1,4,3.5,0),font=2)

coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
structplot(strK3,coloor,effpop,poptiquet,spacepop=4,
           leg_y="K=3",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("Genetic\ncluster\nK=3",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)
#adding sampling years for aerial trap "population"
axis(3,at=c(c(0,temp2$cumu[1:18])+temp2$decal+
              (temp2$cumu-c(0,temp2$cumu[1:18]))/2)[12:19],
     labels=FALSE,pos=1,lwd.ticks=2,lwd=0)
text(c(c(0,temp2$cumu[1:18])+temp2$decal+
               (temp2$cumu-c(0,temp2$cumu[1:18]))/2)[12:19]-5,
     rep(par("usr")[4]+0.1,19)[12:19],
     labels=temp2$Var1[12:19],
     srt=45,xpd=NA,pos=4,cex=1.4,font=2)

structplot(strK4,coloor[c(1,3,2,4,5)],effpop,poptiquet,spacepop=4,
           leg_y="K=4",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("Genetic\ncluster\nK=4",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)

coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],brewer.pal(9,"Greens")[5],"grey80")
structplot(strKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="KDR",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("kdr",side=2,line=-2,cex=1.2,las=1,font=4,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)
structplot(strsKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="sKDR",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("skdr",side=2,line=-2,cex=1.2,las=1,font=4,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)
structplot(strMACE,coloor[c(2,3,4)],effpop,poptiquet,spacepop=4,
           leg_y="MACE",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("MACE",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)
structplot(strNEO,coloor,effpop,poptiquet,spacepop=4,
           leg_y="R81T",cexy=1,mef=c(0,1,1,1,0),colbord=NA,
           cexpop=1.5,distxax=0.3)
mtext("R81T",side=2,line=-2,cex=1.2,las=1,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal)[12:19],
     rep(0,length(temp2$cumu))[12:19],
     c(temp2$cumu+temp2$decal)[12:19],
     rep(1,length(temp2$cumu))[12:19],
     lwd=2)

par(op)

#export to .pdf 12 x 5 inches


##############################################################################/
#END
##############################################################################/

#the plot for the different K values
layout(matrix(c(1,1,1,1,
                2,2,2,2,
                3,4,5,6),12,1,byrow=TRUE))
op<-par(mar=c(0.1,1.1,0.1,0),oma=c(4.1,4,3.5,0),font=2)

coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")
structplot(strK3,coloor,effpop,poptiquet,spacepop=4,
           leg_y="K=3",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("Genetic\ncluster\nK=3",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)
#adding sampling years for aerial trap "population"
axis(3,at=c(c(0,temp2$cumu[1:18])+temp2$decal+
                    (temp2$cumu-c(0,temp2$cumu[1:18]))/2),
     labels=FALSE,pos=1,lwd.ticks=2,lwd=0)
text(c(c(0,temp2$cumu[1:18])+temp2$decal+c(0,0,0,0,0,0,-1,1,2,0,0,0,0,0,0,0,0,0,0)+
               (temp2$cumu-c(0,temp2$cumu[1:18]))/2)-2,
     rep(par("usr")[4]+0.05,19),
     labels=temp2$Var1,
     srt=90,xpd=NA,pos=4,cex=1.1,font=2)

structplot(strK4,coloor[c(1,3,2,4,5)],effpop,poptiquet,spacepop=4,
           leg_y="K=4",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("Genetic\ncluster\nK=4",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)

coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],brewer.pal(9,"Greens")[5],"grey80")
structplot(strKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="KDR",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("kdr",side=2,line=-2,cex=1.2,las=1,font=4,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)
structplot(strsKDR,coloor,effpop,poptiquet,spacepop=4,
           leg_y="sKDR",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("skdr",side=2,line=-2,cex=1.2,las=1,font=4,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)
structplot(strMACE,coloor[c(2,3,4)],effpop,poptiquet,spacepop=4,
           leg_y="MACE",cexy=1,mef=c(0,1,1,0,0),colbord=NA)
mtext("MACE",side=2,line=-2,cex=1.2,las=1,font=2,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)
structplot(strNEO,coloor,effpop,poptiquet,spacepop=4,
           leg_y="R81T",cexy=1,mef=c(0,1,1,1,0),colbord=NA,
           cexpop=1.5,distxax=0.3)
mtext("R81T",side=2,line=-2,cex=1.2,las=1,adj=1)
rect(c(c(0,temp2$cumu)[1:length(temp2$cumu)]+temp2$decal),
     rep(0,length(temp2$cumu)),
     c(temp2$cumu+temp2$decal),
     rep(1,length(temp2$cumu)),
     lwd=2)

par(op)

#export to .pdf 12 x 5 inches



c(0,0,0,0,0,0,-1,1,2,0,0,0,0,0,0,0,0,0,0)