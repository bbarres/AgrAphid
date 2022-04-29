##############################################################################/
##############################################################################/
#Clonal/MLG network analyses
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the data set for network analysis####
##############################################################################/

#preparing the data set
temp<-as.data.table(Aerial)
#reformatting KDR genotypes
temp$KDRg<-temp$KDR
levels(temp$KDRg)<-c("K-RR","K-RS","K-SS")
temp$KDRg<-as.character(temp$KDRg)
temp$KDRg[is.na(temp$KDRg)]<-"K-miss"
#reformatting sKDR genotypes
temp$sKDRg<-temp$sKDR
levels(temp$sKDRg)<-c("sK-RR","sK-RS","sK-RS","sK-SS","sK-RS",
                      "sK-RS","sK-RR")
temp$sKDRg<-as.character(temp$sKDRg)
temp$sKDRg[is.na(temp$sKDRg)]<-"sK-miss"
#reformatting MACE genotypes
temp$MACEg<-temp$MACE
levels(temp$MACEg)<-c("M-SS","M-RS")
temp$MACEg<-as.character(temp$MACEg)
temp$MACEg[is.na(temp$MACEg)]<-"M-miss"
#reformatting R81T genotypes
temp$R81Tg<-temp$R81T
levels(temp$R81Tg)<-c("Neo-RR","Neo-RS","Neo-SS")
temp$R81Tg<-as.character(temp$R81Tg)
temp$R81Tg[is.na(temp$R81Tg)]<-"Neo-miss"


##############################################################################/
#Figure 2: network of the MLG from the aerial trap####
##############################################################################/

#grouped by Cluster K=4
AerTrap_ClustK4<-temp
AerTrap_ClustK4$Clust_K4<-as.character(AerTrap_ClustK4$Clust_K4)
AerTrap_ClustK4$Clust_K4[is.na(AerTrap_ClustK4$Clust_K4)]<-"undef"
AerTrap_ClustK4<-drop.levels(AerTrap_ClustK4)

op<-par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(0,0,0,0))

#converting to genind object
dataNetwo<-df2genind(AerTrap_ClustK4[,c("cor_MP_27","cor_MP_39","cor_MP_44",
                                        "cor_MP_5","cor_MP_7","cor_MP_23",
                                        "cor_MP_45","cor_MP_28","cor_MP_9",
                                        "cor_MP_13","cor_MP_2","cor_MP_38",
                                        "cor_MP_4","cor_MP_46")],
                        ncode=3,
                        ind.names=AerTrap_ClustK4$indiv_ID, 
                        pop=AerTrap_ClustK4$Clust_K4,
                        ploidy=2,NA.char="999")
#pick a set of color
coloor<-c("royalblue4","firebrick","khaki2","chartreuse4","grey80")
#plotting the network
truc<-poppr.msn(dataNetwo,showplot=FALSE,
                diss.dist(dataNetwo),
                include.ties=TRUE)
#5,9,13,17,24,30,320
set.seed(320)
plot_poppr_msn(dataNetwo,
               truc,cex.main=3,
               nodescale=8,pop.leg=FALSE,size.leg=FALSE,
               palette=coloor,
               scale.leg=FALSE,wscale=FALSE,
               mlg=FALSE,
               inds="",nodelab=200)
par(mfg=c(1,1))
mtext("Genetic Cluster (K=4)",side=3,cex=1.8,line=-2.5,font=2)
mtext("  A.",side=3,cex=3,line=-3,font=2,las=0,adj=0)

par(mfg=c(1,2))
#grouped by KDR
dataNetwo<-df2genind(AerTrap_ClustK4[,c("cor_MP_27","cor_MP_39","cor_MP_44",
                                        "cor_MP_5","cor_MP_7","cor_MP_23",
                                        "cor_MP_45","cor_MP_28","cor_MP_9",
                                        "cor_MP_13","cor_MP_2","cor_MP_38",
                                        "cor_MP_4","cor_MP_46")],
                     ncode=3,
                     ind.names=AerTrap_ClustK4$indiv_ID, 
                     pop=as.character(AerTrap_ClustK4$KDRg),
                     ploidy=2,NA.char="999")
#pick a set of color
coloor<-coloor<-c(brewer.pal(9,"Purples")[7],
                  brewer.pal(9,"Oranges")[6],
                  brewer.pal(9,"Greens")[5],
                  "grey80")[c(2,3,4,1)]
#plotting the network
truc<-poppr.msn(dataNetwo,showplot=FALSE,
                diss.dist(dataNetwo),
                include.ties=TRUE)
set.seed(320)
plot_poppr_msn(dataNetwo,
               truc,
               nodescale=8,pop.leg=FALSE,size.leg=FALSE,
               palette=coloor,
               scale.leg=FALSE,wscale=FALSE,
               mlg=FALSE,
               inds="",nodelab=200)
par(mfg=c(1,2))
mtext("kdr",side=3,cex=1.8,line=-2.5,font=4)
mtext("  B.",side=3,cex=3,line=-3,font=2,las=0,adj=0)

par(mfg=c(2,1))
#grouped by sKDR
dataNetwo<-df2genind(AerTrap_ClustK4[,c("cor_MP_27","cor_MP_39","cor_MP_44",
                                        "cor_MP_5","cor_MP_7","cor_MP_23",
                                        "cor_MP_45","cor_MP_28","cor_MP_9",
                                        "cor_MP_13","cor_MP_2","cor_MP_38",
                                        "cor_MP_4","cor_MP_46")],
                     ncode=3,
                     ind.names=AerTrap_ClustK4$indiv_ID, 
                     pop=as.character(AerTrap_ClustK4$sKDRg),
                     ploidy=2,NA.char="999")
#pick a set of color
coloor<-coloor<-c(brewer.pal(9,"Purples")[7],
                  brewer.pal(9,"Oranges")[6],
                  brewer.pal(9,"Greens")[5],
                  "grey80")[c(3,4,2,1)]
#plotting the network
truc<-poppr.msn(dataNetwo,showplot=FALSE,
                diss.dist(dataNetwo),
                include.ties=TRUE)
set.seed(320)
plot_poppr_msn(dataNetwo,
               truc,
               nodescale=8,pop.leg=FALSE,size.leg=FALSE,
               palette=coloor,
               scale.leg=FALSE,wscale=FALSE,
               mlg=FALSE,
               inds="",nodelab=200)
par(mfg=c(2,1))
mtext("skdr",side=3,cex=1.8,line=-2.5,font=4)
mtext("  C.",side=3,cex=3,line=-3,font=2,las=0,adj=0)

par(mfg=c(2,2))
#grouped by MACE
dataNetwo<-df2genind(AerTrap_ClustK4[,c("cor_MP_27","cor_MP_39","cor_MP_44",
                                        "cor_MP_5","cor_MP_7","cor_MP_23",
                                        "cor_MP_45","cor_MP_28","cor_MP_9",
                                        "cor_MP_13","cor_MP_2","cor_MP_38",
                                        "cor_MP_4","cor_MP_46")],
                     ncode=3,
                     ind.names=AerTrap_ClustK4$indiv_ID, 
                     pop=as.character(AerTrap_ClustK4$MACEg),
                     ploidy=2,NA.char="999")
#pick a set of color
coloor<-coloor<-c(brewer.pal(9,"Purples")[7],
                  brewer.pal(9,"Oranges")[6],
                  brewer.pal(9,"Greens")[5],
                  "grey80")[c(3,4,2)]
#plotting the network
truc<-poppr.msn(dataNetwo,showplot=FALSE,
                diss.dist(dataNetwo),
                include.ties=TRUE)
set.seed(320)
plot_poppr_msn(dataNetwo,
               truc,
               nodescale=8,pop.leg=FALSE,size.leg=FALSE,
               palette=coloor,
               scale.leg=FALSE,wscale=FALSE,
               mlg=FALSE,
               inds="",nodelab=200)
par(mfg=c(2,2))
mtext("MACE",side=3,cex=1.8,line=-2.5,font=2)
mtext("  D.",side=3,cex=3,line=-3,font=2,las=0,adj=0)

par(op)

#export to .pdf 13 x 13 inches


##############################################################################/
#END
##############################################################################/