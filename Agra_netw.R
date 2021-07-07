##############################################################################/
##############################################################################/
#Clonal/MLG network analyses
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the dataset for network analysis####
##############################################################################/

#preparing the dataset
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
#Aerial samples####
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
coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],
          brewer.pal(9,"Greens")[5],"grey80")[c(2,3,4,1)]
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
coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],
          brewer.pal(9,"Greens")[5],"grey80")[c(3,4,2,1)]
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
coloor<-c(brewer.pal(9,"YlOrRd")[c(8,6)],
          brewer.pal(9,"Greens")[5],"grey80")[c(3,4,2)]
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

par(op)

#export to .pdf 13 x 13 inches


##############################################################################/
#END
##############################################################################/







##############################################################################/
#All samples####
##############################################################################/

#grouped by host
All_repet<-datAgra
All_repet$Clust_K5<-as.character(All_repet$Clust_K5)
All_repet$Clust_K5[is.na(All_repet$Clust_K5)]<-"undef"
All_repet<-drop.levels(All_repet)
#converting to genind object
dataNetwo<-df2genind(All_repet[,c("MP_27","MP_39","MP_44","MP_5",
                                  "MP_7","MP_23","MP_45","MP_28",
                                  "MP_9","MP_13","MP_2","MP_38",
                                  "MP_4","MP_46")],
                     ncode=3,
                     ind.names=All_repet$indiv_ID, 
                     pop=All_repet$Clust_K5,
                     ploidy=2,NA.char="999")

set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=4,
               palette=coloor[c(1,6,5,3,4,2)],
               scale.leg=FALSE,
               mlg=TRUE,
               label.color="invisible")


##############################################################################/
#All repeated samples####
##############################################################################/

#grouped by host
All_repet<-datAgra[datAgra$repeated==1,]
All_repet$Clust_K3<-as.character(All_repet$Clust_K3)
All_repet$Clust_K3[is.na(All_repet$Clust_K3)]<-"undef"
All_repet<-drop.levels(All_repet)
#converting to genind object
dataNetwo<-df2genind(All_repet[,c("MP_27","MP_39","MP_44","MP_5",
                                  "MP_7","MP_23","MP_45","MP_28",
                                  "MP_9","MP_13","MP_2","MP_38",
                                  "MP_4","MP_46")],
                     ncode=3,
                     ind.names=All_repet$indiv_ID, 
                     pop=All_repet$host,
                     ploidy=2,NA.char="999")

set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=4,
               palette=brewer.pal(5,"Set2"),
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         diss.dist(dataNetwo),
                         include.ties=FALSE),
               nodescale=4,
               palette=brewer.pal(5,"Set2"),
               scale.leg=FALSE,
               mlg=TRUE)

plot(spantree(diss.dist(dataNetwo)))


#group by K=3 STUCTURE genetic clusters
#converting to genind object
dataNetwo<-df2genind(All_repet[,c("MP_27","MP_39","MP_44","MP_5",
                                  "MP_7","MP_23","MP_45","MP_28",
                                  "MP_9","MP_13","MP_2","MP_38",
                                  "MP_4","MP_46")],
                     ncode=3,
                     ind.names=All_repet$indiv_ID, 
                     pop=All_repet$Clust_K3,
                     ploidy=2,NA.char="999")

set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=4,
               palette=coloor[c(1,3,4,2)],
               scale.leg=FALSE,
               mlg=TRUE)


##############################################################################/
#Oilseed rape samples####
##############################################################################/

#grouped by Cluster K=3
oil_ClustK3<-datAgra[datAgra$host=="oilseed_rape",]
oil_ClustK3$Clust_K3<-as.character(oil_ClustK3$Clust_K3)
oil_ClustK3$Clust_K3[is.na(oil_ClustK3$Clust_K3)]<-"undef"
oil_ClustK3<-drop.levels(oil_ClustK3)
dataNetwo<-df2genind(oil_ClustK3[,c("MP_27","MP_39","MP_44","MP_5",
                                    "MP_7","MP_23","MP_45","MP_28",
                                    "MP_9","MP_13","MP_2","MP_38",
                                    "MP_4","MP_46")],
                     ncode=3,
                     ind.names=oil_ClustK3$indiv_ID, 
                     pop=oil_ClustK3$Clust_K3,
                     ploidy=2,NA.char="999")

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=FALSE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)

#grouped by KDR genotype
oil_KDR<-datAgra[datAgra$host=="oilseed_rape",]
oil_KDR$KDR<-as.character(oil_KDR$KDR)
oil_KDR$KDR[is.na(oil_KDR$KDR)]<-"miss"
oil_KDR<-drop.levels(oil_KDR)
dataNetwo<-df2genind(oil_KDR[,c("MP_27","MP_39","MP_44","MP_5",
                                "MP_7","MP_23","MP_45","MP_28",
                                "MP_9","MP_13","MP_2","MP_38",
                                "MP_4","MP_46")],
                     ncode=3,
                     ind.names=oil_KDR$indiv_ID, 
                     pop=oil_KDR$KDR,
                     ploidy=2,NA.char="999")

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=FALSE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)

#grouped by sKDR genotype
oil_sKDR<-datAgra[datAgra$host=="oilseed_rape",]
oil_sKDR$sKDR<-as.character(oil_sKDR$sKDR)
oil_sKDR$sKDR[is.na(oil_sKDR$sKDR)]<-"miss"
oil_sKDR<-drop.levels(oil_sKDR)
dataNetwo<-df2genind(oil_sKDR[,c("MP_27","MP_39","MP_44","MP_5",
                                 "MP_7","MP_23","MP_45","MP_28",
                                 "MP_9","MP_13","MP_2","MP_38",
                                 "MP_4","MP_46")],
                     ncode=3,
                     ind.names=oil_sKDR$indiv_ID, 
                     pop=oil_sKDR$sKDR,
                     ploidy=2,NA.char="999")

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=FALSE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)
