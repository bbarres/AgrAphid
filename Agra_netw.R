##############################################################################/
##############################################################################/
#Clonal/MLG network analyses
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the dataset for genind importation####
##############################################################################/

#we reorganize the levels of the host_corrected column, because the 
#alphabetical order doesn't fit our needs
datAgra$host_corrected<-factor(datAgra$host_corrected,
                               levels=c("peach","oilseed_rape","tobacco",
                                        "other_crops","Aerial_trap",
                                        "several_hosts"))

#pick a set of color
coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")



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
#Aerial samples####
##############################################################################/

#grouped by Cluster K=3
AerTrap_ClustK3<-datAgra[datAgra$host=="Aerial_trap",]
AerTrap_ClustK3$Clust_K3<-as.character(AerTrap_ClustK3$Clust_K3)
AerTrap_ClustK3$Clust_K3[is.na(AerTrap_ClustK3$Clust_K3)]<-"undef"
AerTrap_ClustK3<-drop.levels(AerTrap_ClustK3)
#converting to genind object
dataNetwo<-df2genind(AerTrap_ClustK3[,c("MP_27","MP_39","MP_44","MP_5",
                                        "MP_7","MP_23","MP_45","MP_28",
                                        "MP_9","MP_13","MP_2","MP_38",
                                        "MP_4","MP_46")],
                        ncode=3,
                        ind.names=AerTrap_ClustK3$indiv_ID, 
                        pop=AerTrap_ClustK3$Clust_K3,
                        ploidy=2,NA.char="999")

set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,provesti.dist(dataNetwo),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         dist.asd(genind2loci(dataNetwo)),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         diss.dist(dataNetwo),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4)],
               scale.leg=FALSE,
               mlg=TRUE)


#grouped by Cluster K=5
AerTrap_ClustK5<-datAgra[datAgra$host=="Aerial_trap",]
AerTrap_ClustK5$Clust_K5<-as.character(AerTrap_ClustK5$Clust_K5)
AerTrap_ClustK5$Clust_K5[is.na(AerTrap_ClustK5$Clust_K5)]<-"undef"
AerTrap_ClustK5<-drop.levels(AerTrap_ClustK5)
#converting to genind object
dataNetwo<-df2genind(AerTrap_ClustK5[,c("MP_27","MP_39","MP_44","MP_5",
                                        "MP_7","MP_23","MP_45","MP_28",
                                        "MP_9","MP_13","MP_2","MP_38",
                                        "MP_4","MP_46")],
                        ncode=3,
                        ind.names=AerTrap_ClustK5$indiv_ID, 
                        pop=AerTrap_ClustK5$Clust_K5,
                        ploidy=2,NA.char="999")


set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,provesti.dist(dataNetwo),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4,6,5)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4,6,5)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         dist.asd(genind2loci(dataNetwo)),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4,6,5)],
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         diss.dist(dataNetwo),
                         include.ties=TRUE),
               nodescale=2,
               palette=coloor[c(3,2,1,4,6,5)],
               scale.leg=FALSE,
               mlg=TRUE)


#only repeated clones grouped by year
#grouped by Cluster K=5
AerTrap_rep<-datAgra[datAgra$host=="Aerial_trap" & datAgra$repeated==1,]
AerTrap_ClustK5$Clust_K5<-as.character(AerTrap_ClustK5$Clust_K5)
AerTrap_ClustK5$Clust_K5[is.na(AerTrap_ClustK5$Clust_K5)]<-"undef"
AerTrap_rep<-drop.levels(AerTrap_rep)
#converting to genind object
dataNetwo<-df2genind(AerTrap_rep[,c("MP_27","MP_39","MP_44","MP_5",
                                        "MP_7","MP_23","MP_45","MP_28",
                                        "MP_9","MP_13","MP_2","MP_38",
                                        "MP_4","MP_46")],
                     ncode=3,
                     ind.names=AerTrap_rep$indiv_ID, 
                     pop=AerTrap_rep$year,
                     ploidy=2,NA.char="999")

set.seed(333)
plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         bruvo.dist(dataNetwo,replen=rep(1,14)),
                         include.ties=TRUE),
               nodescale=10,
               palette=brewer.pal(8,"Dark2"),
               scale.leg=FALSE,
               mlg=TRUE)

plot_poppr_msn(dataNetwo,
               poppr.msn(dataNetwo,
                         diss.dist(dataNetwo),
                         include.ties=FALSE),
               nodescale=4,
               palette=brewer.pal(8,"Dark2"),
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


##############################################################################/
#END
##############################################################################/