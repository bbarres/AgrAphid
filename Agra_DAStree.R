##############################################################################/
##############################################################################/
#DAS tree on the clone corrected data
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the data set for network analysis####
##############################################################################/

#preparing the dataset
temp<-as.data.table(Aerial_CC)
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

#grouped by Cluster K=4
AerTrap_ClustK4<-temp
AerTrap_ClustK4$Clust_K4<-as.character(AerTrap_ClustK4$Clust_K4)
AerTrap_ClustK4$Clust_K4[is.na(AerTrap_ClustK4$Clust_K4)]<-"undef"
AerTrap_ClustK4<-drop.levels(AerTrap_ClustK4)

#converting to genind object
datArb<-df2genind(AerTrap_ClustK4[,c("cor_MP_27","cor_MP_39","cor_MP_44",
                                     "cor_MP_5","cor_MP_7","cor_MP_23",
                                     "cor_MP_45","cor_MP_28","cor_MP_9",
                                     "cor_MP_13","cor_MP_2","cor_MP_38",
                                     "cor_MP_4","cor_MP_46")],
                     ncode=3,
                     ind.names=AerTrap_ClustK4$indiv_ID, 
                     pop=AerTrap_ClustK4$Clust_K4,
                     ploidy=2,NA.char="999")
datArb.mean<-missingno(datArb,type="mean")


##############################################################################/
#Figure S5: Aerial samples neighbour joining tree####
##############################################################################/

#with a dissimilarity distance computed with an ape function
distDASArb<-diss.dist(datArb)
treenj<-bionj(distDASArb)
#pick a set of color
coloor<-c("royalblue4","firebrick","khaki2","chartreuse4","grey80")
op<-par(mfrow=c(2,3))
plot(treenj,type="radial",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treenj,type="cladogram",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treenj,type="fan",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treenj,type="unrooted",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treenj,type="phylogram",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
par(op)

#final DAS plot
op<-par(mar=c(0.1,0.1,0.1,0.1))
plot(treenj,type="unrooted",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
par(op)
#export to .pdf 6 x 8 inches


##############################################################################/
#END
##############################################################################/