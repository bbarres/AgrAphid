##############################################################################/
##############################################################################/
#DAS tree on the clone corrected data
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the dataset for network analysis####
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
#Aerial samples####
##############################################################################/

#Shared allele distance computed using population1.2.32
mat<-read.table("data/AgraDAS_mat.txt",header=FALSE,
                check.names=TRUE,row.names=1)
mat<-as.dist(mat)
treeDAS<-bionj(mat)
labord<-as.data.frame(attr(mat,"Labels"))
colnames(labord)="labels"
labord<-merge(labord,AerTrap_ClustK4,by.x="labels",by.y="indiv_ID",sort=FALSE)
#pick a set of color
coloor<-c("khaki2","royalblue4","firebrick",
          "grey80","chartreuse4")
op<-par(mfrow=c(2,3))
plot(treeDAS,type="radial",show.tip=FALSE,main="Radial")
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treeDAS,type="cladogram",show.tip=FALSE,main="Cladogram")
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treeDAS,type="fan",show.tip=FALSE,main="Fan")
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treeDAS,type="unrooted",show.tip=FALSE,main="Unrooted")
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
plot(treeDAS,type="phylogram",show.tip=FALSE,main="Phylogram")
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
par(op)
#export to .pdf 18 x 13 inches

#final DAS plot
op<-par(mar=c(0.1,0.1,0.1,0.1))
plot(treeDAS,type="unrooted",show.tip=FALSE)
tiplabels(pch=20,col=coloor[datArb@pop],cex=2)
par(op)
#export to .pdf 6 x 10 inches


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


#en premier lieu, essayons de représenter les individus sous forme d'arbre 
#basé sur une matrice de distance génétique pour voir comment s'organisent 
#les individus les uns par rapport aux autres

#dans cette exemple on importe la matrice de distance à partir de la sortie 
#de calcul de distance entre individus de TASSEL. La distance calculée se 
#rapproche d'une DAS (Shared Allele Distance, à revérifier qd même)
mat<-read.table("RT_desc_tot_mat.txt",header=T,check.names=TRUE, row.names=1)
mat<-as.dist(mat)

#en plus du fichier de distance, on a besoin d'un fichier qui fait le lien 
#entre les individus et leur appartenance à une famille, on utilise pour 
#cela des colonnes réalisées à partir du fichier RT_verif_fam. Attention à 
#garder le même ordre que dans le fichier de matrice de distance
#D'abord exemple de différentes méthodes de reconstruction des arbres
trenj<-nj(mat)
trebionj<-bionj(mat)
#tremst<-mst(mat)

#ensuite il faut gérer les labels
lab<-read.table("RT_desc_tot_etiq.txt",skip=1)
#en principe pas besoin de remettre dans l'ordre si on a pris le même ordre 
#que pour la matrice de distance, mais pour bionj, tout est mélangé, donc 
# on se débrouille pour que l'ordre soit respecter dans tous les cas
treOrd<-as.matrix(trebionj$tip.label)
labOrd<-merge(treOrd,lab,by="V1",sort=F,suffixes=c(".x",".y"))
colo<-as.numeric(labOrd$V2)
palette(rainbow(length(levels(as.factor(colo)))))

#et puis on réalise le plot ensuite, voir plus bas pour d'autre forme 
#d'arbre et les impression des fichiers en pdf
plot(trebionj,type="unr",show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 2)
tiplabels(trebionj$tip.label, cex = 0.1,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")

#impression pdf d'un arbre non-enraciné
plot(trebionj,type="unr",show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(labOrd$V1,labOrd$V3,labOrd$V4), 
          cex = 0.25,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")

#impression pdf d'un arbre en rond
plot(trebionj,type="radial",show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(labOrd$V1,labOrd$V3,labOrd$V4), 
          cex = 0.5,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")

#impression pdf d'un arbre en arbre
plot(trebionj,show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(labOrd$V1,labOrd$V3,labOrd$V4), 
          cex = 0.5,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")




#ce qui est décrit juste au-dessus pour ranger les individus ne marchent 
#pas bizarement dans le cas de RT_tot. Du coup on fait le travail sur 
#excell pour que l'ordre de trebionj colle avec le fichier RT_etiq.txt
#la raison c'est l'impossibilité de rentrer la matrice de distance avec 
#des noms d'individus identiques pour remédier à ça voilà la démarche
#quand il y a une répétition du nom d'individu, il faut utiliser ça et enlever 
#la première colonne dans la matrice de distance avant l'importation
mat<-read.table("totdesc_dist.mat.txt",header=T,check.names=TRUE)#, row.names=1)
mat<-as.dist(mat)
trebionj<-bionj(mat)
write.table(trebionj$tip.label,file="ordretrebionj.txt")
#ensuite on fait recherche verticale sur le fichier d'étiquette pour remettre
#dans l'ordre de "ordretrebionj.txt", puis on importe
lab<-read.table("etiq_totdesc.txt",header=T)
colo<-as.numeric(lab$population)
palette(rainbow(length(levels(as.factor(colo)))))

#impression pdf d'un arbre non-enraciné
plot(trebionj,type="unr",show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(lab$nom_Ind,lab$population,lab$father), cex = 0.25,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")

#impression pdf d'un arbre en rond
plot(trebionj,type="radial",show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(lab$population,lab$father), cex = 0.5,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")

#impression pdf d'un arbre en arbre
plot(trebionj,show.tip=FALSE)
tiplabels(pch = 20, col = colo, cex = 4)
tiplabels(paste(trebionj$tip.label,lab$population,lab$father), cex = 0.5,frame="none",bg=colo)
title("Neighbour-Joining tree of oaks of the natural stands")


