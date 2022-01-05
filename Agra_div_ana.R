##############################################################################/
##############################################################################/
#MLG and genetic diversity indices computation
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#MLG Diversity indices by host population####
##############################################################################/

#To compute the diversity we need a presence/absence of MLG table
#for each host population
datHost<-datAgra[datAgra$host!="Aerial_trap",]
datHost<-drop.levels(datHost)
hostocc<-table(datHost$host,datHost$MLG_ID)
Host_nb_samples<-rowSums(hostocc)
Host_nb_MLG<-specnumber(hostocc)
Host_GsurN<-specnumber(hostocc)/rowSums(hostocc)
Host_MLG_richness<-rarefy(hostocc,min(rowSums(hostocc)))
Host_simpson_div<-vegan::diversity(hostocc,index="simpson")
Host_pielou_even<-vegan::diversity(hostocc)/log(specnumber(hostocc))


##############################################################################/
#Genetic Diversity indices by host population####
##############################################################################/

#converting data to a genind format for the full data set
compdiv<-datHost #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=3,pop=compdiv$host,ploidy=2,
                  NA.char=c("999"),ind.names=as.character(compdiv$indiv_ID))
#Allelic richness for the full temporal data set
Host_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#Private Allelic richness for the full temporal data set
Host_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#Nei heterozygosity for the full temporal data set
Host_HetNei<-Hs(COMPDI)


##############################################################################/
#building the final table for host population diversity####
##############################################################################/

divHosttable<-cbind(Host_nb_samples,Host_nb_MLG,Host_GsurN,Host_MLG_richness,
                   Host_simpson_div,Host_pielou_even,Host_Ar,Host_PrivAr,
                   Host_HetNei)
write.table(divHosttable,file="output/divHosttable.txt",
            quote=FALSE,sep="\t")


##############################################################################/
#MLG Diversity indices by genetic cluster####
##############################################################################/

#To compute the diversity we need a presence/absence of MLG table
#for each genetic cluster
Aerial_cor<-Aerial[!is.na(Aerial$Clust_K4),]
clustocc<-table(Aerial_cor$Clust_K4,Aerial_cor$MLG_ID)
Clust_nb_samples<-rowSums(clustocc)
Clust_nb_MLG<-specnumber(clustocc)
Clust_GsurN<-specnumber(clustocc)/rowSums(clustocc)
Clust_MLG_richness<-rarefy(clustocc,min(rowSums(clustocc)))
Clust_simpson_div<-vegan::diversity(clustocc,index="simpson")
Clust_pielou_even<-vegan::diversity(clustocc)/log(specnumber(clustocc))


##############################################################################/
#Genetic Diversity indices by Cluster on the complete data set####
##############################################################################/

#converting data to a genind format for the full data set
compdiv<-Aerial_cor #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=3,pop=compdiv$Clust_K4,ploidy=2,
                  NA.char=c("999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-Aerial_CC #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=3,pop=compdiv$Clust_K4,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999"))

#Allelic richness for the full temporal data set
Clust_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#Private Allelic richness for the full temporal data set
Clust_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#Nei heterozygosity for the full temporal data set
Clust_HetNei<-Hs(COMPDI)


##############################################################################/
#building the final table for genetic cluster diversity####
##############################################################################/

divClusttable<-cbind(Clust_nb_samples,Clust_nb_MLG,Clust_GsurN,
                     Clust_MLG_richness,Clust_simpson_div,Clust_pielou_even,
                     Clust_Ar,Clust_PrivAr,Clust_HetNei)
write.table(divClusttable,file="output/divClusttable.txt",
            quote=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/