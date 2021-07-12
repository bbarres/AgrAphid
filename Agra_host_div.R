##############################################################################/
##############################################################################/
#By cluster diversity indices computation
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#MLG Diversity indices by genetic cluster####
##############################################################################/

#To compute the diversity we need a presence/absence of MLG table
#for each genetic cluster
Aerial_cor<-Aerial[!is.na(Aerial$Clust_K4),]
hostocc<-table(Aerial_cor$Clust_K4,Aerial_cor$MLG_ID)
Host_nb_samples<-rowSums(hostocc)
Host_nb_MLG<-specnumber(hostocc)
Host_GsurN<-specnumber(hostocc)/rowSums(hostocc)
Host_MLG_richness<-rarefy(hostocc,min(rowSums(hostocc)))
Host_simpson_div<-vegan::diversity(hostocc,index="simpson")
Host_pielou_even<-vegan::diversity(hostocc)/log(specnumber(hostocc))


##############################################################################/
#Genetic Diversity indices by Cluster on the complete dataset####
##############################################################################/

#converting data to a genind format for the full dataset
compdiv<-Aerial_cor #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=3,pop=compdiv$Clust_K4,ploidy=2,
                  NA.char=c("999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-Aerial_CC #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=3,pop=compdiv$Clust_K4,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999"))

#Allelic richness for the full temporal dataset
Host_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#Private Allelic richness for the full temporal dataset
Host_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#Nei heterozygosity for the full temporal dataset
Host_HetNei<-Hs(COMPDI)


##############################################################################/
#building the final table####
##############################################################################/

divtable<-cbind(Host_nb_samples,Host_nb_MLG,Host_GsurN,Host_MLG_richness,
                Host_simpson_div,Host_pielou_even,Host_Ar,Host_PrivAr,
                Host_HetNei)
write.table(divtable,file="output/divtable.txt",quote=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/