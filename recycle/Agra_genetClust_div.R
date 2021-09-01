##############################################################################/
##############################################################################/
#Diversity indices computation for genetic clusters
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#MLG Diversity indices by genetic clusters####
##############################################################################/

#For K=3 genetic clusters
AgrOcc<-table(TempAgra$Clust_K3,TempAgra$MLG_ID)
K3nb_samples<-rowSums(AgrOcc)
K3nb_MLG<-specnumber(AgrOcc)
K3GsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
K3MLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
K3simpson_div<-diversity(AgrOcc,index="simpson")
K3pielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))

#For K=4 genetic clusters
AgrOcc<-table(TempAgra$Clust_K4,TempAgra$MLG_ID)
K5nb_samples<-rowSums(AgrOcc)
K5nb_MLG<-specnumber(AgrOcc)
K5GsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
K5MLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
K5simpson_div<-diversity(AgrOcc,index="simpson")
K5pielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


##############################################################################/
#Genetic Diversity indices by Cluster K=3####
##############################################################################/

#converting data to a genind format for the full data set
compdiv<-TempAgra[!is.na(TempAgra$Clust_K3),] #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K3,
                  ploidy=2,NA.char=c("999999"),
                  ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc[!is.na(TempAgracc$Clust_K3),] #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K3,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal data set
K3_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 34 (the minimum in the 
#clone-corrected data set)
K3_Ar<-apply(AllRich(COMPDI,34)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal data set
K3_Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal data set
K3_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 34 (the minimum in the 
#clone-corrected dataset)
K3_PrivAr<-apply(PrivAllRich(COMPDI,34)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal data set
K3_PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal data set
K3_HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal data set
K3_HetNeicc<-HeterNei(COMPDIcc)


##############################################################################/
#Genetic Diversity indices by Cluster K=4####
##############################################################################/

#converting data to a genind format for the full data set
compdiv<-TempAgra[!is.na(TempAgra$Clust_K4),] #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K5,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc[!is.na(TempAgracc$Clust_K4),] #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$Clust_K4,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal data set
K5_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 14 (the minimum in the 
#clone-corrected data set)
K5_Ar<-apply(AllRich(COMPDI,14)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal data set
K5_Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal data set
K5_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 14 (the minimum in the 
#clone-corrected data set)
K5_PrivAr<-apply(PrivAllRich(COMPDI,14)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal data set
K5_PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal data set
K5_HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal data set
K5_HetNeicc<-HeterNei(COMPDIcc)


##############################################################################/
#END
##############################################################################/