##############################################################################/
##############################################################################/
#Temporal diversity indices computation
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#MLG Diversity indices by semester####
##############################################################################/

#To compute the diversity indices by semester, we need a presence/absence of 
#MLG table for each semester
AgrOcc<-table(TempAgra$semester,TempAgra$MLG_ID)
nb_samples<-rowSums(AgrOcc)
nb_MLG<-specnumber(AgrOcc)
GsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
MLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
simpson_div<-diversity(AgrOcc,index="simpson")
pielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


##############################################################################/
#MLG Diversity indices by year####
##############################################################################/

#To compute the diversity indices by semester, we need a presence/absence of 
#MLG table for each semester
AgrOcc<-table(TempAgra$year,TempAgra$MLG_ID)
Ynb_samples<-rowSums(AgrOcc)
Ynb_MLG<-specnumber(AgrOcc)
YGsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
YMLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
Ysimpson_div<-diversity(AgrOcc,index="simpson")
Ypielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


##############################################################################/
#Genetic Diversity indices by semester####
##############################################################################/

#converting data to a genind format for the full dataset
compdiv<-TempAgra #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$semester,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$semester,ploidy=2,
                   ind.names=as.character(compdiv$indiv_ID),
                   NA.char=c("999999"))

#Allelic richness for the full temporal dataset
Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 5 (the minimum in the 
#clone-corrected dataset)
Ar<-apply(AllRich(COMPDI,5)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 5 (the minimum in the 
#clone-corrected dataset)
PrivAr<-apply(PrivAllRich(COMPDI,5)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
HetNeicc<-HeterNei(COMPDIcc)


##############################################################################/
#Genetic Diversity indices by year####
##############################################################################/

#converting data to a genind format for the full dataset
compdiv<-TempAgra #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$year,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-TempAgracc #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$year,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal dataset
YAr<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 5 (the minimum in the 
#clone-corrected dataset)
YAr<-apply(AllRich(COMPDI,19)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
YArcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
YPrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 5 (the minimum in the 
#clone-corrected dataset)
YPrivAr<-apply(PrivAllRich(COMPDI,19)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
YPrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
YHetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
YHetNeicc<-HeterNei(COMPDIcc)


##############################################################################/
#END
##############################################################################/