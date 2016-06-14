###############################################################################
###############################################################################
#Diversity indices computation for resistotype
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)
library(vegan)
library(combinat)
library(pegas)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#preparing the datasets
###############################################################################

RezAgra<-TempAgra
levels(RezAgra$KDR)[levels(RezAgra$KDR)!="GAGGAG"] <- "1"
levels(RezAgra$KDR)[levels(RezAgra$KDR)=="GAGGAG"] <- "0"

levels(RezAgra$sKDR)[levels(RezAgra$sKDR)!="CATCAT"] <- "1"
levels(RezAgra$sKDR)[levels(RezAgra$sKDR)=="CATCAT"] <- "0"

levels(RezAgra$MACE)[levels(RezAgra$MACE)!="TCATCA"] <- "1"
levels(RezAgra$MACE)[levels(RezAgra$MACE)=="TCATCA"] <- "0"

RezAgra<-cbind(RezAgra,
               "resisto"=paste(RezAgra$KDR,RezAgra$sKDR,RezAgra$MACE,sep=""))
RezAgra<-RezAgra[!RezAgra$resisto %in% c("1NA0","NANA0","NANANA"),]
RezAgra<-drop.levels(RezAgra)

t(table(RezAgra$resisto,RezAgra$semester))
plot(t(table(RezAgra$resisto,RezAgra$semester))[,1],type="b")
plot(t(table(RezAgra$resisto,RezAgra$semester))[,4],type="b")


#the clone-corrected dataset
RezAgracc<-TempAgracc
levels(RezAgracc$KDR)[levels(RezAgracc$KDR)!="GAGGAG"] <- "1"
levels(RezAgracc$KDR)[levels(RezAgracc$KDR)=="GAGGAG"] <- "0"

levels(RezAgracc$sKDR)[levels(RezAgracc$sKDR)!="CATCAT"] <- "1"
levels(RezAgracc$sKDR)[levels(RezAgracc$sKDR)=="CATCAT"] <- "0"

levels(RezAgracc$MACE)[levels(RezAgracc$MACE)!="TCATCA"] <- "1"
levels(RezAgracc$MACE)[levels(RezAgracc$MACE)=="TCATCA"] <- "0"

RezAgracc<-cbind(RezAgracc,
               "resisto"=paste(RezAgracc$KDR,RezAgracc$sKDR,RezAgracc$MACE,
                               sep=""))
RezAgracc<-RezAgracc[!RezAgracc$resisto %in% c("1NA0","NANA0","NANANA"),]
RezAgracc<-drop.levels(RezAgracc)


###############################################################################
#MLG Diversity indices by resistotypes
###############################################################################

AgrOcc<-table(RezAgra$resisto,RezAgra$MLG_ID)
Reznb_samples<-rowSums(AgrOcc)
Reznb_MLG<-specnumber(AgrOcc)
RezGsurN<-specnumber(AgrOcc)/rowSums(AgrOcc)
RezMLG_richness<-rarefy(AgrOcc,min(rowSums(AgrOcc)))
Rezsimpson_div<-diversity(AgrOcc,index="simpson")
Rezpielou_even<-diversity(AgrOcc)/log(specnumber(AgrOcc))


###############################################################################
#Genetic Diversity indices by Cluster k=3
###############################################################################

#converting data to a genind format for the full dataset
compdiv<-RezAgra[!is.na(RezAgra$resisto),] #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$resisto,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-RezAgracc[!is.na(RezAgracc$resisto),] #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$resisto,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal dataset
Rez_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 1 (the minimum in the 
#clone-corrected dataset)
Rez_Ar<-apply(AllRich(COMPDI,1)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
Rez_Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
Rez_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 1 (the minimum in the 
#clone-corrected dataset)
Rez_PrivAr<-apply(PrivAllRich(COMPDI,1)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
Rez_PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
Rez_HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
Rez_HetNeicc<-HeterNei(COMPDIcc)


###############################################################################
#END
###############################################################################