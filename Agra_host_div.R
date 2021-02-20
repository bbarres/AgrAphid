##############################################################################/
##############################################################################/
#By host diversity indices computation
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#MLG Diversity indices by semester####
##############################################################################/

#To compute the diversity indices by semester, we need a presence/absence of 
#MLG table for each semester
hostocc<-table(datAgra$host,datAgra$MLG_ID)
Host_nb_samples<-rowSums(hostocc)
Host_nb_MLG<-specnumber(hostocc)
Host_GsurN<-specnumber(hostocc)/rowSums(hostocc)
Host_MLG_richness<-rarefy(hostocc,min(rowSums(hostocc)))
Host_simpson_div<-diversity(hostocc,index="simpson")
Host_pielou_even<-diversity(hostocc)/log(specnumber(hostocc))


##############################################################################/
#Genetic Diversity indices by host on the complete dataset####
##############################################################################/

#converting data to a genind format for the full dataset
compdiv<-datAgra #name of the input file
COMPDI<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$host,ploidy=2,
                  NA.char=c("999999"),ind.names=as.character(compdiv$indiv_ID))
#converting data to a genind format
compdiv<-datAgra[datAgra$one_MLG_host==1,] #name of the input file
COMPDIcc<-df2genind(compdiv[,11:24],ncode=6,pop=compdiv$host,ploidy=2,
                    ind.names=as.character(compdiv$indiv_ID),
                    NA.char=c("999999"))

#Allelic richness for the full temporal dataset
Host_Ar<-apply(AllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 16 (the minimum in the 
#clone-corrected dataset)
Host_Ar<-apply(AllRich(COMPDI,16)[[2]],1,mean)
#Allelic richness for the clone-corrected temporal dataset
Host_Arcc<-apply(AllRich(COMPDIcc)[[2]],1,mean)

#Private Allelic richness for the full temporal dataset
Host_PrivAr<-apply(PrivAllRich(COMPDI)[[2]],1,mean)
#we set the minimum number of samples to 16 (the minimum in the 
#clone-corrected dataset)
Host_PrivAr<-apply(PrivAllRich(COMPDI,16)[[2]],1,mean)
#Private Allelic richness for the clone-corrected temporal dataset
Host_PrivArcc<-apply(PrivAllRich(COMPDIcc)[[2]],1,mean)

#Nei heterozygosity for the full temporal dataset
Host_HetNei<-HeterNei(COMPDI)
#Nei heterozygosity for the clone-corrected temporal dataset
Host_HetNeicc<-HeterNei(COMPDIcc)


##############################################################################/
#END
##############################################################################/