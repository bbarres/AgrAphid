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

#converting data to a genind format, first we use only the microsatellite data
temp_ClustK3<-datAgra[datAgra$host=="Aerial_trap",]
temp_ClustK3$Clust_K3<-as.character(temp_ClustK3$Clust_K3)
temp_ClustK3$Clust_K3[is.na(temp_ClustK3$Clust_K3)]<-"undef"
temp_ClustK3<-drop.levels(temp_ClustK3)
temp_ClustK3<-df2genind(temp_ClustK3[,c("MP_27","MP_39","MP_44","MP_5",
                                        "MP_7","MP_23","MP_45","MP_28",
                                        "MP_9","MP_13","MP_2","MP_38",
                                        "MP_4","MP_46")],
                        ncode=6,
                        ind.names=temp_ClustK3$MLG_ID, 
                        pop=temp_ClustK3$Clust_K3,
                        ploidy=2,NA.char="999999")
summary(temp_ClustK3)

bruvo.msn(temp_ClustK3, replen = rep(1, 14),vertex.label="inds", 
          vertex.label.cex=0.7, vertex.label.dist=0.4,cex=0.5)







##############################################################################/
#END
##############################################################################/