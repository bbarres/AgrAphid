###############################################################################
###############################################################################
#Prepare data for Maurice
###############################################################################
###############################################################################

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/AgrAphid_data")


###############################################################################
#data by semester
###############################################################################

datasemester<-cbind("pop_ID"=names(nb_samples),
                    nb_samples,nb_MLG,GsurN,MLG_richness,simpson_div,
                    pielou_even,Ar,Arcc,PrivAr,PrivArcc,HetNei,HetNeicc,
                    table(TempAgra$semester,TempAgra$Clust_K3),
                    table(TempAgra$semester,TempAgra$Clust_K5),
                    table(TempAgra$semester,TempAgra$MLG_ID)
                    [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                    "nb_unik"=rowSums(table(TempAgra$semester,
                                            TempAgra$MLG_ID)
                    [,colSums(table(TempAgra$semester,TempAgra$MLG_ID))<2]))

#export the data table
write.table(datasemester,file="datasemester.txt",sep="\t",
            row.names=FALSE,quote=FALSE)

###############################################################################
#data by year
###############################################################################

datayear<-cbind("pop_ID"=names(Ynb_samples),
                Ynb_samples,Ynb_MLG,YGsurN,YMLG_richness,Ysimpson_div,
                Ypielou_even,YAr,YArcc,YPrivAr,YPrivArcc,YHetNei,YHetNeicc,
                table(TempAgra$year,TempAgra$Clust_K3),
                table(TempAgra$year,TempAgra$Clust_K5),
                table(TempAgra$year,TempAgra$MLG_ID)
                [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                "nb_unik"=rowSums(table(TempAgra$year,
                                        TempAgra$MLG_ID)
                [,colSums(table(TempAgra$year,TempAgra$MLG_ID))<2]))

#export the data table
write.table(datayear,file="datayear.txt",sep="\t",
            row.names=FALSE,quote=FALSE)


###############################################################################
#data by Cluster K3
###############################################################################

dataClustK3<-cbind("pop_ID"=names(K3nb_samples),
                K3nb_samples,K3nb_MLG,K3GsurN,K3MLG_richness,K3simpson_div,
                K3pielou_even,K3_Ar,K3_Arcc,K3_PrivAr,K3_PrivArcc,K3_HetNei,
                K3_HetNeicc,
                table(TempAgra$Clust_K3,TempAgra$year),
                table(TempAgra$Clust_K3,TempAgra$MLG_ID)
                [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                "nb_unik"=rowSums(table(TempAgra$Clust_K3,
                                        TempAgra$MLG_ID)
                 [,colSums(table(TempAgra$Clust_K3,TempAgra$MLG_ID))<2]))

#export the data table
write.table(dataClustK3,file="dataClustK3.txt",sep="\t",
            row.names=FALSE,quote=FALSE)


###############################################################################
#data by Cluster K5
###############################################################################

dataClustK5<-cbind("pop_ID"=names(K5nb_samples),
                   K5nb_samples,K5nb_MLG,K5GsurN,K5MLG_richness,K5simpson_div,
                   K5pielou_even,K5_Ar,K5_Arcc,K5_PrivAr,K5_PrivArcc,K5_HetNei,
                   K5_HetNeicc,
                   table(TempAgra$Clust_K5,TempAgra$year),
                   table(TempAgra$Clust_K5,TempAgra$MLG_ID)
                   [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                   "nb_unik"=rowSums(table(TempAgra$Clust_K5,
                                           TempAgra$MLG_ID)
                    [,colSums(table(TempAgra$Clust_K5,TempAgra$MLG_ID))<2]))

#export the data table
write.table(dataClustK5,file="dataClustK5.txt",sep="\t",
            row.names=FALSE,quote=FALSE)


###############################################################################
#END
###############################################################################