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

datasemester<-cbind(nb_samples,nb_MLG,GsurN,MLG_richness,simpson_div,
                    pielou_even,Ar,Arcc,PrivAr,PrivArcc,HetNei,HetNeicc,
                    table(TempAgra$semester,TempAgra$Clust_K3),
                    table(TempAgra$semester,TempAgra$Clust_K5),
                    table(TempAgra$semester,TempAgra$MLG_ID)
                    [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                    "nb_unik"=rowSums(table(TempAgra$semester,
                                            TempAgra$MLG_ID)
                    [,colSums(table(TempAgra$semester,TempAgra$MLG_ID))<2]))


###############################################################################
#data by year
###############################################################################

datayear<-cbind(Ynb_samples,Ynb_MLG,YGsurN,YMLG_richness,Ysimpson_div,
                Ypielou_even,YAr,YArcc,YPrivAr,YPrivArcc,YHetNei,YHetNeicc,
                table(TempAgra$year,TempAgra$Clust_K3),
                table(TempAgra$year,TempAgra$Clust_K5),
                table(TempAgra$year,TempAgra$MLG_ID)
                [,order(-table(TempAgra$MLG_ID))[1:majMLG]],
                "nb_unik"=rowSums(table(TempAgra$year,
                                        TempAgra$MLG_ID)
                [,colSums(table(TempAgra$year,TempAgra$MLG_ID))<2]))


###############################################################################
#END
###############################################################################
