###############################################################################
###############################################################################
#Script for the study of the repartition of the different genetic clusters
###############################################################################
###############################################################################

#before using the following code, you have to run the Agra_load.R and the 
#Agra_temp_div.R code
#define a set of colors to be consistent across the plots
coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")


###############################################################################
#Distribution of the different genetic cluster for every years combined monthly
###############################################################################

barplot(table(datAgra[datAgra$data_batch=="AgrAphid" & 
                        !is.na(datAgra$data_batch),"Clust_K3"],
              months(datAgra[datAgra$data_batch=="AgrAphid" & 
                               !is.na(datAgra$data_batch),"sampling_date"])
              )[,c(5,4,9,2,8,7,6,1,12,11,10,3)],col=coloor)





###############################################################################
#END
###############################################################################