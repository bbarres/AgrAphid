##############################################################################/
##############################################################################/
#Script for the study of the distribution of the different genetic clusters
##############################################################################/
##############################################################################/

#before using the following code, you have to run the Agra_load.R and the 
#Agra_temp_div.R code
#define a set of colors to be consistent across the plots
coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")


##############################################################################/
#Distribution of the different genetic cluster by years combined monthly####
##############################################################################/

op<-par(mfrow=c(2,1))
barplot(table(datAgra[datAgra$data_batch=="AgrAphid" & 
                        !is.na(datAgra$data_batch),"Clust_K3"],
              months(datAgra[datAgra$data_batch=="AgrAphid" & 
                               !is.na(datAgra$data_batch),"sampling_date"])
              )[,c(5,4,8,1,9,7,6,2,12,11,10,3)],
        col=coloor[c(1,3,2)],beside=TRUE)

barplot(table(datAgra[datAgra$data_batch=="AgrAphid" & 
                        !is.na(datAgra$data_batch),"Clust_K5"],
              months(datAgra[datAgra$data_batch=="AgrAphid" & 
                               !is.na(datAgra$data_batch),"sampling_date"])
              )[c(3,1,2,4,5),c(5,4,8,1,9,7,6,2,12,11,10,3)],
        col=coloor[c(1,3,5,4,2)],beside=TRUE)
par(op)


##############################################################################/
#END
##############################################################################/