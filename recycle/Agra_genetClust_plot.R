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
#Distribution of the different genetic cluster by month in the aerial trap####
##############################################################################/

#all individuals including repeated clones
op<-par(mfrow=c(2,1))

barplot(table(TempAgra[,"Clust_K3"],
              month(TempAgra[,"sampling_date"])),
        col=coloor[c(1,3,2)],beside=TRUE,las=1,ylim=c(0,85),
        names.arg=c("January","February","March","April","May","June",
                    "July","August","September","October","November",
                    "December"))
abline(h=20,lty=2,col=grey(0.7),lwd=2)
abline(h=40,lty=2,col=grey(0.7),lwd=2)
abline(h=60,lty=2,col=grey(0.7),lwd=2)
abline(h=80,lty=2,col=grey(0.7),lwd=2)
box(bty="l",lwd=2)

barplot(table(TempAgra[,"Clust_K4"],
              month(TempAgra[,"sampling_date"]))[c(2,1,3,4),],
        col=coloor[c(1,3,4,2)],beside=TRUE,las=1,ylim=c(0,85),
        names.arg=c("January","February","March","April","May","June",
                    "July","August","September","October","November",
                    "December"),
        xlab=c("Months"))
abline(h=20,lty=2,col=grey(0.7),lwd=2)
abline(h=40,lty=2,col=grey(0.7),lwd=2)
abline(h=60,lty=2,col=grey(0.7),lwd=2)
abline(h=80,lty=2,col=grey(0.7),lwd=2)
box(bty="l",lwd=2)

par(op)

#export to .pdf 14 x 12 inches


##############################################################################/
#END
##############################################################################/