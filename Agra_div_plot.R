###############################################################################
###############################################################################
#Script for the plot of the variation of the diversity indices
###############################################################################
###############################################################################

#before using the following code, you have to run the Agra_analyze.R file


###############################################################################
#Distribution of the number of repetition of the different MLG
###############################################################################

#distribution for the complete dataset
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))]),
        "Clust_K3"])
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))]),
        "Clust_K5"])

#only MLG that are repeated more than once
summary(table(datAgra$MLG_ID)>1)
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74]),
        "Clust_K3"])
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74]),
        "Clust_K5"])

#distribution for the aerial trap samples only
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))]),
        "Clust_K3"])
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))]),
        "Clust_K5"])
#only MLG that are repeated more than once
summary(table(TempAgra$MLG_ID)>1)
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32]),
        "Clust_K3"])
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32],
        col=datAgracc[datAgracc$MLG_ID %in% 
        names(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32]),
        "Clust_K5"])


###############################################################################
#Distribution of the number of individuals in the different genetic clusters
###############################################################################

#plot the of the amount of the different genetic clusters with K=3
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),
        col=rainbow(5)[1:3],beside=TRUE)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),
        col=rainbow(5)[1:3],beside=FALSE)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K3)),
        col=rainbow(5)[1:3],beside=FALSE)

op<-par(mfrow=c(3,1),mar=c(0,0,1,0),oma=c(2,2,0,0))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[1,],col=rainbow(5)[1],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[2,],col=rainbow(5)[2],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[3,],col=rainbow(5)[3],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
par(op)

#plot the of the amount of the different genetic clusters with K=5
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5)),
        col=rainbow(5)[1:5],beside=TRUE)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5)),
        col=rainbow(5)[1:5],beside=FALSE)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K5)),
        col=rainbow(5)[1:5],beside=FALSE)

op<-par(mfrow=c(5,1),mar=c(0,0,1,0),oma=c(2,2,0,0))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[1,],col=rainbow(5)[1],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[2,],col=rainbow(5)[2],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[3,],col=rainbow(5)[3],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[4,],col=rainbow(5)[4],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[5,],col=rainbow(5)[5],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
par(op)


###############################################################################
#Plot diversity indices by semester
###############################################################################

def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,3,4,5),5,1,byrow=TRUE),
       widths=c(2,2,2,2,2),heights=c(1,1,1,1,1))
par(mar = c(2,2,1,1))
barplot(rowSums(AgrOcc),axes=TRUE,axisnames=FALSE,space=0,xlim=c(0.5,13.5))
plot(specnumber(AgrOcc)/rowSums(AgrOcc),type="b",main="G/N",xlim=c(1,14))
plot(rarefy(AgrOcc,min(rowSums(AgrOcc))),type="b",main="MLG richness",
     xlim=c(1,14))
#plot(diversity(AgrOcc,index="shannon"),type="b",main="Shannon",xlim=c(1,14))
plot(diversity(AgrOcc,index="simpson"),type="b",main="Simpson",xlim=c(1,14))
#plot(diversity(AgrOcc,index="invsimpson"),type="b",main="Invert Simpson")
plot(diversity(AgrOcc)/log(specnumber(AgrOcc)),type="b", 
     main="Pielou's evenness",xlim=c(1,14))
par(def.par)


op<-par(mfrow=c(4,2))

#plot of the number of samples for each period of time
barplot(nb_samples,axes=TRUE,axisnames=FALSE,space=0,xlim=c(0.5,13.5), 
        main="Number of samples")

#plot of the evolution of G/N
plot(GsurN,type="b",ylim=c(0.3,1),main="G/N",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(GsurN),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(GsurN,type="b")
box()

#plot of the evolution of MLG richness
plot(MLG_richness,type="b",ylim=c(6,14),main="MLG Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(MLG_richness),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(MLG_richness,type="b")
box()

#plot of the evolution of Simpson index
plot(simpson_div,type="b",ylim=c(0.75,1),main="Simpson Index",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(simpson_div),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(simpson_div,type="b")
box()

#plot of the evolution of Pielou's evenness
plot(pielou_even,type="b",ylim=c(0.75,1),main="Pielou's Evenness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(pielou_even),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(pielou_even,type="b")
box()

#plot of the evolution of allelic richness
plot(Ar,type="b",ylim=c(2.5,5.5),main="Allelic Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(Ar),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(Ar,type="b")
lines(Arcc,type="b",col="red")
box()

#plot of the evolution of private allelic richness
plot(PrivAr,type="b",ylim=c(0.05,0.6),main="Private Allelic Richness",
     axes=FALSE,xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(PrivAr),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(PrivAr,type="b")
lines(PrivArcc,type="b",col="red")
box()

#plot of the evolution of heterozygosity
plot(HetNei,type="b",ylim=c(0.6,0.85),main="Heterozygosity",axes=FALSE,
     xlab="",ylab="")
axis(side=2,las=1)
axis(side=1,labels=names(HetNei),at=1:14,las=2)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(HetNei,type="b")
lines(HetNeicc,type="b",col="red")
box()

par(op)


###############################################################################
#END
###############################################################################