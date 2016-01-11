###############################################################################
###############################################################################
#Script for the plot of the variation of the diversity indices
###############################################################################
###############################################################################

#before using the following code, you have to run the Agra_analyze.R file
#define a set of colors to be consistent across the plots
coloor <- c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")


###############################################################################
#Distribution of the different genetic cluster by sampled host
###############################################################################

#for K=3 genetic clusters
table(datAgra$host,datAgra$Clust_K3)
#distribution of genetic cluster for each sampled host
op<-par(mfrow=c(5,1),mar=c(0,0,3,0),oma=c(4,3,0,0))
barplot((table(datAgra$host,datAgra$Clust_K3))[4,],las=1,
        col=coloor[c(1,3,2)],main="Peach tree samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K3))[2,],las=1,
        col=coloor[c(1,3,2)],main="Oilseed rape samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K3))[5,],las=1,
        col=coloor[c(1,3,2)],main="Tobacco samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K3))[3,],las=1,
        col=coloor[c(1,3,5,4,2)],main="Other crops samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K3))[1,],las=1,
        col=coloor[c(1,3,2)],main="Aerial samples",axisnames=FALSE)
Xsemest<-barplot(table(datAgra$host,datAgra$Clust_K3),plot=FALSE,beside=FALSE)
text(x=Xsemest,y=rep(par("usr")[3],5)-(par("usr")[4]-par("usr")[3])/20,
     labels=c("Cluster\nPrimary","Cluster\nSecondary","Cluster\nWild"),
     srt=0,xpd=NA,pos=1,cex=1.2)
par(op)

#export 3.6 x 18

#distribution of sampled host by genetic cluster
t(table(datAgra$host,datAgra$Clust_K3))
op<-par(mfrow=c(3,1),mar=c(0,0,3,0),oma=c(4,3,0,0))
barplot(t(table(datAgra$host,datAgra$Clust_K3))[1,c(4,2,5,3,1)],las=1,
        col=coloor[1],main="Cluster Primary",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K3))[2,c(4,2,5,3,1)],las=1,
        col=coloor[3],main="Cluster Secondary",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K3))[3,c(4,2,5,3,1)],las=1,
        col=coloor[2],main="Cluster Wild",axisnames=FALSE)
Xsemest<-barplot(t(table(datAgra$host,datAgra$Clust_K5)),plot=FALSE,
                 beside=FALSE)
text(x=Xsemest,y=rep(par("usr")[3],5)-(par("usr")[4]-par("usr")[3])/20,
     labels=c("Peach tree\nsamples","Oilseed rape\nsamples",
              "Tobacco\nsamples","Other crops\nsamples","Aerial\nsamples"),
     srt=0,xpd=NA,pos=1,cex=1.2)
par(op)

#export pdf 6 x 6

#for K=5 genetic clusters
table(datAgra$host,datAgra$Clust_K5)
#distribution of genetic cluster for each sampled host
op<-par(mfrow=c(5,1),mar=c(0,0,3,0),oma=c(4,3,0,0))
barplot((table(datAgra$host,datAgra$Clust_K5))[4,c(3,1,2,4,5)],las=1,
        col=coloor[c(1,3,5,4,2)],main="Peach tree samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K5))[2,c(3,1,2,4,5)],las=1,
        col=coloor[c(1,3,5,4,2)],main="Oilseed rape samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K5))[5,c(3,1,2,4,5)],las=1,
        col=coloor[c(1,3,5,4,2)],main="Tobacco samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K5))[3,c(3,1,2,4,5)],las=1,
        col=coloor[c(1,3,5,4,2)],main="Other crops samples",axisnames=FALSE)
barplot((table(datAgra$host,datAgra$Clust_K5))[1,c(3,1,2,4,5)],las=1,
        col=coloor[c(1,3,5,4,2)],main="Aerial samples",axisnames=FALSE)
Xsemest<-barplot(table(datAgra$host,datAgra$Clust_K5),plot=FALSE,beside=FALSE)
text(x=Xsemest,y=rep(par("usr")[3],5)-(par("usr")[4]-par("usr")[3])/20,
     labels=c("Cluster\nPeach","Cluster\nOilseed1","Cluster\nOilseed2",
              "Cluster\nTobacco","Cluster\nEnigma"),
     srt=0,xpd=NA,pos=1,cex=1.2)
par(op)

#distribution of sampled host by genetic cluster
t(table(datAgra$host,datAgra$Clust_K5))
op<-par(mfrow=c(5,1),mar=c(0,0,3,0),oma=c(4,3,0,0))
barplot(t(table(datAgra$host,datAgra$Clust_K5))[3,c(4,2,5,3,1)],las=1,
        col=coloor[1],main="Cluster Peach",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K5))[1,c(4,2,5,3,1)],las=1,
        col=coloor[3],main="Cluster Oilseed1",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K5))[2,c(4,2,5,3,1)],las=1,
        col=coloor[5],main="Cluster Oilseed2",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K5))[4,c(4,2,5,3,1)],las=1,
        col=coloor[4],main="Cluster Tobacco",axisnames=FALSE)
barplot(t(table(datAgra$host,datAgra$Clust_K5))[5,c(4,2,5,3,1)],las=1,
        col=coloor[2],main="Cluster Enigma",axisnames=FALSE)
Xsemest<-barplot(t(table(datAgra$host,datAgra$Clust_K5)),plot=FALSE,
                 beside=FALSE)
text(x=Xsemest,y=rep(par("usr")[3],5)-(par("usr")[4]-par("usr")[3])/20,
     labels=c("Peach tree\nsamples","Oilseed rape\nsamples",
              "Tobacco\nsamples","Other crops\nsamples","Aerial\nsamples"),
     srt=0,xpd=NA,pos=1,cex=1.2)
par(op)

#export to pdf 6 x 18 inches


###############################################################################
#Distribution of the number of repetition of the different MLG: complete data
###############################################################################

barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))])
#for assigning a color according to the genetic cluster, we 
#build a table of MLG x Cluster belonging
clustbelong<-datAgracc[,c("MLG_ID","Clust_K3","Clust_K5")]
levels(clustbelong$Clust_K3)<-c("firebrick","chartreuse4","royalblue4")
levels(clustbelong$Clust_K5)<-c("chartreuse4","darkorange","firebrick","khaki2",
                                "royalblue4")
row.names(clustbelong)<-clustbelong$MLG_ID
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))],
        col=as.character(clustbelong[names(table(datAgra$MLG_ID)
                              [order(-table(datAgra$MLG_ID))]),"Clust_K3"]))
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))],
        col=as.character(clustbelong[names(table(datAgra$MLG_ID)
                              [order(-table(datAgra$MLG_ID))]),"Clust_K5"]))

#only with MLG that are repeated more than once
summary(table(datAgra$MLG_ID)>1)
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74],
        col=as.character(clustbelong[names(table(datAgra$MLG_ID)
                              [order(-table(datAgra$MLG_ID))]),
                        "Clust_K3"][1:74]),cex.names=0.8,las=2)
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(datAgra$MLG_ID)[order(-table(datAgra$MLG_ID))][1:74],
        col=as.character(clustbelong[names(table(datAgra$MLG_ID)
                              [order(-table(datAgra$MLG_ID))]),
                        "Clust_K5"][1:74]),cex.names=0.8,las=2)


###############################################################################
#Distribution of the number of repetition of the different MLG: aerial samples
###############################################################################

barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))])
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))],
        col=as.character(clustbelong[names(table(TempAgra$MLG_ID)
                              [order(-table(TempAgra$MLG_ID))]),"Clust_K3"]))
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))],
        col=as.character(clustbelong[names(table(TempAgra$MLG_ID)
                              [order(-table(TempAgra$MLG_ID))]),"Clust_K5"]))
#only MLG that are repeated more than once
summary(table(TempAgra$MLG_ID)>1)
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32],
        cex.names=0.8,las=2)
#same figure but MLG are colored according to the cluster to which they belong
#for K=3
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32],
        col=as.character(clustbelong[names(table(TempAgra$MLG_ID)
                              [order(-table(TempAgra$MLG_ID))])[1:32],
                        "Clust_K3"]),cex.names=0.8,las=2)
#same figure but MLG are colored according to the cluster to which they belong
#for K=5
barplot(table(TempAgra$MLG_ID)[order(-table(TempAgra$MLG_ID))][1:32],
        col=as.character(clustbelong[names(table(TempAgra$MLG_ID)
                              [order(-table(TempAgra$MLG_ID))])[1:32],
                        "Clust_K5"]),cex.names=0.8,las=2)


###############################################################################
#Distribution of the number of individuals in the different genetic clusters
###############################################################################

#plot the of the amount of the different genetic clusters with K=3
#all aerial trapped individuals
Xsemest<- barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),plot=FALSE,
                  col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),
        col=coloor[c(1,3,2)],beside=TRUE,las=1)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),
        col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE,las=1)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/30,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)

#clone-corrected by year aerial trapped individuals
#the x-coordinates of the 'by semester' barplot
Xsemest<- barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),plot=FALSE,
                  col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K3)),
        col=coloor[c(1,3,2)],beside=TRUE,las=1)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K3)),
        col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE,las=1)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/30,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)

op<-par(mfrow=c(3,1),mar=c(0,0,1,0),oma=c(4,3,0,0))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[1,],col=coloor[1],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[2,],col=coloor[3],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[3,],col=coloor[2],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K3))))
Xsemest<-barplot(t(table(TempAgra$semester,TempAgra$Clust_K3))[3,],space=0,
                 plot=FALSE,beside=TRUE)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/10,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)
par(op)

#plot the of the amount of the different genetic clusters with K=5
#all aerial trapped individuals
Xsemest<- barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),plot=FALSE,
                  col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[c(3,1,2,4,5),],
        col=coloor[c(1,3,5,4,2)],beside=TRUE,las=1)
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[c(3,1,2,4,5),],
        col=coloor[c(1,3,5,4,2)],beside=FALSE,axisnames=FALSE,las=1)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/30,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)

#clone-corrected by year aerial trapped individuals
Xsemest<- barplot(t(table(TempAgra$semester,TempAgra$Clust_K3)),plot=FALSE,
                  col=coloor[c(1,3,2)],beside=FALSE,axisnames=FALSE)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K5))[c(3,1,2,4,5),],
        col=coloor[c(1,3,5,4,2)],beside=TRUE,las=1)
barplot(t(table(TempAgracc$semester,TempAgracc$Clust_K5))[c(3,1,2,4,5),],
        col=coloor[c(1,3,5,4,2)],beside=FALSE,axisnames=FALSE,las=1)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/30,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)

op<-par(mfrow=c(5,1),mar=c(0,0,1,0),oma=c(4,3,0,0))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[3,],col=coloor[1],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[1,],col=coloor[3],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[2,],col=coloor[5],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[4,],col=coloor[4],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[5,],col=coloor[2],
        beside=TRUE,axisnames=FALSE,ann=FALSE,axes=TRUE,space=0,las=1,
        ylim=c(0,max(table(TempAgra$semester,TempAgra$Clust_K5))))
Xsemest<-barplot(t(table(TempAgra$semester,TempAgra$Clust_K5))[3,],space=0,
                 plot=FALSE,beside=TRUE)
text(x=Xsemest-0.4,y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/10,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1)
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
axis(side=1,labels=FALSE,at=1:14,las=2)
text(x=(1:14),y=rep(par("usr")[3],14)-(par("usr")[4]-par("usr")[3])/15,
     labels=names(HetNei),srt=45,xpd=NA,pos=1,cex=1,adj=0)
barplot(c(rep(20,7)),axes=FALSE,axisnames=FALSE,space=c(0.5,1,1,1,1,1,1),
        add=TRUE,xpd=FALSE,border=NA,offset=-1)
lines(HetNei,type="b")
lines(HetNeicc,type="b",col="red")
box()

par(op)


###############################################################################
#END
###############################################################################