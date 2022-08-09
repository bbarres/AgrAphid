#this is a script which was written to produce a figure for a poster
#communication in ESEB 2022. The comment of the code are therefore 
#sloppy and some details are lackin. The code will be improved
#if I found time to do so...

##############################################################################/
#defining additional function for the mapping####
##############################################################################/

#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}


##############################################################################/
#Producing the map####
##############################################################################/


colovec<-c(brewer.pal(9,"Set1")[c(8,6,7,2,1)])

dataccc<-datAgra[!is.na(datAgra$longitude),]
oldSA.wgs<-SpatialPointsDataFrame(coords=dataccc[,c("longitude","latitude")],
                                  data=dataccc,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
)
oldSA<-spTransform(oldSA.wgs,CRS("+init=epsg:2154"))
#reorder the factors
oldSA$host<-factor(oldSA$host,levels=c("peach","oilseed_rape","tobacco",
                                       "other_crops","Aerial_trap"))

op<-par(mar=c(0,0,0,0))
plot(DEP_SHP,main="",border="grey70")
plot(REG_SHP,lwd=2,add=TRUE)
points(
  x=as.numeric(oldSA$longitude),
  y=as.numeric(oldSA$latitude),
  bg=colovec[as.numeric(oldSA$host)],  #colors of the points
  pch=21,                  #plotting character
  cex=1.3                           #size of the points
)
points(
  x=as.numeric(oldSA[oldSA$host=="Aerial_trap",]$longitude),
  y=as.numeric(oldSA[oldSA$host=="Aerial_trap",]$latitude),
  bg=colovec[5], col="white",  #colors of the points
  pch=24,                  #plotting character
  cex=1.4                           #size of the points
)
legend(100000,7150000,title="",
       legend=c("Peach tree","Oilseed rape","Tobacco",
                "Other crops","Aerial trap"),cex=0.8,
       pt.cex=c(1.4,1.4,1.4,1.4,1.5),
       y.intersp=1.2,x.intersp=0.9,
       pch=c(21,21,21,21,24),title.adj=0.3,
       pt.bg=colovec,col=c(rep("black",4),"white"),
       bty="n")
scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
par(op)


#export to .pdf 6 x 6 inches