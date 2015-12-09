###############################################################################
###############################################################################
#STRUCTURE-like plot function
###############################################################################
###############################################################################


###############################################################################
#Defining a function to make structure-like plot
###############################################################################

structplot<-function(qmat,coolcol,effP,nameP,leg_y="",cexpop=1,cexy=2,
                     mef=c(1,1,1,1,1),colbord="NA",angl=0,distxax=0.005,
                     spacepop=0)
  #'qmat': the q-matrix like matrix
  #'coolcol': a vector of colors, for the different genetic clusters
  #'effP': a vector giving the number of individuals in each population
  #'nameP': a list giving the names of the different populations
  #'leg_y': a characters string used for the Y-axis legend
  #'cexpop': the cex factor for the population name on the x-axis
  #'cexy': the cex factor for the Y-legend
  #'mef': a vector of length 5 to pimp the graph. Each 1 value add a feature
  #the first is for the external rectangle, the second is to deliminate the 
  #different populations, the third is for adding an x-axis with tick, the 
  #fourth is for adding the name of the different populations and the fifth 
  #is for adding a Y legend
  #'colbord': the color of the line between individuals
  #'angl': the angle of the tag of the x-axis
  #'distxax': control the distance of the tag to the x-axis
  #'spacepop': space between populations, in number of bars

{
  vecspace<-c(rep(0,effP[1]))
  for (i in 1:(length(effP)-1)) {
    vecspace<-c(vecspace,spacepop,rep(0,effP[i+1]-1))
  }
  
  barplot(qmat,col=coolcol,beside=FALSE,border=colbord,
          space=vecspace,ylim=c(-0.03,1.03),axisnames=FALSE,axes=FALSE)
  
  if(mef[1]==1) {
    #drawing an external rectangle
    rect(0-1/dim(qmat)[2],
         0-1/500,
         dim(qmat)[2]+spacepop*(length(effP)-1)+1/dim(qmat)[2],
         1+1/500,
         lwd=3)
  }
  
  if(mef[2]==1) {
    decal<-c(0,cumsum(rep(spacepop,length(effP)-1)))
    #deliminated the different populations
    rect(c(0,cumsum(effP))[1:length(effP)]+decal,
         rep(0,length(effP)),
         cumsum(effP)+decal,
         rep(1,length(effP)),
         lwd=2)
  }
  
  if(mef[3]==1) {
    decal<-c(0,cumsum(rep(spacepop,length(effP)-1)))
    #add an x-axis
    axis(1,at=c(0,cumsum(effP))[1:length(effP)]+decal+
           (cumsum(effP)-c(0,cumsum(effP))[1:length(effP)])/2,
         labels=FALSE,pos=0,lwd.ticks=2,lwd=0)
  }
  
  if(mef[4]==1) {
    decal<-c(0,cumsum(rep(spacepop,length(effP)-1)))
    #add the name of the different populations
    text(c(0,cumsum(effP))[1:length(effP)]+decal+
           (cumsum(effP)-c(0,cumsum(effP))[1:length(effP)])/2,
         rep(par("usr")[3]-distxax,length(effP)),labels=nameP,srt=angl,
         xpd=NA,pos=1,cex=cexpop)
  }
  
  if(mef[5]==1) {
    #add some legend on the Y-axis
    mtext(leg_y,side=2,las=1,cex=cexy,adj=0.5,line=1)
  }
  
}


###############################################################################
#END
###############################################################################