##############################################################################/
##############################################################################/
#Delta-K method plotting functions
##############################################################################/
##############################################################################/


##############################################################################/
#Computing Delta-k####
##############################################################################/

chooseK<-function(str_out,nb_K,nb_rep) {
  #'str_out': the summary of simulations file exported from STRUCTURE, 
  #before importing the file to R, replace the white space in the column 
  #header by underscores, you might also have to replace "?1" by "alpha" 
  #and remove double white space
  #'nb_K': number of different K values considered
  #'nb_rep':number of repetition for each K
  datatable<-data.frame("K"=c(rep(1:nb_K,each=nb_rep)),"Ln(Pd)"=str_out[,4])
  Lprim<-c(rep("NA",nb_rep))
  for (i in ((nb_rep+1):(nb_K*nb_rep))) {
    Lprim<-c(Lprim,str_out[i,4]-str_out[i-nb_rep,4])
  }
  datatable<-data.frame(datatable,as.numeric(Lprim))
  Lsecond<-c(rep("NA",nb_rep))
  for (i in (((2*nb_rep)+1):(nb_K*nb_rep))) {
    Lsecond<-c(Lsecond,abs(datatable[i,3]-datatable[i-nb_rep,3]))
  }
  Lsecond<-c(Lsecond,rep("NA",nb_rep))
  datatable<-data.frame(datatable,as.numeric(Lsecond))
  reztable<-data.frame("K"=c(1:nb_K))
  meanL<-c()
  sdL<-c()
  for (i in (1:nb_K)) {
    meanL<-c(meanL,mean(datatable[datatable$K==i,2]))
    sdL<-c(sdL,sd(datatable[datatable$K==i,2]))
  }
  reztable<-data.frame(reztable,meanL,sdL)
  meanLprime<-c()
  sdLprime<-c()
  for (i in (1:nb_K)) {
    meanLprime<-c(meanLprime,mean(as.numeric(datatable[datatable$K==i,3])))
    sdLprime<-c(sdLprime,sd(datatable[datatable$K==i,3]))
  }
  reztable<-data.frame(reztable,meanLprime,sdLprime)
  meanLsecond<-c()
  sdLsecond<-c()
  for (i in (1:nb_K)) {
    meanLsecond<-c(meanLsecond,mean(as.numeric(datatable[datatable$K==i,4])))
    sdLsecond<-c(sdLsecond,sd(datatable[datatable$K==i,4]))
  }
  reztable<-data.frame(reztable,meanLsecond,sdLsecond)
  deltaK<-c()
  for (i in (1:nb_K)) {
    deltaK<-c(deltaK,reztable[reztable$K==i,6]/reztable[reztable$K==i,3])
  }
  reztable<-data.frame(reztable,deltaK)
  return(reztable)
}


##############################################################################/
#function to plot variation of Delta K and Ln(P(X|K)) as a function of K####
##############################################################################/

plotdeltaK<-function(datadeltaK,nb_K,titre){
  #'datadeltak': the output file of 'chooseK' function
  #'nb_K': the number of different K considered
  #'titre': the title of the plot you want to be displayed
  op<-par(pty="s")
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,cex=2.5,lwd=4,lty=1,
       col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,bty="n",xaxt="n",
       yaxt="n",ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Delta K",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}

#the same function using log(deltaK), just in order to see smaller variation 
#of deltaK

plotlogdeltaK<-function(datadeltaK,nb_K,titre){
  #'datadeltak': the output file of 'chooseK' function
  #'nb_K': the number of different K considered
  #'titre': the title of the plot you want to be displayed
  op<-par(pty="s")
  plot(log(datadeltaK[1:(nb_K-2),8]+1),type="b",pch=24,cex=2.5,lwd=4,
       lty=1,col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(log(datadeltaK[1:(nb_K-2),8]+1),type="b",pch=24,bty="n",xaxt="n",
       yaxt="n",ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Ln(Delta K+1)",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}


##############################################################################/
#END
##############################################################################/