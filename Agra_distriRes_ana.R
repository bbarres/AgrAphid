##############################################################################/
##############################################################################/
#Plotting the results of STRUCTURE post CLUMPP
##############################################################################/
##############################################################################/

source("Agra_load.R")
#defining a vector of colors
coloor<-c(brewer.pal(9,"Purples")[7],
          brewer.pal(9,"Oranges")[6],
          brewer.pal(9,"Greens")[5])

##############################################################################/
#preparing the data set for the comparison of resistance genotype####
##############################################################################/

#preparing the data set
temp<-as.data.table(datAgraccHost)
#reformatting KDR genotypes
temp$KDRg<-temp$KDR
levels(temp$KDRg)<-c("K-RR","K-RS","K-SS")
temp$KDRg<-as.character(temp$KDRg)
temp$KDRg[is.na(temp$KDRg)]<-"K-miss"
temp<-spread(temp,KDRg,KDRg)
#reformatting sKDR genotypes
temp$sKDRg<-temp$sKDR
levels(temp$sKDRg)<-c("sK-RR","sK-RR","sK-RR","sK-RS","sK-RS",
                      "sK-SS","sK-RS","sK-RS","sK-RS")
temp$sKDRg<-as.character(temp$sKDRg)
temp$sKDRg[is.na(temp$sKDRg)]<-"sK-miss"
temp<-spread(temp,sKDRg,sKDRg)
#reformatting MACE genotypes
temp$MACEg<-temp$MACE
levels(temp$MACEg)<-c("M-SS","M-RS")
temp$MACEg<-as.character(temp$MACEg)
temp$MACEg[is.na(temp$MACEg)]<-"M-miss"
temp<-spread(temp,MACEg,MACEg)
#reformatting R81T genotypes
temp$R81Tg<-temp$R81T
levels(temp$R81Tg)<-c("Neo-RR","Neo-RS","Neo-SS")
temp$R81Tg<-as.character(temp$R81Tg)
temp$R81Tg[is.na(temp$R81Tg)]<-"Neo-miss"
temp<-spread(temp,R81Tg,R81Tg)

#turning the resistance data into a 0-1 matrix 
temp<-cbind(temp,ifelse(is.na(temp[,c(64:78)]),0,1))
temp<-temp[,-c(64:78)]
#reordering the levels of the host
temp$host<-factor(temp$host,levels(temp$host)[c(4,2,5,3,1)])
#limiting the data set to necessary columns
temp<-temp[,c(1,7,43,64:78)]


##############################################################################/
#testing distribution of kdr genotypes by hosts####
##############################################################################/

temp2<-temp[temp$host!="Aerial_trap" & temp$host!="other_crops" 
            & temp$`K-miss`!=1,
            c(2,5:7)]
temp2<-drop.levels(temp2)

dd<-temp2 %>% 
  group_by(host) %>% 
  summarise(KRR=sum(`K-RR`),KRS=sum(`K-RS`),KSS=sum(`K-SS`))
dd

zz<-as.matrix(dd[c(2,1,3),2:4])
dimnames(zz)[[1]]<-c("Peach","Oilseed rape","tobacco")
zz
barplot(t(zz),beside=TRUE,col=coloor,main="KDR by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
fisher.test(zz)


##############################################################################/
#testing distribution of skdr genotypes by hosts####
##############################################################################/

temp2<-temp[temp$host!="Aerial_trap" & temp$host!="other_crops" 
            & temp$`sK-miss`!=1,
            c(2,9:11)]
temp2<-drop.levels(temp2)
dd<-temp2 %>% 
  group_by(host) %>% 
  summarise(sKRR=sum(`sK-RR`),sKRS=sum(`sK-RS`),sKSS=sum(`sK-SS`))
dd

xx<-as.matrix(dd[c(2,1,3),2:4])
dimnames(xx)[[1]]<-c("Peach","Oilseed rape","tobacco")
xx
barplot(t(xx),beside=TRUE,col=coloor,main="sKDR by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
fisher.test(xx)


##############################################################################/
#testing distribution of MACE genotypes by hosts####
##############################################################################/

temp2<-temp[temp$host!="Aerial_trap" & temp$host!="other_crops" 
            & temp$`M-miss`!=1,
            c(2,13:14)]
temp2<-drop.levels(temp2)
dd<-temp2 %>% 
  group_by(host) %>% 
  summarise(MRS=sum(`M-RS`),MSS=sum(`M-SS`))
dd

yy<-as.matrix(dd[c(2,1,3),2:3])
dimnames(yy)[[1]]<-c("Peach","Oilseed rape","tobacco")
yy
barplot(t(yy),beside=TRUE,col=coloor[c(2:3)],main="MACE by host")
legend(legend=c("RS","SS"),fill=coloor[c(2:3)],x="topright")
fisher.test(yy)


##############################################################################/
#testing distribution of R81T genotypes by hosts####
##############################################################################/

temp2<-temp[temp$host!="Aerial_trap" & temp$host!="other_crops" 
            & temp$`Neo-miss`!=1,
            c(2,16:18)]
temp2<-drop.levels(temp2)
dd<-temp2 %>% 
  group_by(host) %>% 
  summarise(NeoRR=sum(`Neo-RR`),NeoRS=sum(`Neo-RS`),NeoSS=sum(`Neo-SS`))
dd

ww<-as.matrix(dd[c(2,1,3),2:4])
dimnames(ww)[[1]]<-c("Peach","Oilseed rape","tobacco")
ww
barplot(t(ww),beside=TRUE,col=coloor,main="R81T by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
fisher.test(ww)


##############################################################################/
#combined plots of R genotype distribution loci by host####
##############################################################################/

op<-par(mfrow=c(2,2))
barplot(t(zz),beside=TRUE,col=coloor,main="KDR by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
barplot(t(xx),beside=TRUE,col=coloor,main="sKDR by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
barplot(t(yy),beside=TRUE,col=coloor[c(2:3)],main="MACE by host")
legend(legend=c("RS","SS"),fill=coloor[c(2:3)],x="topright")
barplot(t(ww),beside=TRUE,col=coloor,main="R81T by host")
legend(legend=c("RR","RS","SS"),fill=coloor,x="topright")
par(op)
#export to .pdf 8 x 7 inches


##############################################################################/
#END
##############################################################################/