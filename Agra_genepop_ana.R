##############################################################################/
##############################################################################/
#DAS tree on the clone corrected data
##############################################################################/
##############################################################################/

source("Agra_load.R")


##############################################################################/
#Formatting the dataset for network analysis####
##############################################################################/


onebyhost<-"data/onebyhost.txt"

basic_info(onebyhost,"output/onebyhost.txt.INF")
test_LD(onebyhost,"output/onebyhost.txt.DIS")
clean_workdir()

#the output file have to be edited so it can be used for representation
LDbyhost<-readLines("output/onebyhost.txt.DIS")
#removing the first lines that are non relevant
LDcolnames<-unlist(strsplit(LDbyhost[13],"\\s+"))
LDbyhost<-LDbyhost[-c(1:14)]
#for 14 markers and 4 populations there are 
14*13/2*4 #LD tests for markers pairs
LDbyhost<-LDbyhost[-c(365:length(LDbyhost))]
LDbyhost<-as.data.table(matrix(unlist(strsplit(LDbyhost,"\\s+")),
                               nrow=length(strsplit(LDbyhost,"\\s+")),
                               byrow=TRUE))
colnames(LDbyhost)<-LDcolnames
LDbyhost$`P-Value`<-as.numeric(LDbyhost$`P-Value`)
LDpeach<-as.matrix(LDbyhost[c(1:91),])
row.names(LDpeach)<-paste(LDpeach[,2],LDpeach[,3])
LDpeach<-rbind(LDpeach[,c(2,3,4)],LDpeach[,c(3,2,4)])
LDpeach<-spread(data.frame(LDpeach),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDpeach)<-as.character(LDpeach[,1])
#removing the first unnecessary column
LDpeach<-LDpeach[,c(-1)]
#reordering the columns and turning the object into a matrix
LDpeach<-as.matrix(LDpeach[c(2,8,12,13,14,1,3:7,9:11),
                        c(2,8,12,13,14,1,3:7,9:11)])
#scaling the p-value so it is easily usable with the LDheatm

temp<-LDpeach
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.05 & temp<0.8]<-0.6
temp[temp>0.01 & temp<0.6]<-0.4
temp[temp>0.001 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
chaudemap<-LDheatmap(temp,title=NULL,
                     add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                     color=c(rep(grey(0.8),3),
                             brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                     name="CHR",flip=TRUE,add.key=FALSE)
grid.edit(gPath("CHR","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("CHR","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("B",x=unit(0.048,"npc"),y=unit(0.955,"npc"),gp=gpar(fontsize=50),
          check=TRUE)
grid.text("Peach population",x=unit(0.5,"npc"),y=unit(0.75,"npc"),gp=gpar(fontsize=40),
          check=TRUE)
#export to pdf 8 x 8 inches



##############################################################################/
#END
##############################################################################/