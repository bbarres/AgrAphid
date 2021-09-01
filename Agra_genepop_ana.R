##############################################################################/
##############################################################################/
#Classic population genetic computations
##############################################################################/
##############################################################################/

source("Agra_load.R")
coloor<-c("firebrick","royalblue4","chartreuse4","khaki2","darkorange")

##############################################################################/
#Linkage disequilibrium in Host populations: computing and formatting####
##############################################################################/

onebyhost<-"data/onebyhost.txt"
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

#Linkage disequilibrium in the peach population
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

#Linkage disequilibrium in the oil seed population
LDoils<-as.matrix(LDbyhost[c(92:182),])
row.names(LDoils)<-paste(LDoils[,2],LDoils[,3])
LDoils<-rbind(LDoils[,c(2,3,4)],LDoils[,c(3,2,4)])
LDoils<-spread(data.frame(LDoils),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDoils)<-as.character(LDoils[,1])
#removing the first unnecessary column
LDoils<-LDoils[,c(-1)]
#reordering the columns and turning the object into a matrix
LDoils<-as.matrix(LDoils[c(2,8,12,13,14,1,3:7,9:11),
                         c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the tobaco population
LDtaba<-as.matrix(LDbyhost[c(183:273),])
row.names(LDtaba)<-paste(LDtaba[,2],LDtaba[,3])
LDtaba<-rbind(LDtaba[,c(2,3,4)],LDtaba[,c(3,2,4)])
LDtaba<-spread(data.frame(LDtaba),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDtaba)<-as.character(LDtaba[,1])
#removing the first unnecessary column
LDtaba<-LDtaba[,c(-1)]
#reordering the columns and turning the object into a matrix
LDtaba<-as.matrix(LDtaba[c(2,8,12,13,14,1,3:7,9:11),
                         c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the other host population
LDothe<-as.matrix(LDbyhost[c(274:364),])
row.names(LDothe)<-paste(LDothe[,2],LDothe[,3])
LDothe<-rbind(LDothe[,c(2,3,4)],LDothe[,c(3,2,4)])
LDothe<-spread(data.frame(LDothe),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDothe)<-as.character(LDothe[,1])
#removing the first unnecessary column
LDothe<-LDothe[,c(-1)]
#reordering the columns and turning the object into a matrix
LDothe<-as.matrix(LDothe[c(2,8,12,13,14,1,3:7,9:11),
                         c(2,8,12,13,14,1,3:7,9:11)])


##############################################################################/
#Figure S6: LD between markers within host populations####
##############################################################################/

#peach
temp<-LDpeach
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP1<-viewport(x=0,y=0.5,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp1")
pushViewport(VP1)
PEACH<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="peach",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("peach","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("peach","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Peach",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40),check=TRUE)
upViewport()

#oilseed rape
temp<-LDoils
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP2<-viewport(x=0.5,y=0.5,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp2")
pushViewport(VP2)
OILS<-LDheatmap(temp,title=NULL,
                add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                color=c(rep(grey(0.8),3),
                        brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                name="oils",flip=TRUE,add.key=FALSE,
                newpage=FALSE)
grid.edit(gPath("oils","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("oils","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Oilseed",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40),check=TRUE)
upViewport()

#tobaco
temp<-LDtaba
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP3<-viewport(x=0,y=0,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp3")
pushViewport(VP3)
TABA<-LDheatmap(temp,title=NULL,
                add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                color=c(rep(grey(0.8),3),
                        brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                name="taba",flip=TRUE,add.key=FALSE,
                newpage=FALSE)
grid.edit(gPath("taba","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("taba","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Tobacco",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40),check=TRUE)
upViewport()

#other hosts
temp<-LDothe
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP4<-viewport(x=0.5,y=0,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp4")
pushViewport(VP4)
OTHE<-LDheatmap(temp,title=NULL,
                add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                color=c(rep(grey(0.8),3),
                        brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                name="othe",flip=TRUE,add.key=FALSE,
                newpage=FALSE)
grid.edit(gPath("othe","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("othe","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Other hosts",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40),check=TRUE)
upViewport()

#export to pdf 15 x 15 inches


##############################################################################/
#Linkage disequilibrium in genetic clusters: computing and formatting####
##############################################################################/

oneclustrap4<-"data/oneclustrap4.txt"
test_LD(oneclustrap4,"output/oneclustrap4.txt.DIS")
clean_workdir()

#the output file have to be edited so it can be used for representation
LDclust4<-readLines("output/oneclustrap4.txt.DIS")
#removing the first lines that are non relevant
LDcolnames<-unlist(strsplit(LDclust4[13],"\\s+"))
LDclust4<-LDclust4[-c(1:14)]
#for 14 markers and 4 populations there are 
14*13/2*4 #LD tests for markers pairs
LDclust4<-LDclust4[-c(365:length(LDclust4))]
LDclust4<-as.data.table(matrix(unlist(strsplit(LDclust4,"\\s+")),
                               nrow=length(strsplit(LDclust4,"\\s+")),
                               byrow=TRUE))
colnames(LDclust4)<-LDcolnames
LDclust4$`P-Value`<-as.numeric(LDclust4$`P-Value`)

#Linkage disequilibrium in the red cluster population
LDcluRed<-as.matrix(LDclust4[c(1:91),])
row.names(LDcluRed)<-paste(LDcluRed[,2],LDcluRed[,3])
LDcluRed<-rbind(LDcluRed[,c(2,3,4)],LDcluRed[,c(3,2,4)])
LDcluRed<-spread(data.frame(LDcluRed),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDcluRed)<-as.character(LDcluRed[,1])
#removing the first unnecessary column
LDcluRed<-LDcluRed[,c(-1)]
#reordering the columns and turning the object into a matrix
LDcluRed<-as.matrix(LDcluRed[c(2,8,12,13,14,1,3:7,9:11),
                             c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the green cluster population
LDcluGre<-as.matrix(LDclust4[c(92:182),])
row.names(LDcluGre)<-paste(LDcluGre[,2],LDcluGre[,3])
LDcluGre<-rbind(LDcluGre[,c(2,3,4)],LDcluGre[,c(3,2,4)])
LDcluGre<-spread(data.frame(LDcluGre),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDcluGre)<-as.character(LDcluGre[,1])
#removing the first unnecessary column
LDcluGre<-LDcluGre[,c(-1)]
#reordering the columns and turning the object into a matrix
LDcluGre<-as.matrix(LDcluGre[c(2,8,12,13,14,1,3:7,9:11),
                             c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the yellow cluster population
LDcluYel<-as.matrix(LDclust4[c(183:273),])
row.names(LDcluYel)<-paste(LDcluYel[,2],LDcluYel[,3])
LDcluYel<-rbind(LDcluYel[,c(2,3,4)],LDcluYel[,c(3,2,4)])
LDcluYel<-spread(data.frame(LDcluYel),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDcluYel)<-as.character(LDcluYel[,1])
#removing the first unnecessary column
LDcluYel<-LDcluYel[,c(-1)]
#reordering the columns and turning the object into a matrix
LDcluYel<-as.matrix(LDcluYel[c(2,8,12,13,14,1,3:7,9:11),
                             c(2,8,12,13,14,1,3:7,9:11)])

#Linkage disequilibrium in the blue cluster population
LDcluBlu<-as.matrix(LDclust4[c(274:364),])
row.names(LDcluBlu)<-paste(LDcluBlu[,2],LDcluBlu[,3])
LDcluBlu<-rbind(LDcluBlu[,c(2,3,4)],LDcluBlu[,c(3,2,4)])
LDcluBlu<-spread(data.frame(LDcluBlu),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(LDcluBlu)<-as.character(LDcluBlu[,1])
#removing the first unnecessary column
LDcluBlu<-LDcluBlu[,c(-1)]
#reordering the columns and turning the object into a matrix
LDcluBlu<-as.matrix(LDcluBlu[c(2,8,12,13,14,1,3:7,9:11),
                             c(2,8,12,13,14,1,3:7,9:11)])


##############################################################################/
#Figure S7: LD between markers within genetic cluster populations####
##############################################################################/

#Red cluster
temp<-LDcluRed
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP1<-viewport(x=0,y=0.5,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp1")
pushViewport(VP1)
CLRED<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="clred",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("clred","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("clred","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Red Cluster",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40,col=coloor[1]),check=TRUE)
upViewport()

#oilseed rape
temp<-LDcluGre
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP2<-viewport(x=0.5,y=0.5,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp2")
pushViewport(VP2)
CLGRE<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="clgre",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("clgre","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("clgre","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Green Cluster",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40,col=coloor[3]),check=TRUE)
upViewport()

#tobaco
temp<-LDcluYel
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP3<-viewport(x=0,y=0,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp3")
pushViewport(VP3)
CLYEL<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="clyel",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("clyel","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("clyel","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Yellow Cluster",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40,col=coloor[4]),check=TRUE)
upViewport()

#other hosts
temp<-LDcluBlu
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.00055 & temp<0.8]<-0.6
temp[temp>0.00011 & temp<0.6]<-0.4
temp[temp>0.000011 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
VP4<-viewport(x=0.5,y=0,width=0.5,height=0.5,just=c("left","bottom"),
              name="vp4")
pushViewport(VP4)
CLBLU<-LDheatmap(temp,title=NULL,
                 add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                 color=c(rep(grey(0.8),3),
                         brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                 name="clblu",flip=TRUE,add.key=FALSE,
                 newpage=FALSE)
grid.edit(gPath("clblu","heatMap","heatmap"),gp=gpar(col="white",lwd=2.5))
grid.edit(gPath("clblu","SNPnames"),
          gp=gpar(col="black",rot="0",cex=1.2,font=2),
          rot=0,hjust=0.8)
grid.text("Blue Cluster",x=unit(0.5,"npc"),y=unit(0.75,"npc"),
          gp=gpar(fontsize=40,col=coloor[2]),check=TRUE)
upViewport()

#export to pdf 15 x 15 inches


##############################################################################/
#END
##############################################################################/