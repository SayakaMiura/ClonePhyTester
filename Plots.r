library(reshape2)
library(ggplot2)
require(gdata)

jpeg('plot.jpg')
#par(mfrow=c(3,1))
MetLs=c("YourMethod","CloneFinder","MACHINA","TreeOmics",  "LICHeE", "MixPhy", "PhyloWGS",  "Cloe")
#C:\\Users\\tuf78332\\Documents\\GitHub\\ClonePhyTester\\
Seq00 = read.xls("All.xlsx", sheet = 'Sequential_MutPair_count', header = TRUE)
Seq0=Seq00[is.na(Seq00$YourMethod)!=TRUE,]
Len=length(Seq0$Sim)
Seq=cbind(Seq0,rep('SeqMutOrder',Len))
colnames(Seq)[colnames(Seq)=="rep(\"SeqMutOrder\", Len)"] <- "ErrorScore"
Seq1=melt(data = Seq, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

Par00 = read.xls("All.xlsx", sheet = 'parallel_MutPair_count', header = TRUE)
Par0=Par00[is.na(Par00$YourMethod)!=TRUE,]
Len=length(Par0$Sim)
Par=cbind(Par0,rep('ParaMutOrder',Len))
colnames(Par)[colnames(Par)=="rep(\"ParaMutOrder\", Len)"] <- "ErrorScore"
Par1=melt(data = Par, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

Con00 = read.xls("All.xlsx", sheet = 'concurrent_MutPair_count', header = TRUE)
Con0=Con00[is.na(Con00$YourMethod)!=TRUE,]
Len=length(Con0$Sim)
Con=cbind(Con0,rep('ConcuMutOrder',Len))
colnames(Con)[colnames(Con)=="rep(\"ConcuMutOrder\", Len)"] <- "ErrorScore"
Con1=melt(data = Con, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

ML00 = read.xls("All.xlsx", sheet = 'MLTED', header = TRUE)
ML0=ML00[is.na(ML00$YourMethod)!=TRUE,]
Len=length(ML0$Sim)
ML=cbind(ML0,rep('MLTED',Len))
colnames(ML)[colnames(ML)=="rep(\"MLTED\", Len)"] <- "ErrorScore"
ML1=melt(data = ML, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

TV00 = read.xls("All.xlsx", sheet = 'TreeVec', header = TRUE)
TV0=TV00[is.na(TV00$YourMethod)!=TRUE,]
Len=length(TV0$Sim)
TV=cbind(TV0,rep('TreeVec',Len))
colnames(TV)[colnames(TV)=="rep(\"TreeVec\", Len)"] <- "ErrorScore"
TV1=melt(data = TV, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

RF00 = read.xls("All.xlsx", sheet = 'RF', header = TRUE)
RF0=RF00[is.na(RF00$YourMethod)!=TRUE,]
Len=length(RF0$Sim)
RF=cbind(RF0,rep('RF',Len))
colnames(RF)[colnames(RF)=="rep(\"RF\", Len)"] <- "ErrorScore"
RF1=melt(data = RF, id.vars = c("Sim","Data","ErrorScore"), measure.vars = MetLs)

Ta=rbind(Seq1,Par1,Con1,ML1,TV1,RF1)
colnames(Ta)[colnames(Ta)=="variable"] <- "Method"
colnames(Ta)[colnames(Ta)=="value"] <- "ErrorRate"
Ta$Sim <- factor(Ta$Sim, levels = c("G7", "G12","P10",  "MA","MA50","G7cna","TGlinear","TGstep","TGconst"))
Ta$Method <- factor(Ta$Method, levels = MetLs)
ggplot(Ta, aes(x=Method, y=as.numeric(ErrorRate), color=Method)) +geom_boxplot()+facet_grid(ErrorScore~Sim, scales = "free_y")+ scale_x_discrete(limits=MetLs)+scale_color_manual(values=c("black","#6b9b87", "#5871ba", "#91c452", "#7f57b0", "#d0ba4a", "#e4616e", "#d29a44"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")

dev.off()
