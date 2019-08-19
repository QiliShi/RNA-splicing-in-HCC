### Main figures of RNA splicing in HCC - Qili Shi - 2019-03-25

#AS pie

pie.data <- as.data.frame(table(LIHC.PI.DS.data$splice_type))

mycolor<-c('#ff5f58','#ffeb6d','#ffae42','#ff7cb3','#00b2e3','#00fadc','#9fe98c')

library(extrafont)

tiff(filename = 'E:/LIHC_splicing/pie_LIHC.tif', res=600, compression = "lzw", height=8, width=8, units="in",pointsize = 24)

par(family='Times New Roman')

pie(pie.data$Freq,labels=pie.data$Freq,col=mycolor,border = NA)

dev.off()

# module bar plot


modules.splicing <- mapply(function(x,y){tmp <- as.data.frame(table(x));tmp[,2]<-tmp[,2]/sum(tmp[,2]);tmp$type <- y;return(tmp)},splicing.type,names(splicing.type),SIMPLIFY = F)

modules.splicing.plot <- do.call(rbind,modules.splicing)

modules.splicing.plot <-modules.splicing.plot[modules.splicing.plot$type!='grey',]

library(ggplot2)

library(extrafont)

modules.splicing.plot$type <- factor(modules.splicing.plot$type,levels = names(splicing.type)[order(names(splicing.type),decreasing = T)],ordered = T)

tiff("E:/LIHC_splicing/modules_bar.tif", res=600, compression = "lzw", height=9, width=10, units="in")

bp<- ggplot(modules.splicing.plot , aes(x=type, y=Freq*100, fill=x))+
  geom_bar(stat = "identity")+scale_fill_manual(values=c('#ff5f58','#ffeb6d','#ffae42','#ff7cb3','#00b2e3','#00fadc','#9fe98c'))
bp+theme_bw()+theme(text=element_text(family="Times New Roman",face = 'bold'),axis.title=element_text(size = 36),axis.text=element_text(size =28,colour = 'black'),
                    legend.text =element_text(size =28),legend.title=element_blank(),legend.key = element_rect(size = 8),
                    legend.key.size = unit(3, 'lines'))+labs(x ='',y='Percent(%)')+ coord_flip()
#+geom_text(aes(y=c(41/79*100,100),label=label), vjust=1.2, size=20)


dev.off()

## AT GO

modules.GOs.df.mito <- modules.GOs.df[grep('mitoch',modules.GOs.df$Description),]


modules.GOs.df.mito$module <- sapply(strsplit(row.names(modules.GOs.df.mito),"\\."),'[',1)

modules.GOs.df.mito <- modules.GOs.df.mito[!(modules.GOs.df.mito$module=='red'|modules.GOs.df.mito$module=='blue'|
                                               modules.GOs.df.mito$module=='grey'),]
modules.GOs.df.mito$GeneRatio <- sapply(modules.GOs.df.mito$GeneRatio ,function(x) eval(parse(text =x)))



#fix(modules.GOs.df.mito)

modules.GOs.df.mito <-modules.GOs.df.mito[rev(order(modules.GOs.df.mito$Description)),]

modules.GOs.df.mito$Description <- factor(modules.GOs.df.mito$Description,
                                          levels = unique(modules.GOs.df.mito$Description),ordered = T)

modules.GOs.df.mito$p.adjust <-(-log10(as.numeric(modules.GOs.df.mito$p.adjust )))



library(ggplot2)

library(extrafont)

p <- ggplot(modules.GOs.df.mito, aes_(x = ~module, y = ~Description, 
                               size = ~GeneRatio))+geom_point() + aes_string(color ='p.adjust') + scale_size_continuous(range = c(4,10))+
  scale_colour_gradient(low = "#80e1f1", high ="#ff6362" ,name='-log10(p)')+theme_bw()+labs(y = "",x="")

tiff("E:/LIHC_splicing/AT_modules_GO.tif", res=600, compression = "lzw", height=8, width=13, units="in")

p+theme(text=element_text(family="Times New Roman",face = 'bold'),axis.text.y=element_text(size =18,colour = 'black'),axis.text.x = element_text(colour = 'black',angle=45, size=20,hjust = 1),
                    axis.title.x=element_text(size =20),legend.text =element_text(size =20),legend.title =element_text(size =20))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))
dev.off()




# RBP counts

library(ggplot2)

RBP.counts <- apply(LIHC.RBP.p.value,2,function(x) sum(x<0.05))

RBP.counts <- data.frame(counts=RBP.counts,modules=substring(names(RBP.counts),3))

RBP.counts <- RBP.counts[RBP.counts$modules!='grey',]
tiff("E:/LIHC_splicing//RBP_signiticant_new.tif", res=600, compression = "lzw", height=8, width=10, units="in")

p<- ggplot(RBP.counts,aes(x=modules, y=counts))+
  geom_bar(stat = "identity",fill='#ff7473',width = 0.75)
p+theme_bw()+theme(text=element_text(family="Times New Roman",face = 'bold'),
                   axis.title=element_text(size = 28),axis.text=element_text(size =28,colour = 'black'),axis.text.x=element_text(angle=45, size=24,hjust = 1))+
  labs(x ='',y='Number of RBPs')
#  geom_text(aes(y=Percent*100+3,label =round(Percent*100,1)),size=7)
#+geom_text(aes(y=c(41/79*100,100),label=label), vjust=1.2, size=20)

dev.off()


RBP.GOs.df <- do.call(rbind,RBP.GOs)

RBP.GOs.table <- as.data.frame(table(RBP.GOs.df$Description))

choose.GO <- c('GO:0008380','GO:0050658','GO:0031123','GO:0009451','GO:0043631','GO:0001510','GO:0000966','GO:0000184','GO:0000288','GO:0006368')

RBP.GOs.plot <- RBP.GOs.df[RBP.GOs.df$ID%in%choose.GO,]

RBP.GOs.plot$modules <- gsub('ME(.+)\\.GO:.+','\\1',row.names(RBP.GOs.plot))

RBP.GOs.plot$Description <-gsub('nuclear-transcribed mRNA catabolic process,(.+)','\\1',RBP.GOs.plot$Description)

RBP.GOs.plot$GeneRatio <- sapply(RBP.GOs.plot$GeneRatio ,function(x) eval(parse(text =x)))

RBP.GOs.plot$Description[RBP.GOs.plot$ID=='GO:0006368'] <- 'transcription elongation from\n RNA polymerase II promoter'

RBP.GOs.plot <- RBP.GOs.plot[RBP.GOs.plot$modules!='grey',]

RBP.GOs.plot$p.adjust <- -log10(RBP.GOs.plot$p.adjust)

library(ggplot2)

library(extrafont)

p <- ggplot(RBP.GOs.plot, aes_(x = ~modules, y = ~Description, 
                              size = ~GeneRatio))+geom_point() + aes_string(color ='p.adjust') + scale_size_continuous(range = c(4,10))+
  scale_colour_gradient(low = "#80e1f1", high = "#ff6362",name='-log10(p)')+theme_bw()+labs(y = "")

tiff("E:/LIHC_splicing//RBP_modules_GO2.tif", res=600, compression = "lzw", height=8, width=13, units="in")

p+labs(x ='')+theme(text=element_text(family="Times New Roman",face = 'bold'),axis.text.y=element_text(size =18,colour = 'black'),axis.text.x = element_text(colour = 'black',angle=45, size=20,hjust = 1),
                    axis.title.x=element_text(size =20),legend.text =element_text(size =20),legend.title =element_text(size =20))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))
dev.off()

#module figure 

library(WGCNA)

tiff("E:/LIHC_splicing/module_cor2.tif", res=600, compression = "lzw", height=8, width=4.5,units="in")

# Will display correlations and their p-values

opar<-par(no.readonly = T)

par(mar=c(6.5, 10.4, 2.7, 1)+0.3)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(LIHC.RBP.cor[1:5,]),
               colorLabels = FALSE,
               xLabels = row.names(LIHC.RBP.cor)[1:5],
               yLabels = colnames(LIHC.RBP.cor),
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab = 1.3,
               zlim = c(-1,1))
dev.off()

par <- opar


# autoregulation

# box

tiff(filename = 'E:/LIHC_splicing/PI_DEGS_RBPs_new.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)

p <- list('ARSP-RBPs'=rowMeans(log2(LIHCMatrix.RBP[row.names(LIHCMatrix.RBP)%in%all.genes[all.genes%in%row.names(LIHC.DEGs.RBP)],]+1)),
          'AR-RBPs'=rowMeans(log2(LIHCMatrix.RBP[row.names(LIHCMatrix.RBP)%in%LIHC.DEGs.RBP.PI,]+1)),
          'NAR-RBPs'=rowMeans(log2(LIHCMatrix.RBP[!row.names(LIHCMatrix.RBP)%in%LIHC.DEGs.RBP.PI,]+1)))


boxplot(p,col=c("#ff6362","#80e1f1","#80debb"),outline=F,boxwex=0.5,ylim=c(10.5,30),ylab='mRNA expression')

dev.off()

#bar

library(ggplot2)

library(extrafont)

RBP.bar <- LIHC.DEGs.RBP[row.names(LIHC.DEGs.RBP)%in%LIHC.DEGs.RBP.PI&LIHC.DEGs.RBP$logFC>0,]


RBP.bar <- RBP.bar[order(RBP.bar$logFC,decreasing = T)[1:10],]

RBP.bar$symbol <- factor(row.names(RBP.bar),levels = row.names(RBP.bar),order=T)

RBP.bar$q <- -log10(RBP.bar$P.ajust)
 
tiff("E:/LIHC_splicing/RBP_overlap.tif", res=600, compression = "lzw", height=8, width=8, units="in")

p<- ggplot(RBP.bar,aes(x=symbol,y=logFC,fill=q))+
  geom_bar(stat = "identity",width = 0.75)+scale_fill_gradient(low="#80e1f1",high='#ff6362',name='-log10(p)')

  
p+theme_minimal()+theme(text=element_text(family="Times New Roman",face = 'bold'),
                   axis.title=element_text(size = 28),axis.text=element_text(size =28,colour = 'black'),axis.text.x=element_text(angle=45, size=24,hjust = 1),
                   legend.text =element_text(size =20),legend.title =element_text(size =20))+
  labs(x ='',y='log2(FC)')
#  geom_text(aes(y=Percent*100+3,label =round(Percent*100,1)),size=7)
#+geom_text(aes(y=c(41/79*100,100),label=label), vjust=1.2, size=20)

dev.off()



#venn

library(VennDiagram)

venn.diagram(list(row.names(LIHC.DEGs.RBP),LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL)]),
                  filename = "E:/LIHC_splicing/RBP_DASEs_DEGs_new.tif",fill =c("#80e1f1",'#ff6362'),cex=0,lwd=2,resolution = 800,category.names = c(1,1))

library(VennDiagram)

venn.diagram(list(RBP.self,LIHC.DEGs.RBP.PI),
             filename = "E:/LIHC_splicing/AR-RBP_eclip_new.tif",fill =c("#8484FB", "#CBC9FB"),cex=0,lwd=2,resolution = 800,category = rep("", 2))


#autoregulation survival

tiff(filename = 'E:/LIHC_splicing/PI_sur_RBPs.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)

p <- list('RBPs-Sig-onco'=rowMeans(LIHC.RBP.sig.onco),
          'RBPs-Sig'=rowMeans(LIHC.RBP.sig.exp[LIHC.RBP.anti,]),
          'RBPs-NS'=rowMeans(LIHC.RBP.sig.exp[!row.names(LIHCMatrix.RBP)%in%LIHC.RBP.anti,]))


boxplot(p,col=c("#ff6362","#80e1f1","#80debb"),outline=F,boxwex=0.5,ylim=c(10,30),ylab='mRNA expression')

dev.off()


stage




venn.diagram(list(eCLIP=RBP.self,'RBP-Sig'=LIHC.RBP.sig,'AS-Sig'=LIHC.PI.RBP.sig),filename = "E:/LIHC_splicing/SURRBPVenn.tif")


# hUb MMS

tiff(filename = 'E:/LIHC_splicing/MMS_sur_RBPs.tif', res=600, compression = "lzw", height=10, width=8, units="in")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)


LIHC.GMM.list <- list('RBPs-AS-Sig'=abs(LIHC.GMM$MMS[LIHC.PI.pvalue < 0.05&LIHC.GMM$module!='grey'&LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL)]),
  'AS-Sig'=abs(LIHC.GMM$MMS[LIHC.PI.pvalue < 0.05&LIHC.GMM$module!='grey']),
  'AS-NS'=abs(LIHC.GMM$MMS[LIHC.PI.pvalue>=0.05&LIHC.GMM$module!='grey']))

boxplot(LIHC.GMM.list,col=c("#ff6362","#80e1f1","#80debb"),outline=F,boxwex=0.5,ylim=c(0.25,1.15),ylab='Module membership')

dev.off()

#cytoscape

library(WGCNA)


MEs0 <- moduleEigengenes(t(LIHC.PI.DS), mergedColors.LIHC)$eigengenes

MEs <- orderMEs(MEs0)

TOM <- TOMsimilarityFromExpr(t(LIHC.PI.DS), power = 3,networkType = "unsigned", TOMType = "unsigned")

modules <-  c('green')

inModule <- is.finite(match(mergedColors.LIHC, modules))

modgenes <- LIHC.PI.DS.data$as_id[inModule]

modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modgenes, modgenes)

library(WGCNA)

cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste('E:/LIHC_splicing/cytoscape_edge.txt', sep=""),
                                nodeFile = paste('E:/LIHC_splicing/cytoscape_node.txt', sep=""),
                                weighted = TRUE,
                                threshold =0.05,
                                nodeNames = modgenes,
                                nodeAttr = mergedColors.LIHC[inModule]);
edge <-cyt$edgeData

node <- cyt$nodeData

edge <- edge[order(abs(edge$weight),decreasing = T)[1:(0.05*nrow(edge))],]

node <- node[node$nodeName%in%c(edge$fromNode,edge$toNode),]


node$symbol <-LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$as_id%in%node$nodeName]


library(igraph)

g.green <- graph.data.frame(cyt$edgeData,directed="FALSE")

node$strength <- graph.strength(g.green)[as.character(node$nodeName)]

node$RBP <- ifelse(node$symbol%in%c(RBP$SYMBOL,RDBP),'Y','N')

node$MMS <- abs(LIHC.GMM[LIHC.PI.DS.data$as_id%in%node$nodeName,1])

node$type <- 'NS'

node$type[node$nodeName%in%LIHC.PI.DS.data$as_id[LIHC.PI.pvalue<0.05]] <- 'significant'


node$type[node$nodeName%in%LIHC.PI.DS.data$as_id[LIHC.PI.pvalue<0.05&LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL)]] <-'RBP-Sig'

write.table(edge,'E:/LIHC_splicing/cytoscape_green_edge.txt',quote = F,row.names = F,sep='\t')

write.table(node,'E:/LIHC_splicing/cytoscape_green_node.txt',quote = F,row.names = F,sep='\t')

# RNA splicing modules

modules <- 'blue'

inModule <- mergedColors.LIHC=='blue'& LIHC.PI.DS.data$symbol%in%blue.genes

modgenes <- LIHC.PI.DS.data$as_id[inModule]

modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modgenes, modgenes)


cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste('E:/LIHC_splicing/cytoscape_edge_blue.txt', sep=""),
                                nodeFile = paste('E:/LIHC_splicing/cytoscape_node_blue.txt', sep=""),
                                weighted = TRUE,
                                threshold =0.02,
                                nodeNames = modgenes,
                                nodeAttr = mergedColors.LIHC[inModule]);

node <- cyt$nodeData

node$symbol <-LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$as_id%in%node$nodeName]


library(igraph)

g.blue <- graph.data.frame(cyt$edgeData,directed="FALSE")

node$strength <- graph.strength(g.blue)[as.character(node$nodeName)]

node$RBP <- ifelse(node$symbol%in%c(RBP$SYMBOL,RDBP),'Y','N')

write.table(node,'E:/LIHC_splicing/cytoscape_node_blue.txt',quote = F,row.names = F,sep='\t')


#magenta

modules <- 'magenta'

inModule <- mergedColors.LIHC=='magenta'& LIHC.PI.DS.data$symbol%in%magenta.genes

modgenes <- LIHC.PI.DS.data$as_id[inModule]

modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modgenes, modgenes)


cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste('E:/LIHC_splicing/cytoscape_edge_magenta.txt', sep=""),
                                nodeFile = paste('E:/LIHC_splicing/cytoscape_node_magenta.txt', sep=""),
                                weighted = TRUE,
                                threshold =0.02,
                                nodeNames = modgenes,
                                nodeAttr = mergedColors.LIHC[inModule]);

node <- cyt$nodeData

node$symbol <-LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$as_id%in%node$nodeName]


library(igraph)

g.magenta <- graph.data.frame(cyt$edgeData,directed="FALSE")

node$strength <- graph.strength(g.magenta)[as.character(node$nodeName)]

node$RBP <- ifelse(node$symbol%in%c(RBP$SYMBOL,RDBP),'Y','N')

write.table(node,'E:/LIHC_splicing/cytoscape_node_magenta.txt',quote = F,row.names = F,sep='\t')

#tan

modules <- 'tan'

inModule <- mergedColors.LIHC=='tan'& LIHC.PI.DS.data$symbol%in%tan.genes

modgenes <- LIHC.PI.DS.data$as_id[inModule]

modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modgenes, modgenes)


cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste('E:/LIHC_splicing/cytoscape_edge_tan.txt', sep=""),
                                nodeFile = paste('E:/LIHC_splicing/cytoscape_node_tan.txt', sep=""),
                                weighted = TRUE,
                                threshold =0.02,
                                nodeNames = modgenes,
                                nodeAttr = mergedColors.LIHC[inModule]);

node <- cyt$nodeData

node$symbol <-LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$as_id%in%node$nodeName]


library(igraph)

g.tan <- graph.data.frame(cyt$edgeData,directed="FALSE")

node$strength <- graph.strength(g.tan)[as.character(node$nodeName)]

node$RBP <- ifelse(node$symbol%in%c(RBP$SYMBOL,RDBP),'Y','N')

write.table(node,'E:/LIHC_splicing/cytoscape_node_tan.txt',quote = F,row.names = F,sep='\t')

# green GO
green <- modules.GOs$green
all.genes <- strsplit(green$geneID,'/')

hub.rbp <- unique(node$symbol[node$type=='RBP-Sig'&node$MMS>0.9])

green$'HubRatio' <- sapply(all.genes,function(x) sum(hub.rbp%in%x))

green$'HubRatio' <- green$'HubRatio'/green$Count

go.plot <- green[sapply(all.genes,function(x) sum(hub.rbp%in%x)>5),]

go.plot$Description <-gsub('nuclear-transcribed mRNA catabolic process,(.+)','\\1',go.plot$Description)

go.plot$GeneRatio <- sapply(go.plot$GeneRatio ,function(x) eval(parse(text =x)))

go.plot <- go.plot[c(3,6,7,8,10:13,17),]

go.plot$Description <- factor(go.plot$Description,levels = rev(go.plot$Description),ordered = T)

go.plot$p.adjust <-  -log10(go.plot$p.adjust)


library(ggplot2)

library(extrafont)

go.plot$p.adjust <- signif(go.plot$p.adjust,1)

p <- ggplot(go.plot , aes(x = HubRatio, y = Description, 
                               size = GeneRatio,color=p.adjust))+geom_point() + scale_size_continuous(range = c(4,10))+
  scale_colour_gradient(low = "#80e1f1", high = "#ff6362",name='-log10(p)')+theme_bw()+labs(y='')

tiff("E:/LIHC_splicing/green_GO.tif", res=600, compression = "lzw", height=10, width=10, units="in")

p+theme(text=element_text(family="Times New Roman",face = 'bold'),axis.text.y=element_text(size =20,colour = 'black'),axis.text.x = element_text(colour = 'black',angle=45, size=20,hjust = 1),
                    axis.title=element_text(size =24),legend.text =element_text(size =20),legend.title =element_text(size =20))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))+scale_x_continuous(limits = c(0.1, 0.18))
dev.off()

# RBP expression survival

library(ggplot2)

library(extrafont)


tmp <- data.frame('Percent'=c(231/302,1-231/302),group=c('Anti RBPs-sig','All RBPs-sig'),label=c(231,302))

tiff("E:/LIHC_splicing/RBPbar.tif", res=600, compression = "lzw", height=8, width=6.5, units="in")

bp<- ggplot(tmp, aes(x="", y=Percent*100, fill=group))+
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values=c("#80e1f1","#ff6362"))
bp+theme_bw()+theme(text=element_text(family="Times New Roman",face = 'bold'),axis.title=element_text(size = 36),axis.text=element_text(size =36,colour = 'black'),
                    legend.text =element_text(size =28),legend.title=element_blank(),legend.key = element_rect(size = 8),
                    legend.key.size = unit(3, 'lines'))+labs(x ='',y='Percent(%)')+
  geom_text(aes(y=c(231/302*100,100),label=label), vjust=1.2, size=20)

dev.off()

library("survminer")

library(survival)

df <- data.frame(OS=LIHC.clinical$OS,groups=as.character(ifelse(t(LIHCMatrix.RBP['SOX11',])>
                                                     median(t(LIHCMatrix.RBP['SOX11',])),'high','low')),
                          status=LIHC.clinical$vital_status=='dead',stringsAsFactors = F)

SOX11.fit <- survfit(Surv(as.numeric(OS)/30,status)~groups,data =df )


tiff("E:/LIHC_splicing/SOX11_survial.tif", res=600, compression = "lzw", height=10, width=10, units="in")

ggsurvplot(SOX11.fit,data=df ,xlab = "Time in months ",
           legend.title = "SOX11",legend.labs = c("high", "low"),pval = T,
           conf.int = F,palette = c('#e7220e','#354c9f'),
           ggtheme = theme_survminer(base_family = 'Times New Roman',font.x =28,font.y =28,font.legend = 28,font.tickslab=c(24,"plain", "black")))
graph2ppt(file='E:/LIHC_splicing/SOX11_survial.pptx',width=12,height=12)
dev.off()

df <- data.frame(OS=LIHC.clinical$OS,groups=as.character(ifelse(t(LIHCMatrix.RBP['KPNA2',])>
                                                                  median(t(LIHCMatrix.RBP['KPNA2',])),'high','low')),
                 status=LIHC.clinical$vital_status=='dead',stringsAsFactors = F)

KPNA2.fit <- survfit(Surv(as.numeric(OS)/30,status)~groups,data =df )


tiff("E:/LIHC_splicing/KPNA2_survial.tif", res=600, compression = "lzw", height=10, width=10, units="in")

ggsurvplot(KPNA2.fit,data=df ,xlab = "Time in months ",
           legend.title = "KPNA2",legend.labs = c("high", "low"),
           conf.int = F,palette = c('#e7220e','#354c9f'),
           ggtheme = theme_survminer(base_family = 'Times New Roman',font.x =28,font.y =28,font.legend = 28,font.tickslab=c(24,"plain", "black")))
graph2ppt(file='E:/LIHC_splicing/KPNA2_survial.pptx',width=12,height=12)
dev.off()


library(VennDiagram)

venn.diagram(list(LIHC.RBP.anti,row.names(LIHC.RBP.DEG)),
             filename = "E:/LIHC_splicing/RBP_sig.tif",cex=0,category = c("155", '140'))

#stage box

tiff(filename = 'E:/LIHC_splicing/stages_RBPs.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)



boxplot(results.stage.sig,col=c("#80debb","#80e1f1","#ff6362"),outline=F,boxwex=0.5,ylim=c(4,17),ylab='mRNA expression')

legend("topleft", inset=.05,legend=c("Stage1","Stage2","Stage3"),
       fill=c("#ff6362","#80e1f1","#80debb"), cex=2,text.font=2,
       box.lty=0)

dev.off()




# survival

library("survminer")

library(extrafont)

clusters.df <- data.frame(OS=LIHC.clinical$OS,stage=LIHC.clinical$tumor_stage,Cluster=clusters,
                          status=LIHC.clinical$vital_status=='dead',stringsAsFactors = F)

fit <- survfit(Surv(as.numeric(OS)/30,status)~Cluster,data = clusters.df)

tiff("E:/LIHC_splicing/clusters_survial.tif", res=600, compression = "lzw", height=8, width=8, units="in")

ggsurvplot(fit,data=clusters.df,xlab = "Time in months",palette = c("#ff6362","#80e1f1","#80debb",'#ffae00'), pval = TRUE,
           ggtheme = theme_survminer(base_family = 'Times New Roman',font.x =20,font.y =20,font.legend = 20,font.tickslab=c(20,"plain", "black")))
dev.off()

#hub
  
clusters.hub.df <- data.frame(OS=LIHC.clinical$OS,stage=LIHC.clinical$tumor_stage,Cluster=clusters.hub,
                          status=LIHC.clinical$vital_status=='dead',stringsAsFactors = F)

fit <- survfit(Surv(as.numeric(OS)/30,status)~Cluster,data = clusters.hub.df)

tiff("E:/LIHC_splicing/clusters_survial_hub.tif", res=600, compression = "lzw", height=8, width=8, units="in")

ggsurvplot(fit,data=clusters.hub.df,xlab = "Time in months",palette = c("#ff7cb3","#80e1f1","#80debb",'#ffae00',"#00b2e3","#ff6362"), pval = TRUE,
           ggtheme = theme_survminer(base_family = 'Times New Roman',font.x =20,font.y =20,font.legend = 20,font.tickslab=c(20,"plain", "black")))
dev.off()

#box

LIHC.RBP <-log2(LIHCMatrix.RBP+1)

results <- sapply(unique(clusters.hub),function(x) colMeans(LIHC.RBP[,names(clusters.hub)[clusters.hub==x]]),simplify = F)

tiff(filename = 'E:/LIHC_splicing/clusters_RBPs.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)


boxplot(results,col=c("#ff7cb3","#80e1f1","#80debb",'#ffae00',"#00b2e3","#ff6362"),outline=F,boxwex=0.5,ylim=c(16.5,17.8),ylab='RBPs mRNA expression')

graph2ppt(file='E:/LIHC_splicing/clusters_RBPs.pptx',width=8,height=10)

dev.off()

#

LIHC.RBP <-log2(LIHCMatrix.RBP+1)

results <- sapply(unique(clusters.hub),function(x) colMeans(log2(LIHCMatrix.RBP[,names(clusters.hub)[clusters.hub==x]]+1)),simplify = F)

tiff(filename = 'E:/LIHC_splicing/clusters_RBPs-DEG.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)


boxplot(results,col=c("#ff7cb3","#80e1f1","#80debb",'#ffae00',"#00b2e3","#ff6362"),outline=F,boxwex=0.5,ylim=c(16.5,17.5),
        ylab='RBPs-DEG mRNA expression')

dev.off()

#

LIHC.RBP <-log2(LIHCMatrix.RBP+1)

results <- sapply(unique(clusters.hub),function(x) colMeans(LIHC.RBP[LIHC.RBP.sig,names(clusters.hub)[clusters.hub==x]]),simplify = F)

tiff(filename = 'E:/LIHC_splicing/clusters_RBPs-sig.tif', res=600, compression = "lzw", height=10, width=8, units="in")

# windowsFonts(A = windowsFont("Times New Roman"))

# par(family="A")

par(cex.lab=2,cex.axis=2,font.axis=2,font.lab=2,lwd=3,family='Times New Roman',mar=c(7, 5, 2.7, 1)+0.3)


boxplot(results,col=c("#ff7cb3","#80e1f1","#80debb",'#ffae00',"#00b2e3","#ff6362"),outline=F,boxwex=0.5,ylim=c(7,11),ylab='RBPs-DEG mRNA expression')

dev.off()

#subtype 6 heatmap

RBP.subtype6.heatmap <- t(apply(LIHCMatrix.RBP[RBP.subtype6,],1,scale))


RBP.subtype6.heatmap <-RBP.subtype6.heatmap[,order(clusters.hub,decreasing = F)]

cbPalette <- c("#ff7cb3","#80e1f1","#80debb",'#ffae00',"#00b2e3","#ff6362")

table(clusters.hub.df$Cluster)

colors=c(rep(cbPalette[1],40),rep(cbPalette[2],66),rep(cbPalette[3],145),rep(cbPalette[4],55),rep(cbPalette[5],22),rep(cbPalette[6],43))


library(gplots)
tiff("E:/LIHC_splicing/heatmap_colorkey.tif", res=600, compression = "lzw", height=10, width=16, units="in")

#pdf("E:/LIHC_splicing/heatmap3.pdf", height=16, width=10)

heatmap.2(RBP.subtype6.heatmap,Colv=F,col=greenred(75),scale=NULL,margins = c(11,13),
          ColSideColors=colors,key=TRUE,keysize=1,symkey=FALSE, density.info="none", trace="none", cexRow=2.2,cexCol =0.00001,font=2)
dev.off()

graph2ppt(file='E:/LIHC_splicing/heatmap.pptx',width=10,height=16)

P <- heatmap.2(RBP.subtype6.heatmap,Colv=F,col=greenred(75),scale=NULL,
               ColSideColors=colors,key=F,symkey=FALSE, density.info="none", trace="none", cexRow=0.1,cexCol =0.00001)

row.names(RBP.subtype6.heatmap)[rev(P$rowInd)]

# RBPs events
library(ggplot2)

library(extrafont)



RBPs_hub_RBPs <- data.frame(percent=100*rowSums(t(rbp.events.pvalue<0.05))/94,genes=factor(colnames(rbp.events.pvalue),
                            levels = rev(row.names(RBP.subtype6.heatmap)[rev(P$rowInd)]),ordered = T))

#RBPs_hub_RBPs <- RBPs_hub_RBPs[RBPs_hub_RBPs$percent>50,]
tiff("E:/LIHC_splicing/RBPs_RBPs_new.tif", res=600, compression = "lzw", height=10, width=8, units="in")

p<- ggplot(RBPs_hub_RBPs ,aes(x=genes,y=percent,fill='#ff6362'))+
  geom_bar(stat = "identity",width = 0.75)


p+theme_minimal()+theme(text=element_text(family="Times New Roman",face = 'bold'),
                        axis.title=element_text(size = 16),axis.text=element_text(size=16,colour = 'black'),axis.text.x=element_text(angle=45, size=20,hjust = 1,colour = 'black'),
                        legend.text =element_text(size =20),legend.title =element_text(size =24))+
  labs(x ='',y='Percent(%)')+coord_flip()
#  geom_text(aes(y=Percent*100+3,label =round(Percent*100,1)),size=7)
#+geom_text(aes(y=c(41/79*100,100),label=label), vjust=1.2, size=20)

dev.off()

rev(P$rowInd)


# clip track 

auto <- intersect(intersect(row.names(LIHC.DEGs.RBP),RBP.self),intersect(RBP.self,LIHC.PI.DS.data$symbol))

ILF3.AS <- LIHC.PI.DS.data[LIHC.PI.DS.data$symbol=='ILF3',]

