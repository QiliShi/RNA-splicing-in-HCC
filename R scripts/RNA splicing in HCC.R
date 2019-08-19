### RNA splicing in HCC - Qili Shi - 2019-03-06

options(stringsAsFactors = F)

LIHC.PI.data <-  read.delim('E:/LIHC_splicing/PSI_download_LIHC.txt',sep='\t',comment.char = "#",na.strings = 'null')

LIHC.sample <- read.delim('E:/LIHC_splicing/clinical.tsv',sep='\t',comment.char = "#")

LIHC.PI <- LIHC.PI.data[,-(1:10)]

LIHC.PI <- LIHC.PI[,colnames(LIHC.PI)!='X']

colnames(LIHC.PI) <- gsub('_','-',colnames(LIHC.PI))

groups <- ifelse(grepl('Norm',colnames(LIHC.PI)),'N','T')


groups <- factor(groups,levels = c('T','N'),ordered = T)

library(limma)

design <- model.matrix(~groups,ref='T')

fit <- lmFit(LIHC.PI,design)

fit <- eBayes(fit)

LIHC.PI.DS <- topTable(fit,number=nrow(LIHC.PI),coef = 2,adjust='BH',p.value = 0.05) ##top DE genes, number--max number / p.value

#DPI <- abs(rowMeans(LIHC.PI[,grepl('Norm',colnames(LIHC.PI))],na.rm = T)-rowMeans(LIHC.PI[,!grepl('Norm',colnames(LIHC.PI))],na.rm = T))


LIHC.PI.DS.data <- LIHC.PI.data[row.names(LIHC.PI.DS),] 

LIHC.PI.DS <- LIHC.PI[row.names(LIHC.PI.DS),]

library(impute)

LIHC.PI.DS <- impute.knn(as.matrix(LIHC.PI.DS))$data

LIHC.PI.DS <- LIHC.PI.DS[,!grepl('Norm',colnames(LIHC.PI))]

library(preprocessCore)

LIHC.PI.DS <- normalize.quantiles(LIHC.PI.DS)

colnames(LIHC.PI.DS) <- colnames(LIHC.PI)[!grepl('Norm',colnames(LIHC.PI))]

#LIHC.PI.DS.data <- t(apply(LIHC.PI.DS.data,1,scale))

#LIHC.PI.PCA <- princomp(LIHC.PI.DS.data)


# get WGCNA power

source('F:/myR/myscript/WGCNApower.R')

LIHC.power <- WGCNApower(LIHC.PI.DS)

LIHC.power <- 3

LIHC.net <- blockwiseModules(t(LIHC.PI.DS) , power = 3,
                             TOMType = "unsigned", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.1,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = "LIHC.PI",
                             verbose = 3)

save(LIHC.net,file='F:/myR/mydata/LIHC_splicing.RData')

library(WGCNA)

mergedColors.LIHC <- labels2colors(LIHC.net$colors)

modules <- unique(mergedColors.LIHC)

splicing.type <- sapply(modules, function(x) LIHC.PI.DS.data[mergedColors.LIHC==x,3])

sapply(splicing.type,table)

modules.genes <- sapply(modules, function(x) unique(LIHC.PI.DS.data[mergedColors.LIHC==x,1]))

file <- 'E:/LIHC_splicing/GO/'

modules.GOs <- lapply(modules.genes,GOPanel)

modules.GOs.df <-do.call(rbind,modules.GOs)


mapply(function(x,y) write.table(x,file=paste(file,y,'GO.txt',sep=''),quote = F,row.names = F,sep='\t'),modules.GOs,names(modules.GOs))

source('F:/myR/myscript/GOPanel.R')



library(extrafont)

tiff("E:/LIHC_splicing/splicing_module.tif", res=600, compression = "lzw", height=8, width=16, units="in",pointsize = 28)

par(family='Times New Roman',cex.lab=1.5,cex.axis=1.5,font.lab=2,lwd=2)

plotDendroAndColors(LIHC.net$dendrograms[[1]], mergedColors.LIHC [LIHC.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##down load LIHC expression data

library(TCGAbiolinks)

library(SummarizedExperiment)


query <- GDCquery(project = "TCGA-LIHC",##对应癌症
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")

GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query,summarizedExperiment = F)

data$X1 <-sapply(strsplit(data$X1,'\\.'),'[',1)

LIHCMatrix <- data[grepl('ENSG',data$X1),]


#datatable(as.data.frame(colData(data)), options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), rownames = FALSE)


ensemble <- read.delim('E:/LIHC_splicing/human_ensembl.txt')

LIHCMatrix <- merge(LIHCMatrix,ensemble,by.x='X1',by.y='Gene.stable.ID')



LIHCMatrix <- aggregate(LIHCMatrix[,grepl('TCGA',colnames(LIHCMatrix))],by=list(LIHCMatrix[,'Gene.name']),mean)

LIHCMatrix <- LIHCMatrix[,!grepl('-11A',colnames(LIHCMatrix))]

row.names(LIHCMatrix) <- LIHCMatrix$Group.1

colnames(LIHCMatrix) <- gsub('(.+)-0[12][AB]-.+','\\1',colnames(LIHCMatrix))

LIHCMatrix <- LIHCMatrix[,colnames(LIHC.PI.DS)]

LIHCMatrix <- LIHCMatrix[apply(LIHCMatrix,1,sd)!=0,]

LIHCMatrix.RBP <- LIHCMatrix[row.names(LIHCMatrix)%in%c(RBP$SYMBOL,RDBP),]


LIHC.MEs<- moduleEigengenes(t(LIHC.PI.DS), mergedColors.LIHC)$eigengenes

row.names(LIHC.MEs) <- colnames(LIHCMatrix)

library(WGCNA)

LIHC.RBP.cor <- t(cor(LIHC.MEs,t(LIHCMatrix.RBP),use = 'p'))

LIHC.RBP.cor <- LIHC.RBP.cor[,colnames(LIHC.RBP.cor)!='MEgrey']

LIHC.RBP.p.value <- corPvalueStudent(LIHC.RBP.cor,371)

LIHC.RBP.p.value <- apply(LIHC.RBP.p.value,2,function(x) p.adjust(x,method = 'BH'))

RBP.GOs <- apply(LIHC.RBP.p.value,2,function(x) data.frame(enrichGO(gene=row.names(LIHC.RBP.p.value)[x<0.05],OrgDb="org.Hs.eg.db",keyType="SYMBOL",ont ='BP')))

sapply(RBP.GOs, function(x) sum(grepl('polyadenylation',x$Description)))

# survival analysis

LIHC.clinical <- read.delim('E:/LIHC_splicing/clinical.tsv',stringsAsFactors = F)

LIHC.clinical$OS <- LIHC.clinical$days_to_last_follow_up

LIHC.clinical$OS[LIHC.clinical$vital_status=='dead'] <- LIHC.clinical$days_to_death[LIHC.clinical$vital_status=='dead'] 

LIHC.clinical <- LIHC.clinical[LIHC.clinical$submitter_id%in%colnames(LIHCMatrix),]

row.names(LIHC.clinical) <- LIHC.clinical$submitter_id

LIHC.clinical <- LIHC.clinical[colnames(LIHCMatrix),]


library(survival)


my.surv <- Surv(as.numeric(LIHC.clinical$OS),LIHC.clinical$vital_status=='dead')

LIHC.chisq <- sapply(colnames(LIHC.MEs),function(x) survdiff(my.surv~ifelse(LIHC.MEs[,x]>median(as.numeric(LIHC.MEs[,x])),'high','low'))$chisq)

LIHC.pvalue <- pchisq(LIHC.chisq, 1, lower.tail = FALSE)

LIHC.pvalue[LIHC.pvalue<0.05]

LIHC.PI.chisq <- apply(LIHC.PI.DS,1,function(x) survdiff(my.surv~ifelse(x>median(x),'high','low'))$chisq)

LIHC.PI.pvalue <- pchisq(LIHC.PI.chisq , 1, lower.tail = FALSE)

sum(LIHC.PI.pvalue<0.05)

#RBP survival
LIHC.RBP.chisq <- sapply(rownames(LIHCMatrix.RBP),function(x) survdiff(my.surv~t(ifelse(LIHCMatrix.RBP[x,]>median(as.numeric(LIHCMatrix.RBP[x,])),'high','low')))$chisq)

LIHC.RBP.pvalue <- pchisq(LIHC.RBP.chisq, 1, lower.tail = FALSE)

LIHC.RBP.risk<- sapply(rownames(LIHCMatrix.RBP),function(x){tmp<- survdiff(my.surv~t(ifelse(LIHCMatrix.RBP[x,]>median(as.numeric(LIHCMatrix.RBP[x,])),'high','low')));return(tmp$obs[1]>tmp$exp[1])})

LIHC.RBP.sig <- names(LIHC.RBP.pvalue)[LIHC.RBP.pvalue<0.05]

LIHC.RBP.anti <- names(LIHC.RBP.pvalue)[LIHC.RBP.pvalue<0.05&LIHC.RBP.risk]

LIHC.PI.RBP <- LIHC.PI.DS[LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL),]

LIHC.PI.RBP.chisq <- apply(LIHC.PI.RBP,1,function(x) survdiff(my.surv~ifelse(x>median(x),'high','low'))$chisq)

LIHC.PI.RBP.pvalue <- pchisq(LIHC.PI.RBP.chisq , 1, lower.tail = FALSE)

names(LIHC.PI.RBP.pvalue) <- LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL)]

LIHC.PI.RBP.sig <- unique(names(LIHC.PI.RBP.pvalue)[LIHC.PI.RBP.pvalue<0.05])

LIHC.PI.RBP2.sig <-LIHC.PI.RBP.sig[LIHC.PI.RBP.sig%in%LIHC.RBP.sig]

LIHC.RBP.risk <- sapply(rownames(LIHCMatrix),function(x){tmp<- survdiff(my.surv~t(ifelse(LIHCMatrix[x,]>median(as.numeric(LIHCMatrix[x,])),'high','low')));
                                                          return(tmp$obs[1]>tmp$exp[1])})
names(LIHC.RBP.pvalue)[LIHC.RBP.pvalue<0.05&LIHC.RBP.risk]
# cox

library(survival)

cox.data <- as.data.frame(t(log2(LIHCMatrix.RBP+1)))

cox.data$OS <- as.numeric(LIHC.clinical$OS)

cox.data$status <- LIHC.clinical$vital_status=='dead'

colnames(cox.data) <- gsub('-',"_",colnames(cox.data))

covariates <- colnames(cox.data)[1:1644]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = cox.data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(do.call(rbind,univ_results), stringsAsFactors=F)

res$HR <- as.numeric(gsub('(.+) \\(.+','\\1',res$`HR (95% CI for HR)`))

res$lower <- as.numeric(gsub('.+\\((.+)-.+','\\1',res$`HR (95% CI for HR)`))

res$higher <- as.numeric(gsub('.+-(.+)\\)','\\1',res$`HR (95% CI for HR)`))

res <- res[order(res$beta,decreasing = T)[1:10],]
res$Genes <- row.names(res)

rbind(colnames(res)[c(8,2,4)],res[,c(8,2,4)])

library(forestplot)

forestplot(rbind(colnames(res)[c(8,2,4)],res[,c(8,2,4)]),  #显示的文本
           c(NA,res$HR), #误差条的均值(此处为差值的中值)
           c(NA,res$lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,res$higher), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.35, ##误差条中的圆心点大小
           col=fpColors(line = "#ff6362", #误差条的线的颜色
                        box="#00b2e3"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 4,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 1.2), xlab = gpar(cex = 1), cex = 1.2), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Differences in assessment indicators between relevant pairs" #x轴的标题
)


graph2ppt(file='E:/LIHC_splicing/RBPs_forestplot.pptx',width=12,height=10)



#DEG

query <- GDCquery(project = "TCGA-LIHC",##对应癌症
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query, method = "api", files.per.chunk = 10)

LIHCcount  <- GDCprepare(query,summarizedExperiment = F)

LIHCcount$X1 <-sapply(strsplit(LIHCcount$X1,'\\.'),'[',1)

LIHCcount <- LIHCcount[grepl('ENSG',LIHCcount$X1),]


LIHCcount <- merge(LIHCcount,ensemble,by.x='X1',by.y='Gene.stable.ID')


LIHCcount <- aggregate(LIHCcount[,2:425],by=list(LIHCcount[,426]),mean)

row.names(LIHCcount) <- LIHCcount$Group.1

LIHCcount <- LIHCcount[,-1]

group <- !grepl('-11A',colnames(LIHCcount))

edgeRf <- function(data,gp){
  
  library("edgeR")
  
  data <- DGEList(counts=data,group=gp)
  
  keep <- rowSums(cpm(data)>1) >= 2
  
  data <- data[keep, keep.lib.sizes=FALSE]
  
  data <- calcNormFactors(data)
  
  fit <- estimateDisp(data)
  
  et <- exactTest(fit)
  
  DEG <- et$table
  
  DEG$P.ajust <- p.adjust(DEG$PValue,method = 'BH')
  
  return(DEG)
}

LIHC.DEGs <- edgeRf(LIHCcount,group)

LIHC.DEGs.RBP <- LIHC.DEGs[row.names(LIHC.DEGs)%in%c(RBP$SYMBOL,RDBP)&LIHC.DEGs$P.ajust<0.05,]

LIHC.DEGs <- LIHC.DEGs[abs(LIHC.DEGs$logFC)>1&LIHC.DEGs$P.ajust<0.05,]

LIHC.DEGs.PI <- LIHC.DEGs[row.names(LIHC.DEGs)%in%LIHC.PI.DS.data$symbol,]

#table(LIHC.PI.DS.data$splice_type[LIHC.PI.DS.data$symbol%in%row.names(LIHC.DEGs)[LIHC.DEGs$logFC<0]])

#table(LIHC.PI.DS.data$splice_type[LIHC.PI.DS.data$symbol%in%row.names(LIHC.DEGs)[LIHC.DEGs$logFC>0]])

LIHC.DEGs2 <- edgeRf(LIHCcount,group)

LIHC.DEGs2 <- LIHC.DEGs2[LIHC.DEGs2$P.ajust<0.05,]

library(clusterProfiler)

up.GO <- enrichGO(gene=row.names(LIHC.DEGs2)[LIHC.DEGs2$logFC>0],OrgDb="org.Hs.eg.db",keyType="SYMBOL",ont ='BP',minGSSize = 30)

a<-data.frame(up.GO)

down.GO <- enrichGO(gene=row.names(LIHC.DEGs)[LIHC.DEGs$logFC<0],OrgDb="org.Hs.eg.db",keyType="SYMBOL",ont ='BP',minGSSize = 30)

down.GO <- as.data.frame(down.GO)

write.table(unique(LIHC.PI.DS.data$symbol),'E:/LIHC_splicing/DEASs.txt',quote = F,row.names = F)

PI.GO <-  GOPanel(unique(LIHC.PI.DS.data$symbol))

library(clusterProfiler)

PI.GO.enrichmap <- enrichGO(gene=unique(LIHC.PI.DS.data$symbol),OrgDb="org.Hs.eg.db",keyType="SYMBOL",ont ='BP',minGSSize = 30)

PI.GO.enrichmap <- as.data.frame(PI.GO.enrichmap)

PI.GO.enrichmap <- PI.GO.enrichmap[,c(1:2,5:8)]

PI.GO.enrichmap$geneID <- gsub('/',',',PI.GO.enrichmap$geneID)

PI.GO.types <- sapply(unique(LIHC.PI.DS.data$splice_type),function(x){enrichGO(
                      gene=unique(LIHC.PI.DS.data$symbol[LIHC.PI.DS.data$splice_type==x]),OrgDb="org.Hs.eg.db",keyType="SYMBOL",ont ='BP',minGSSize = 30)})

PI.GO.types <- sapply(PI.GO.types,function(x){tmp <- as.data.frame(x);tmp <- tmp[,c(1:2,5:8)];tmp$geneID <- gsub('/',',',tmp$geneID);return(tmp)},simplify = F)

mapply(function(x,y) write.table(x,paste('E:/LIHC_splicing/PI_GO/',y,'GO.txt'),quote = F,row.names = F,sep='\t'),
       PI.GO.types,names(PI.GO.types))


write.table(PI.GO.enrichmap,'E:/LIHC_splicing/PI_GO/PI_GO_enrichmap.txt',quote = F,row.names = F,sep='\t')


library(VennDiagram)


## venn plot for overlaps between DEGs and DNMs

venn.diagram(list(row.names(LIHC.DEGs),LIHC.PI.DS.data$symbol),filename = "E:/LIHC_splicing/DEGs&DESs.tif",fill =c("#8484FB", "#CBC9FB"),cex=0,category = rep("", 2))

## venn plot for overlaps in DNMTs associated genes

#venn.diagram(list(overlap.mrna.methy,unique(Reduce(c,input))),
#             filename = "E:/colorectal cancer/WGCNA_output/vennDNMTs.tif",fill =c("#8484FB", "#CBC9FB"),cex=0,category = rep("", 2))


# splicing regulation

blue <- modules.GOs$blue$geneID[grep("splic",modules.GOs$blue$Description)]

blue.genes <- unique(unlist(strsplit(blue,'/')))

blue.RBPs <- row.names(LIHC.RBP.p.value)[LIHC.RBP.p.value[,'MEblue']<0.05]

blue.genes <- blue.genes[blue.genes%in%blue.RBPs]

magenta <- modules.GOs$magenta$geneID[grep("splic",modules.GOs$magenta$Description)]

magenta.genes <- unique(unlist(strsplit(magenta,'/')))

magenta.RBPs <- row.names(LIHC.RBP.p.value)[LIHC.RBP.p.value[,'MEmagenta']<0.05]

magenta.genes <- magenta.genes[magenta.genes%in%magenta.RBPs]

tan <- modules.GOs$tan$geneID[grep("splic",modules.GOs$tan$Description)]

tan.genes <- unique(unlist(strsplit(tan,'/')))

tan.RBPs <- row.names(LIHC.RBP.p.value)[LIHC.RBP.p.value[,'MEtan']<0.05]

tan.genes <- tan.genes[tan.genes%in%tan.RBPs]

combined <- c()

all <- PI.GO.enrichmap$geneID[grep("splic",PI.GO.enrichmap$Description)]

all.genes <- unique(unlist(strsplit(all,',')))

#all.RBPs <- names(LIHC.PI.RBP.pvalue)[LIHC.PI.RBP.pvalue<0.05]

all.RBPs <- row.names(LIHC.RBP.p.value)[LIHC.RBP.p.value<0.05]

all.genes <- all.genes[all.genes%in%all.RBPs]


library(VennDiagram)

library(gplots)

input <-list(blue=blue.genes ,magenta=magenta.genes,tan=tan.genes)

out.list <- attr(venn(input), "intersections") # intersections results





#### clusters

library(ConsensusClusterPlus)


cc <- ConsensusClusterPlus(LIHC.PI.DS,maxK=9,reps=1000,pItem=0.8,pFeature=1,title="1000pamLIHC",distance="euclidean",
                           clusterAlg="pam",plot="png")

icl <-  calcICL(cc,title="1000pamLIHC",plot="png")


icl[["clusterConsensus"]]

clusters <- cc[[4]]$consensusClass 

names(clusters) <- colnames(LIHC.PI.DS)



plot(fit,col=c('red','blue','green','yellow'))

cc.rbp <- ConsensusClusterPlus(as.matrix(LIHC.PI.DS[LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL),]),maxK=9,reps=1000,pItem=0.8,pFeature=1,title="1000pamLIHCRBP",distance="euclidean",
                           clusterAlg="pam",plot="png")

icl.rbp  <-  calcICL(cc.rbp ,title="1000pamLIHCRBP",plot="png")

icl.rbp[["clusterConsensus"]]

clusters.rbp <- cc.rbp[[8]]$consensusClass 

names(clusters.rbp) <- colnames(LIHC.PI.RBP)

clusters.rbp.df <- data.frame(OS=LIHC.clinical$OS,stage=LIHC.clinical$tumor_stage,Cluster=clusters.rbp,
                          status=LIHC.clinical$vital_status=='dead',RBP <-colMeans(LIHCMatrix),stringsAsFactors = F)

fit.rbp <- survfit(Surv(as.numeric(OS),status)~Cluster,data = clusters.rbp.df)

ggsurvplot(fit.rbp,data=clusters.rbp.df,pval = T)

#hub

cc.hub <- ConsensusClusterPlus(LIHC.PI.DS[abs(LIHC.GMM$MMS)>0.8&LIHC.GMM$module!='grey'&LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL),
                                          ],maxK=9,reps=1000,pItem=0.8,pFeature=1,title="1000pamLIHChub",distance="euclidean",
                               clusterAlg="pam",plot="pdf")

icl.hub  <-  calcICL(cc.hub ,title="1000pamLIHChub",plot="pdf")

icl.hub[["clusterConsensus"]]

clusters.hub <- cc.hub[[6]]$consensusClass 

names(clusters.hub) <- colnames(LIHC.PI.RBP)

clusters.hub.df <- data.frame(OS=LIHC.clinical$OS,stage=LIHC.clinical$tumor_stage,Cluster=clusters.hub,
                              status=LIHC.clinical$vital_status=='dead',RBP <-colMeans(LIHCMatrix),stringsAsFactors = F)

fit.hub <- survfit(Surv(as.numeric(OS),status)~Cluster,data = clusters.hub.df)

ggsurvplot(fit.hub,data=clusters.hub.df,pval = T)




library("survminer")

library(survival)


LIHC.clinical$tumor_stage <- gsub('[a-c]','',LIHC.clinical$tumor_stage)

LIHC.clinical$tumor_stage<- gsub('stge ','',LIHC.clinical$tumor_stage)

LIHC.clinical$tumor_stage <- as.numeric(as.roman(LIHC.clinical$tumor_stage))



clusters.df <- data.frame(OS=LIHC.clinical$OS,stage=LIHC.clinical$tumor_stage,Cluster=clusters,
                          status=LIHC.clinical$vital_status=='dead',stringsAsFactors = F)

fit <- survfit(Surv(as.numeric(OS)/30,status)~Cluster,data = clusters.df)

res.cox <- coxph(Surv(as.numeric(OS), status) ~Cluster,data =  clusters.df)

summary(res.cox)

fit <- survfit(Surv(time,status) ~ sex, data = lung)

ggsurvplot(fit,data=clusters.df,pval = T)

ggsurvplot(fit, data = lung)

survdiff(Surv(as.numeric(OS),stage)~lustes,data = clusters.df)

boxplot(results)

results <- sapply(unique(clusters.hub),function(x) colMeans(log2(LIHCMatrix.RBP+1)[,names(clusters.hub)[clusters.hub==x]]),simplify = F)

boxplot(results)


p.matrix <- sapply(1:5,function(group) apply(LIHCMatrix.RBP,1,function(x) t.test(x[clusters.hub==group],x[clusters.hub==6],alternative = 'less')$p.value))


p.adjust.matrix <- apply(p.matrix,2,function(x) p.adjust(x,method='BH'))

RBP.subtype6 <- row.names(p.adjust.matrix)[apply(p.adjust.matrix <0.01,1,all)]

RBP.subtype6.heatmap <- t(apply(LIHCMatrix.RBP[RBP.subtype6,],1,scale))

RBP.subtype6.heatmap <-RBP.subtype6.heatmap[,order(clusters.hub,decreasing = T)]

library(ggplot2)

library(gplots)

heatmap.2(RBP.subtype6.heatmap,col=greenred(75),scale=NULL,margins = c(11,5),Colv = F,
          key=TRUE,keysize=0.5,symkey=FALSE, density.info="none", trace="none", cexRow=0.0001,cexCol =1.5)

hub.rbp.events <- LIHC.PI.DS[abs(LIHC.GMM$MMS)>0.8&LIHC.GMM$module!='grey'&LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL),]


rbp.events.cor <- cor(t(LIHCMatrix.RBP[RBP.subtype6,]),t(hub.rbp.events))

rbp.events.pvalue <- corPvalueStudent(rbp.events.cor,371)

rbp.events.pvalue <- apply(rbp.events.pvalue,1,function(x) p.adjust(x,method='BH'))

rowSums(t(rbp.events.pvalue<0.05))

names(results) <- c('Cluster1','Cluster2','Cluster3','Cluster4')


LIHC.RBP.sig.exp <-log2(LIHCMatrix.RBP[rownames(LIHCMatrix.RBP)%in%LIHC.RBP.sig,]+1)

LIHC.RBP.sig.onco <-log2(LIHCMatrix.RBP[rownames(LIHCMatrix.RBP)%in%LIHC.RBP.sig&
                                        rownames(LIHCMatrix.RBP)%in%row.names(LIHC.DEGs.RBP)[LIHC.DEGs.RBP$logFC>0],]+1)

results.stage. <- sapply(1:3,function(x){
  colMeans(LIHC.RBP.sig.notanti[,LIHC.clinical$submitter_id[!is.na(LIHC.clinical$tumor_stage)][LIHC.clinical$tumor_stage[!is.na(LIHC.clinical$tumor_stage)]==x]])},
  simplify = F)

results.stage.sig <- sapply(1:3,function(x){
  colMeans(LIHC.RBP.sig.exp[,LIHC.clinical$submitter_id[!is.na(LIHC.clinical$tumor_stage)][LIHC.clinical$tumor_stage[!is.na(LIHC.clinical$tumor_stage)]==x]])},
  simplify = F)

results.stage.onco <- sapply(1:3,function(x){
  colMeans(LIHC.RBP.sig.onco[,LIHC.clinical$submitter_id[!is.na(LIHC.clinical$tumor_stage)][LIHC.clinical$tumor_stage[!is.na(LIHC.clinical$tumor_stage)]==x]])},
  simplify = F)

results.stage.all <-c(results.stage,results.stage.sig,results.stage.onco)

legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       box.lty=0)

boxplot(results.stage.all)

results <- sapply(unique(clusters),function(x) colMeans(LIHC.RBP.sig.exp[,names(clusters)[clusters==x]]),simplify = F)

names(results) <- c('Cluster1','Cluster2','Cluster3','Cluster4')

boxplot(results)

LIHC.RBP.DEG <-log2(LIHCMatrix[rownames(LIHCMatrix)%in%row.names(LIHC.DEGs.RBP),]+1)

results <- sapply(unique(clusters),function(x) colMeans(LIHC.RBP.DEG[,names(clusters)[clusters==x]]),simplify = F)

names(results) <- c('Cluster1','Cluster2','Cluster3','Cluster4')

boxplot(results)

results <- sapply(unique(clusters),function(x) colMeans(LIHC.RBP.sig.onco[,names(clusters)[clusters==x]]),simplify = F)

names(results) <- c('Cluster1','Cluster2','Cluster3','Cluster4')

boxplot(results)

LIHC.RBP.auto <-log2(LIHCMatrix[rownames(LIHCMatrix)%in%LIHC.DEGs.RBP.PI,]+1)

results <- sapply(unique(clusters),function(x) colMeans(LIHC.RBP.auto[,names(clusters)[clusters==x]]),simplify = F)

names(results) <- c('Cluster1','Cluster2','Cluster3','Cluster4')

boxplot(results)



#mutation 
library(dplyr)
library(DT)


maf <- GDCquery_Maf("LIHC", pipelines = "muse")

maf$Tumor_Sample_Barcode <-gsub('(TCGA-.+-.+)-.+-.+-.+-.+','\\1',maf$Tumor_Sample_Barcode)


mutation_analysis <- function(data)  # function to analyze mutation data for Network comparison
{ 
  data <- data[data$Hugo_Symbol%in%c(RBP$SYMBOL,RDBP),c('Hugo_Symbol','Tumor_Sample_Barcode')];
  data<- sapply(unique(data$Tumor_Sample_Barcode),function(x){
    tmp <- data[data$Tumor_Sample_Barcode==x,];
    return(t(unique(tmp[,1])))})
  return(data)
}
LIHC.mutation <- mutation_analysis(maf)

LIHC.mutation.table <- sapply(unique(clusters.hub),function(x) as.data.frame(table(unlist(LIHC.mutation[names(LIHC.mutation.length)%in%names(clusters.hub)[clusters.hub==x]]))),simplify = F)

LIHC.mutation.subtype <- sapply(unique(clusters.hub),function(x) LIHC.mutation.length[names(LIHC.mutation.length)%in%names(clusters.hub)[clusters.hub==x]])

boxplot(LIHC.mutation.subtype)

#CNV

query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")

GDCdownload(query)

CNV.data <- GDCprepare(query)


LIHC.mutation.matrix <- matrix(data=0,nrow=sum(c(RDBP,RBP$SYMBOL)%in%maf$Hugo_Symbol),ncol=length(LIHC.mutation))

row.names(LIHC.mutation.matrix) <- c(RDBP,RBP$SYMBOL)[c(RDBP,RBP$SYMBOL)%in%maf$Hugo_Symbol]

colnames(LIHC.mutation.matrix) <- names(LIHC.mutation)

LIHC.mutation.matrix[t(sapply(row.names(LIHC.mutation.matrix),function(y) sapply(LIHC.mutation,function(x) y%in%x)))]<-1

LIHC.mutation.matrix <-LIHC.mutation.matrix[,colnames(LIHC.mutation.matrix)%in%names(clusters)[clusters==1]]

events <- discover.matrix(LIHC.mutation.matrix)

subset <- rowSums(LIHC.mutation.matrix) > 5

result.mutex <- pairwise.discover.test(events[subset,])

print(result.mutex, fdr.threshold=0.3)

as.data.frame(result.mutex)

# hub genes in the network #MMS 

geneModuleMembership <- as.data.frame(cor(t(LIHC.PI.DS),LIHC.MEs, use = "p"))

colnames(geneModuleMembership) <- substring(colnames(geneModuleMembership),3)

LIHC.GMM <- mapply(function(module,MMS,modulecolor){out <- data.frame(MMS[modulecolor==module,module]);
                                    out$'module'<-module;row.names(out)<- row.names(MMS)[modulecolor==module];
                                    colnames(out)[1]<-'MMS';return(out)},
                                     unique(mergedColors.LIHC),list(geneModuleMembership),
                                     list(mergedColors.LIHC),SIMPLIFY = F)

LIHC.GMM <- do.call(rbind,LIHC.GMM)

rownames(LIHC.GMM) <- as.numeric(sapply(strsplit(rownames(LIHC.GMM),'\\.'),'[',2))

LIHC.GMM <- LIHC.GMM[order(as.numeric(rownames(LIHC.GMM))),]

LIHC.GMM.list <- list(
  RBP=abs(LIHC.GMM$MMS[LIHC.PI.pvalue < 0.05&LIHC.GMM$module!='grey'&LIHC.PI.DS.data$symbol%in%c(RDBP,RBP$SYMBOL)]),
     sig=abs(LIHC.GMM$MMS[LIHC.PI.pvalue < 0.05&LIHC.GMM$module!='grey']),
     NS=abs(LIHC.GMM$MMS[LIHC.PI.pvalue>=0.05&LIHC.GMM$module!='grey']))

boxplot(LIHC.GMM.list)

# RBP self regulation


LIHC.DEGs.RBP.PI <- unique(rownames(LIHC.DEGs.RBP)[rownames(LIHC.DEGs.RBP)%in%LIHC.PI.DS.data$symbol])


RDBP.eCLIP.sample <- read.delim('E:/RBP/RBP_Clip.tsv',stringsAsFactors = F)

RDBP.eCLIP.sample <- RDBP.eCLIP.sample[RDBP.eCLIP.sample$Assembly=='GRCh38'&RDBP.eCLIP.sample$Biosample.term.name=='HepG2',]

RDBP.eCLIP.sample$Experiment.target <- gsub('(.+)-human','\\1',RDBP.eCLIP.sample$Experiment.target)

RBP.self <- lapply(unique(RDBP.eCLIP.sample$Experiment.target),
              function(x){tmp <- RDBP.eCLIP.sample$File.accession[RDBP.eCLIP.sample$Experiment.target%in%x];
              Tar <- unlist(RBPs.targets[tmp]);if (x%in%Tar) return(x)})

RBP.self <-unlist(RBP.self)

LIHC.DEGs.RBP.PI[LIHC.DEGs.RBP.PI%in%RBP.self]

modules.splicing.genes <- unique(c(tan.genes,blue.genes,magenta.genes))



log2(rowMeans(LIHCMatrix[row.names(LIHCMatrix)%in%LIHC.DEGs.RBP.PI[LIHC.DEGs.RBP.PI%in%RBP.self],])+1),
log2(rowMeans(LIHCMatrix[!row.names(LIHCMatrix)%in%LIHC.DEGs.RBP.PI[LIHC.DEGs.RBP.PI%in%RBP.self],])+1)


boxplot(log2(rowMeans(LIHCMatrix[row.names(LIHCMatrix)%in%LIHC.DEGs.RBP.PI,])+1),log2(rowMeans(LIHCMatrix.RBP[!row.names(LIHCMatrix.RBP)%in%LIHC.DEGs.RBP.PI,])+1),
        log2(rowMeans(LIHCMatrix[!row.names(LIHCMatrix)%in%rownames(LIHCMatrix.RBP),])+1))

IU <- LIHC.PI.DS.data[LIHC.PI.DS.data$symbol=='ILF3'|LIHC.PI.DS.data$symbol=='UPF1',,drop = F]


write.csv(t(IU[,grepl('Norm',colnames(IU))]),file = 'E:/LIHC_splicing/IU_normal.csv',quote = F)

normal <- t(IU[,grepl('Norm',colnames(IU))])


tumor <- t(IU[,!grepl('Norm',colnames(IU))])[-c(1:10,382),]


write.csv(tumor,file = 'E:/LIHC_splicing/IU_tumor.csv',quote = F)

var.test(as.numeric(normal[,4]),as.numeric(tumor[,4]))

t.test(as.numeric(normal[,2]),as.numeric(tumor[,2]),var.equal = F)
