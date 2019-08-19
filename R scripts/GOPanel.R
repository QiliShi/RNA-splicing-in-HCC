GOPanel <- function(genename,Organism="org.Hs.eg.db",ko='hsa'){

  require(Organism,character.only = T)

  library(clusterProfiler)
  
  cat(paste('gene size is',length(genename),sep=' '))
  
  cat('\n')
  
  if (length(genename) >= 10 & length(genename) <= 10000 ){
    
    cat('GO analysing ...')

    go_bp <- enrichGO(gene=genename,OrgDb=Organism,keyType="SYMBOL",ont ='BP')
    
    go_mf <- enrichGO(gene=genename,OrgDb=Organism,keyType="SYMBOL",ont ='MF')
    
    go_cc <- enrichGO(gene=genename,OrgDb=Organism,keyType="SYMBOL",ont ='CC')
    
    eg <- bitr(genename, fromType="SYMBOL", toType="ENTREZID",Organism)

    kk<- enrichKEGG(gene=eg$ENTREZID,organism=ko,pvalueCutoff = 0.05)
    
    go_result<-do.call(rbind,list(data.frame(go_bp),data.frame(go_mf),data.frame(go_cc),data.frame(kk)))
  
    return(go_result)
    
  }
    
  else  return('Gene size is too small or too large')
}