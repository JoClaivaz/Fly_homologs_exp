
directory <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/"
setwd(directory)

library(biomaRt)
library(BgeeDB)

## for orthologs: combine mouse and human databases!!!

load(file = "Data/base_data2017-02-22.RData") # biomart, bgee and protein domain informations
## common human - mouse experiment = GSE30352(3 for human, 5 for mouse!), GSE43520 (2 in human, 7 in mouse!!)

orthologs <- read.table("Data/Output/mouse_human_orthologs2017-02-22.txt", h = T)

## Combine databases

gse <- list()
for (tiss in unique(data_bgee_mouse[[7]]$Anatomical.entity.name)[unique(data_bgee_mouse[[7]]$Anatomical.entity.name) %in% unique(data_bgee_human[[2]]$Anatomical.entity.name)]){
  tiss_H <- data_bgee_human[[2]][data_bgee_human[[2]]$`Anatomical.entity.name` == tiss,][,c("Gene.ID", "RPKM")]
  tiss_M <- data_bgee_mouse[[7]][data_bgee_mouse[[7]]$`Anatomical.entity.name` == tiss,][,c("Gene.ID", "RPKM")]
  # for human
  rpkm_tiss_h <- NULL
  for(gene in unique(tiss_H$Gene.ID)){
    rpkm_mean_h <- mean(tiss_H$RPKM[tiss_H$Gene.ID==gene])
    rpkm_tiss_h <- c(rpkm_tiss_h, rpkm_mean_h)
  }
  tiss.h <- data.frame(Gene.ID = unique(tiss_H$Gene.ID), RPKM = rpkm_tiss_h)
  
  # for mice
  rpkm_tiss_m <- NULL
  for(gene in unique(tiss_M$Gene.ID)){
    rpkm_mean_m <- mean(tiss_M$RPKM[tiss_M$Gene.ID==gene])
    rpkm_tiss_m <- c(rpkm_tiss_m, rpkm_mean_m)
  }
  tiss.m <- data.frame(Gene.ID = unique(tiss_M$Gene.ID), RPKM = rpkm_tiss_m)
  
  ## info for long protein domain
  tiss_long <- merge(orthologs, tiss.h, by.x = "Human_Ensembl_Gene_ID_long", by.y = "Gene.ID")
  colnames(tiss_long)[colnames(tiss_long) =="RPKM"] <- "Human_RPKM_long"
  tiss_long <- merge(tiss_long, tiss.m, by.x = "Mouse_Ensembl_Gene_ID_long", by.y = "Gene.ID")
  colnames(tiss_long)[colnames(tiss_long) =="RPKM"] <- "Mouse_RPKM_long"
  
  ## info for short protein domain
  tiss_short <- merge(orthologs, tiss.h, by.x = "Human_Ensembl_Gene_ID_short", by.y = "Gene.ID")
  colnames(tiss_short)[colnames(tiss_short) =="RPKM"] <- "Human_RPKM_short"
  tiss_short <- merge(tiss_short, tiss.m, by.x = "Mouse_Ensembl_Gene_ID_short", by.y = "Gene.ID")
  colnames(tiss_short)[colnames(tiss_short) =="RPKM"] <- "Mouse_RPKM_short"
  
  tissue <- merge(tiss_long, tiss_short, all = T)
  
  ## for human
  rpkm_par_h <- NULL
  for (gene_long in unique(tissue$Human_Ensembl_Gene_ID_long)){
    paral.rpkm.norm <- tissue$Human_RPKM_short[tissue$Human_Ensembl_Gene_ID_long==gene_long]/tissue$Human_RPKM_long[tissue$Human_Ensembl_Gene_ID_long==gene_long]
    rpkm_par_h <- c(rpkm_par_h, paral.rpkm.norm)
  }
  tissue$Human_RPKM_short_norm <- rpkm_par_h
  
  ## for mouse
  rpkm_par_m <- NULL
  for (gene_long in unique(tissue$Mouse_Ensembl_Gene_ID_long)){
    paral.rpkm.norm <- tissue$Mouse_RPKM_short[tissue$Mouse_Ensembl_Gene_ID_long==gene_long]/tissue$Mouse_RPKM_long[tissue$Mouse_Ensembl_Gene_ID_long==gene_long]
    rpkm_par_m <- c(rpkm_par_m, paral.rpkm.norm)
  }
  tissue$Mouse_RPKM_short_norm <- rpkm_par_m
  
  gse[[tiss]] <- tissue
}


# GSE30352_hm <- gse
GSE43520_hm <- gse

##  calculate the correlation between mouse and human orthologs
# load("Data/orthologs_data.RData")
# gse <- GSE43520_hm

corr_test <- NULL
for(prot.length in c("long", "short", "short_norm")){
  test <- NULL
  for (tissue in names(gse)){
    rho <- as.matrix(cor.test(log(gse[[tissue]][paste0("Human_RPKM_", prot.length)][,1]), 
                              log(gse[[tissue]][paste0("Mouse_RPKM_", prot.length)][,1]), method = "spearman")$estimate)
    rownames(rho) <- tissue
    test <- rbind(test, rho)
  }
corr_test <- cbind(corr_test, test)
}
colnames(corr_test) <- c("long", "short", "short_norm")

  
write.table(corr_test, paste0(directory, "Data/Output/corr_result_orthologs20161223.txt"), quote = F, sep = "\t", col.names = F)

## plots for presentation
## better correlation between orthologs: short - brain
smoothScatter(log(GSE43520_hm[["brain"]]$Human_RPKM_short) ~ log(GSE43520_hm[["brain"]]$Mouse_RPKM_short), pch = 20, cex = 0.3, nrpoints = Inf, 
              main = "Short orhologs expression in brain", xlab = "ln(RPKM) of human orthologs", 
              ylab = "ln(RPKM) of mouse orthologs")
text(x = -5, y = 7.5, "rho = 0.87", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE43520_hm[["brain"]]$Human_RPKM_short), (GSE43520_hm[["brain"]]$Mouse_RPKM_short), method = "spearman")

## long brain
smoothScatter(log(GSE43520_hm[["brain"]]$Human_RPKM_long) ~ log(GSE43520_hm[["brain"]]$Mouse_RPKM_long), pch = 20, cex = 0.3, nrpoints = Inf, 
              main = "Long orthologs expression in brain", xlab = "ln(RPKM) of human orthologs", 
              ylab = "ln(RPKM) of mouse orthologs")
text(x = -5, y = 8.5, "rho = 0.61", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE43520_hm[["brain"]]$Human_RPKM_long), (GSE43520_hm[["brain"]]$Mouse_RPKM_long), method = "spearman")


  