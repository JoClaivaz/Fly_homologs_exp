

directory <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/"
setwd(directory)


# Load library and data ---------------------------------------------------

library(vioplot)

load("Data/Human/GSE30611.RData")
GSE30611 <- tissues.list
load("Data/Human/GSE30352.RData")
GSE30352 <- tissues.list
load("Data/Human/GSE43520.RData")
GSE43520 <- tissues.list
rm(tissues.list)

gse_list <- list(GSE30611, GSE30352, GSE43520)


# mouse data
load("Data/Mouse/GSE30352.RData")
GSE30352_mouse <- tissues.list
load("Data/Mouse/GSE30617.RData")
GSE30617_mouse <- tissues.list
load("Data/Mouse/GSE36026.RData")
GSE36026_mouse <- tissues.list
load("Data/Mouse/GSE41338.RData")
GSE41338_mouse <- tissues.list
load("Data/Mouse/GSE41637.RData")
GSE41637_mouse <- tissues.list
load("Data/Mouse/GSE43520.RData")
GSE43520_mouse <- tissues.list
load("Data/Mouse/GSE43721.RData")
GSE43721_mouse <- tissues.list
rm(tissues.list)

gse_list_mouse <- list(GSE30352_mouse, GSE30617_mouse, GSE36026_mouse, GSE41338_mouse, GSE41637_mouse, GSE43520_mouse, GSE43721_mouse)


# Kruskal Wallis test -----------------------------------------------------

## make a funtion for testing  (Kruskal or correlation test)

# Kruskal Wallis test for testing variation of gene expression in function of the domain difference:
# 

kruskal_result <- NULL
for (gse in gse_list) {
  test <- NULL
  for (tissue in names(gse)){
    p <- as.matrix(kruskal.test(log(gse[[tissue]]$RPKM_short_norm) ~ domain_difference, data = gse[[tissue]])$p.value)
    rownames(p) <- tissue
    test <- rbind(test, p)
  }
  kruskal_result <- rbind(kruskal_result, test)
  colnames(kruskal_result) <- "p_value"
}
write.table(round(kruskal_result,2), paste0(directory, "Data/Output/kruskal_result_human.txt"), quote = F, col.names = F, sep = "\t")


# Spearmann correlation test: correlation between paralogs and orthologs expression
# 

corr_result <- NULL
for (gse in gse_list) {
  test <- NULL
  for (tissue in names(gse)){
    rho <- as.matrix(cor.test(log(gse[[tissue]]$RPKM_long), log(gse[[tissue]]$RPKM_short), method = "spearman")$estimate)
    rownames(rho) <- tissue
    test <- rbind(test, rho)
  }
  corr_result <- rbind(corr_result, test)
  colnames(corr_result) <- "p_value"
}
write.table(round(corr_result,2), paste0(directory, "Data/Output/rpkm_corr_result_human.txt"), quote = F, col.names = F, sep = "\t")


## mouse
# Kruskal Wallis test for testing variation of gene expression in function of the domain difference:
# 

kruskal_result <- NULL
for (gse in gse_list_mouse) {
  test <- NULL
  for (tissue in names(gse)){
    p <- as.matrix(kruskal.test(log(gse[[tissue]]$RPKM_short_norm) ~ domain_difference, data = gse[[tissue]])$p.value)
    rownames(p) <- tissue
    test <- rbind(test, p)
  }
  kruskal_result <- rbind(kruskal_result, test)
  colnames(kruskal_result) <- "p_value"
}
write.table(round(kruskal_result,2), paste0(directory, "Data/Output/kruskal_result_mouse.txt"), quote = F, col.names = F, sep = "\t")


# Spearmann correlation test: correlation between paralogs and orthologs expression
# 

corr_result <- NULL
for (gse in gse_list_mouse) {
  test <- NULL
  for (tissue in names(gse)){
    rho <- as.matrix(cor.test(log(gse[[tissue]]$RPKM_long), log(gse[[tissue]]$RPKM_short), method = "spearman")$estimate)
    rownames(rho) <- tissue
    test <- rbind(test, rho)
  }
  corr_result <- rbind(corr_result, test)
  colnames(corr_result) <- "p_value"
}
write.table(round(corr_result,2), paste0(directory, "Data/Output/rpkm_corr_result_mouse.txt"), quote = F, col.names = F, sep = "\t")

## comparing human-mouse orthologs
# Spearman correlation test: correlation between human and mouse orthologs expression
# 

orthologs <- read.table("Data/Output/mouse_human_orthologs2017-02-22.txt", h = T)

test <- NULL
for (tissue in names(GSE30352_mouse)[names(GSE30352_mouse) %in% names(GSE30352)]){
  gse_mouse <- unique(rbind(GSE30352_mouse[[tissue]][GSE30352_mouse[[tissue]]$short_domain %in% orthologs$Mouse_Ensembl_peptide_ID_short,],
                            GSE30352_mouse[[tissue]][GSE30352_mouse[[tissue]]$long_domain %in% orthologs$Mouse_Ensembl_peptide_ID_long,]))
  gse_human <- unique(rbind(GSE30352[[tissue]][GSE30352[[tissue]]$short_domain %in% orthologs$short_domain,],
                            GSE30352[[tissue]][GSE30352[[tissue]]$long_domain %in% orthologs$long_domain,]))
  
  es_short <- merge(orthologs, GSE30352_mouse[[tissue]][c("Ensembl_Gene_ID_short", "RPKM_short")], by.x = "Mouse_Ensembl_Gene_ID_short",
                    by.y = "Ensembl_Gene_ID_short")
  rho <- as.matrix(cor.test(log(GSE30352_mouse[[tissue]]$RPKM_long), log(GSE30352[[tissue]]$RPKM_long), method = "spearman")$estimate)
  rownames(rho) <- tissue
  test <- rbind(test, rho)
}
gse_mouse <- unique(rbind(GSE30352_mouse[[1]][GSE30352_mouse[[1]]$short_domain %in% orthologs$Mouse_Ensembl_peptide_ID_short,],
                   GSE30352_mouse[[1]][GSE30352_mouse[[1]]$long_domain %in% orthologs$Mouse_Ensembl_peptide_ID_long,]))

test <- NULL
for (tissue in names(GSE30352_mouse)[names(GSE30352_mouse) %in% names(GSE30352)]){
  rho <- as.matrix(cor.test(log(GSE30352_mouse[[tissue]]$RPKM_long), log(GSE30352[[tissue]]$RPKM_long), method = "spearman")$estimate)
  rownames(rho) <- tissue
  test <- rbind(test, rho)
}

# write.table(corr_result, paste0(directory, "Output/corr_result_mouse.txt"), quote = F, col.names = F)


smoothScatter(log(GSE30352$`frontal cortex`$RPKM_long), log(GSE30352$`frontal cortex`$RPKM_short), pch = 20, cex = 0.3, nrpoints = Inf,
              main = "Human protein expression in .. ", xlab = "RPKM of protein with longer protein domain", ylab = "RPKM of protein with shorter domain")

boxplot(log(GSE30352$`frontal cortex`$RPKM_long[is.finite(log(GSE30352$`frontal cortex`$RPKM_long))]), 
        log(GSE30352$`frontal cortex`$RPKM_short)[is.finite(log(GSE30352$`frontal cortex`$RPKM_short))], col = c(2,7), 
        main = "Protein domain expression (RPKM) in ..", names = c("Long protein domain", "Short protein domain"))

boxplot(log(GSE30352$`frontal cortex`$RPKM_short_norm)[is.finite(log(GSE30352$`frontal cortex`$RPKM_short_norm))]~GSE30352$`frontal cortex`$domain_difference[is.finite(log(GSE30352$`frontal cortex`$RPKM_short_norm))], col = c(2,2,7,4,5,6,7,4,6,3),
        main = "protein domain expression")


vioplot(log(GSE30352$`frontal cortex`$RPKM_short_norm)[is.finite(log(GSE30352$`frontal cortex`$RPKM_short_norm))])

vioplot(log(GSE30352$`frontal cortex`$RPKM_short_norm)[is.finite(log(GSE30352$`frontal cortex`$RPKM_short_norm))]~GSE30352$`frontal cortex`$domain_difference[is.finite(log(GSE30352$`frontal cortex`$RPKM_short_norm))], col = c(2,2,7,4,5,6,7,4,6,3),
        main = "protein domain expression")


# calculate median expression for every short domain, long domain  --------

domain_med <- data.frame(matrix(NA, 10, 0))
rownames(domain_med) <- c(as.character(unique(GSE30611[[tissue]]$domain_difference)))
for (tissue in names(GSE30611)){
  short <- aggregate(log(RPKM_short_norm)[is.finite(log(RPKM_short_norm))] ~ domain_difference[is.finite(log(RPKM_short_norm))], data = GSE30611[[tissue]], median)
  # long  <- median(log(GSE30611[[tissue]]$RPKM_long)[is.finite(log(GSE30611[[tissue]]$RPKM_short_norm))], na.rm = T)
  rownames(short) <- short$`domain_difference[is.finite(log(RPKM_short_norm))]`
  short$`domain_difference[is.finite(log(RPKM_short_norm))]` <- NULL
  # long_short <- rbind(long, short)
  # rownames(long_short)[1] <- "long"
  domain_med[tissue] <- short[match(rownames(domain_med), rownames(short)),]
}

par(mfrow = c(1,1))
plot(domain_med$`adipose tissue`, type = "n", ylim = c(min(domain_med), max(domain_med)), xlab = "", 
     ylab = "RPKM value (log)", xaxt = "n", main = "Median of the protein domain expression for different tissues")
for (i in 1:length(domain_med)){
  points(domain_med[,i], col = rainbow(17)[i], pch = 19)
}
text(x = seq(1, 10, by=1), par("usr")[3] - 0.1, labels = rownames(domain_med), srt = 45, pos = 1, xpd = TRUE)
legend(6.5,0.5, colnames(domain_med), col = c(rainbow(17)), pch = 19, ncol = 2, 
       text.font = 6, bty = "n")


par(mfrow=c(4,5))
for(i in 1:length(domain_med)){
  barplot(domain_med[,i], ylim = c(0, 2.5), col = c(rainbow(17)[i]), main = colnames(domain_med)[i],
          names.arg = rownames(domain_med),las = 2)
}


par(mfrow=c(4,5))
for(i in 1:length(domain_med)){
  barplot(domain_med[,i], ylim = c(0, 2.5), col = c(rainbow(17)[i]), main = colnames(domain_med)[i],
          names.arg = rownames(domain_med),las = 2)
}

par(mfrow = c(4,5))
for(i in 1:ncol(domain_med)){
  vioplot(domain_med[,i], col = c(rainbow(17)[i]), names = colnames(domain_med)[i] )
}

# Visualize the data ------------------------------------------------------
path <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/Output/figure/tissues/"

for (tissue in names(GSE30611)){
  rpkm <- GSE30611[[tissue]][is.finite(log(GSE30611[[tissue]]$RPKM_short_norm)),]
    ggplot(rpkm, aes(factor(domain_difference), log(RPKM_short_norm))) + geom_violin(aes(fill = domain_difference)) + 
    guides(fill=FALSE) + theme(axis.text=element_text(size=30), 
                               axis.title=element_text(size=30,face="bold"), 
                               plot.title = element_text(size=50)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle(paste0("Violin plot of normalized expression for different short protein \ndomain in human's ", tissue))
ggsave(file = paste0(path, tissue, "_vioplot_domain_expression", Sys.Date(), ".png"), width = 20, height = 20)
}

ggplot(rpkm, aes(factor(domain_difference), log(RPKM_short_norm))) + geom_violin(aes(fill = domain_difference)) + 
  ggtitle("violin plot of normalized expression for different short protein domain in human") + theme_bw()

par(mfrow=c(4,5))
for(i in 1:length(domain_med)){
  barplot(domain_med[,i], ylim = c(0, 2.5), col = c(rainbow(17)[i]), main = colnames(domain_med)[i],
          names.arg = rownames(domain_med),las = 2)
}


## better correlation between paralogs in human: testis
pdf()
smoothScatter(log(GSE43520[["testis"]]$RPKM_long) ~ log(GSE43520[["testis"]]$RPKM_short), pch = 20, cex = 0.3, nrpoints = Inf, 
              main = "Human paralogs expression in testis", xlab = "ln(RPKM) of longer paralogs", ylab = "ln(RPKM) shorter paralogs")
text(x = -5, y = 8.5, "rho = 0.39", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE43520[["testis"]]$RPKM_long), (GSE43520[["testis"]]$RPKM_short), method = "spearman")

## worst correlation: lymp node
smoothScatter(log(GSE30611[["lymph node"]]$RPKM_long) ~ log(GSE30611[["lymph node"]]$RPKM_short), pch = 20, cex = 0.3, nrpoints = Inf, 
              main = "Human paralogs expression in lymph node", xlab = "ln(RPKM) of longer paralogs", ylab = "ln(RPKM) shorter paralogs")
text(x = -6, y = 10, "rho = 0.17", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE30611[["lymph node"]]$RPKM_long), (GSE30611[["lymph node"]]$RPKM_short), method = "spearman")

## mouse testis
pdf()
smoothScatter(log(GSE43520_mouse[["testis"]]$RPKM_long) ~ log(GSE43520_mouse[["testis"]]$RPKM_short), 
              pch = 20, cex = 0.3, nrpoints = Inf, main = "Mouse paralogs expression in testis", 
              xlab = "ln(RPKM) of longer protein", ylab = "ln(RPKM) shorter protein")
text(x = -4.7, y = 8.5, "rho = 0.36", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE43520_mouse[["testis"]]$RPKM_long), (GSE43520_mouse[["testis"]]$RPKM_short), method = "spearman")

## worst correlation: colon
smoothScatter(log(GSE41637_mouse[["colon"]]$RPKM_long) ~ log(GSE41637_mouse[["colon"]]$RPKM_short), pch = 20, cex = 0.3, nrpoints = Inf, 
              main = "Mouse paralogs expression in colon", xlab = "ln(RPKM) of longer paralogs", ylab = "ln(RPKM) shorter paralogs")
text(x = -6, y = 9, "rho = 0.21", xpd=TRUE, pos = 1, cex = 1.2)
cor.test((GSE41637_mouse[["colon"]]$RPKM_long), (GSE41637_mouse[["colon"]]$RPKM_short), method = "spearman")

