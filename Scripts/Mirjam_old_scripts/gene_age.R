
directory <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/"
setwd(directory)


library(ggplot2)

## taxonomic node data:
taxa_df <- read.table("Data/taxonAgeTable.txt", h=T, sep = "\t")
## tissue specificity score
score_human <- read.table("Data/Human/tissue_spec_score/GSE30611spec_score2017-02-22.txt") # GSE30611 experience with more tissues for human
colnames(score_human) <- c("gene.id", "score")
score_mouse <- read.table("Data/Mouse/tissue_spec_score/GSE36026spec_score2017-02-22.txt") # GSE36026 expreience with more tissues for mouse
colnames(score_mouse) <- c("gene.id", "score")

## gene - taxa data
gene_age <- read.table("Data/geneAge.txt", h=F, sep="\t")
colnames(gene_age) <- c("gene_id", "Taxon")
gene_age <- merge(gene_age, taxa_df)

# gene_age_human <- gene_age[grep("ENSG0", gene_age$gene_id),]
# gene_age_human <- merge(gene_age_human, taxa_df)
# gene_age_mouse <- gene_age[grep("ENSMUSG", gene_age$gene_id),]
# gene_age_mouse <- merge(gene_age_mouse, taxa_df)

## human-mouse paralogs
hu_mu_paralogs <- read.table("Data/Output/mouse_human_orthologs2017-02-22.txt", h= T)

es <- merge(hu_mu_paralogs, gene_age, by.x = "Human_Ensembl_Gene_ID_short", by.y = "gene_id")
colnames(es)[colnames(es) == "Taxon"] <- "Human_taxon_short"
colnames(es)[colnames(es) == "Age"] <- "Human_age_short"

es1 <- merge(es, gene_age, by.x = "Mouse_Ensembl_Gene_ID_short", by.y = "gene_id")
colnames(es1)[colnames(es1) == "Taxon"] <- "Mouse_taxon_short"
colnames(es1)[colnames(es1) == "Age"] <- "Mouse_age_short"

es2 <- merge(es1, gene_age, by.x = "Human_Ensembl_Gene_ID_long", by.y = "gene_id")
colnames(es2)[colnames(es2) == "Taxon"] <- "Human_taxon_long"
colnames(es2)[colnames(es2) == "Age"] <- "Human_age_long"

es3 <- merge(es2, gene_age, by.x = "Mouse_Ensembl_Gene_ID_long", by.y = "gene_id")
colnames(es3)[colnames(es3) == "Taxon"] <- "Mouse_taxon_long"
colnames(es3)[colnames(es3) == "Age"] <- "Mouse_age_long"

es4 <- merge(es3, score_human, by.x = "Human_Ensembl_Gene_ID_short", by.y = "gene.id")
colnames(es4)[colnames(es4) == "score"] <- "Human_score_short"

es5 <- merge(es4, score_mouse, by.x = "Mouse_Ensembl_Gene_ID_short", by.y = "gene.id")
colnames(es5)[colnames(es5) == "score"] <- "Mouse_score_short"

es6 <- merge(es5, score_human, by.x = "Human_Ensembl_Gene_ID_long", by.y = "gene.id")
colnames(es6)[colnames(es6) == "score"] <- "Human_score_long"

es7 <- unique(merge(es6, score_mouse, by.x = "Mouse_Ensembl_Gene_ID_long", by.y = "gene.id"))
colnames(es7)[colnames(es7) == "score"] <- "Mouse_score_long"

# es8 <- es7[es7$Human_taxon_long == es7$Human_taxon_short,]

# pearson correlation for tissue specificity score function ---------------

corrPearson <- function(taxa.list, data, taxa.info, 
                        type = "Correlation between paralogs gene expressing protein with different length (length) or orthologs", specie, specie2 = NULL, protein.length ){
  corr_pearson <- NULL
  for (taxa in unique(taxa.list)){
    if (type == "length") {
      data <- data[data[paste(specie, "taxon_long", sep = "_")] == data[paste(specie, "taxon_short", sep = "_")],]
      gene.taxa  <- data[data[paste(specie, "taxon", protein.length, sep = "_")] == taxa ,][c(paste(specie, "Ensembl_Gene_ID_long", sep = "_"), 
                                                                                              paste(specie, "Ensembl_Gene_ID_short", sep = "_"),
                                                                                              paste(specie, "score_long", sep = "_"), 
                                                                                              paste(specie, "score_short", sep = "_"))]
      gene_sub <- NULL
      for (gene in unique(gene.taxa[paste(specie, "Ensembl_Gene_ID_long", sep = "_")])[,1]){
        id <- gene.taxa[gene.taxa[paste(specie, "Ensembl_Gene_ID_long", sep = "_")][,1] == gene,
                        ] [sample(nrow(gene.taxa[gene.taxa[paste(specie, "Ensembl_Gene_ID_long", sep = "_")][,1] == gene,]),1),]
        gene_sub <- rbind(gene_sub, id)
      }
      if(nrow(gene_sub) <= 3){
        cor <- NA
        corr <- NA
      } else {
        cor <- cor.test(gene_sub[paste(specie, "score_short", sep = "_")][,1], gene_sub[paste(specie, "score_long", sep = "_")][,1])$estimate
        corr <- cbind(cor, nrow(gene_sub))
        rownames(corr) <- taxa
        corr_pearson <- rbind(corr_pearson, corr)
      }
    }
    
    if (type == "paralogs") {
      data <- data[data[paste(specie, "taxon", protein.length, sep = "_")] == data[paste(specie2, "taxon", protein.length, sep = "_")],]
      gene.taxa  <- data[data[paste(specie, "taxon", protein.length, sep = "_")] == taxa ,][c(paste(specie, "Ensembl_Gene_ID", protein.length, sep = "_"), 
                                                                                              paste(specie2, "Ensembl_Gene_ID", protein.length, sep = "_"),
                                                                                              paste(specie, "score", protein.length, sep = "_"), 
                                                                                              paste(specie2, "score", protein.length, sep = "_"))]
      gene_sub <- NULL
      for (gene in unique(gene.taxa[paste(specie, "Ensembl_Gene_ID", protein.length, sep = "_")])[,1]){
        id <- gene.taxa[gene.taxa[paste(specie, "Ensembl_Gene_ID", protein.length, sep = "_")][,1] == gene,
                        ] [sample(nrow(gene.taxa[gene.taxa[paste(specie, "Ensembl_Gene_ID", protein.length, sep = "_")] == gene,]),1),]
        gene_sub <- rbind(gene_sub, id)
      }
      if(nrow(gene_sub) <= 3){
        cor <- NA
        corr <- NA
      } else {
        cor <- cor.test(gene_sub[paste(specie, "score", protein.length, sep = "_")][,1], gene_sub[paste(specie2, "score", protein.length, sep = "_")][,1])$estimate
        corr <- cbind(cor, nrow(gene_sub))
        rownames(corr) <- taxa
        corr_pearson <- rbind(corr_pearson, corr)
      }
    }
  }
  corr.test <- data.frame(corr_pearson)
  colnames(corr.test)[2] <- "n"
  corr.test$age <- taxa_df$Age[taxa.info$Taxon %in% rownames(corr.test)]
  return(corr.test)
}



# calculate correlation ---------------------------------------------------

## for human 
human_corr <- corrPearson(taxa.list = es7$Human_taxon_long, taxa.info = taxa_df, data = es7, type = "length", specie = "Human", protein.length = "long")

## for mouse 
mouse_corr <- corrPearson(taxa.list = es7$Mouse_taxon_long, taxa.info = taxa_df, data = es7, type = "length", specie = "Mouse", protein.length = "long")


## orthologs
long_corr <- corrPearson(taxa.list = es7$Human_taxon_long, taxa.info = taxa_df, data = es7, type = "paralogs", specie = "Human", specie2 = "Mouse", protein.length = "long")
short_corr <- corrPearson(taxa.list = es7$Human_taxon_short, taxa.info = taxa_df, data = es7, type = "paralogs", specie = "Human", specie2 = "Mouse", protein.length = "short")


# Plots -------------------------------------------------------------------

# long and short orthologs plot
short_corr$homologs <- rep("Short orthologs", nrow(short_corr))
long_corr$homologs <- rep("Long orthologs", nrow(long_corr))
corr_long_short <- rbind(long_corr, short_corr)

ggplot(corr_long_short, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for gene expressiong short or long protein domain") +
  scale_size("Number of\ngenes")

# Human paralogs and long orthologs plot
human_corr$homologs <- rep("Paralogs", nrow(human_corr))
long_corr$homologs <- rep("Long orthologs", nrow(long_corr))
corr_par_ort_human <- rbind(human_corr, long_corr)

ggplot(corr_par_ort_human, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for paralogs and orthologs genes in human") +
  scale_size("Number of\ngenes")

# Human paral vs short orthologs plot
human_corr$homologs <- rep("Paralogs", nrow(human_corr))
short_corr$homologs <- rep("Short orthologs", nrow(short_corr))
corr_par_short_ort_human <- rbind(human_corr, short_corr)

ggplot(corr_par_short_ort_human, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for human paralogs and short orthologs") +
  scale_size("Number of\ngenes")


# Mouse paralogs vs long orthologs plot
mouse_corr$homologs <- rep("Paralogs", nrow(human_corr))
long_corr$homologs <- rep("Long orthologs", nrow(long_corr))
corr_par_ort_mouse <- rbind(long_corr, mouse_corr)

ggplot(corr_par_ort_mouse, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for paralogs and orthologs genes in mouse") +
  scale_size("Number of\ngenes")

# Mouse paralogs vs short orthologs plot
mouse_corr$homologs <- rep("Paralogs", nrow(human_corr))
short_corr$homologs <- rep("Short orthologs", nrow(short_corr))
corr_par_short_ort_mouse <- rbind(mouse_corr, short_corr)

ggplot(corr_par_short_ort_mouse, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for paralogs and orthologs genes in mouse") +
  scale_size("Number of\ngenes")

# human paral vs mouse paralogs
human_corr$homologs <- rep("Human paralogs", nrow(human_corr))
mouse_corr$homologs <- rep("Mouse paralogs", nrow(mouse_corr))
corr_par_human_mouse <- rbind(human_corr, mouse_corr)

ggplot(corr_par_human_mouse, aes(age, cor, size=n, col=homologs)) + geom_point() + theme_bw() +
  xlab("age (milion of years)") + ylab("Pearson correlation value") + 
  ggtitle("Pearson correlation of tissue specificity score for human and mouse paralogs") +
  scale_size("Number of\ngenes")
