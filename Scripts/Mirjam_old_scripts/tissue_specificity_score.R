
## Tissue specificity score calculated for each experiment (GSE..)
# for future analysis only the tissue specificity score calculated with the experiment 
# containing the higher number of tissues are used.The better way to calculate the tissue specificity score
# using different experiment should be tested in the future

directory <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/"
setwd(directory)


# Tissue specific score function ----------------------------------------
## tau calculation function
tau <- function(rpkm){
  score <- sum(1 - (rpkm / max(rpkm))) / (length(rpkm) - 1)
  return(score)
}

## computation function
scoreCalculation <- function(experiment.list, directory){
  for (gse in experiment.list){
    load(paste0(directory, gse))
    short_gene <- unique(tissues.list[[1]]$Ensembl_Gene_ID_short)
    long_gene <- unique(tissues.list[[1]]$Ensembl_Gene_ID_long)
    
    tissue_spec_score_short <- NULL
    for (gene in short_gene){
      gene_expression <- NULL
      for (tissue in names(tissues.list)){
        rpkm <- mean(tissues.list[[tissue]]$RPKM_short[tissues.list[[tissue]]$Ensembl_Gene_ID_short==gene])
        gene_expression <- c(gene_expression,rpkm)
      }
      score <- tau(gene_expression)
      tissue_spec_score_short <- c(tissue_spec_score_short, score)
    }
    
    tissue_spec_score_long <- NULL
    for (gene in long_gene){
      gene_expression <- NULL
      for (tissue in names(tissues.list)){
        rpkm <- mean(tissues.list[[tissue]]$RPKM_long[tissues.list[[tissue]]$Ensembl_Gene_ID_long==gene])
        gene_expression <- c(gene_expression,rpkm)
      }
      score <- tau(gene_expression)
      tissue_spec_score_long <- c(tissue_spec_score_long, score)
    }
    if(!dir.exists(paste0(directory, "tissue_spec_score"))){
      dir.create(paste0(directory, "tissue_spec_score"))
    }
    tissue_spec_score <- data.frame(Gene.ID = c(long_gene,short_gene), 
                                    c(tissue_spec_score_long, tissue_spec_score_short)
                                    ) [complete.cases(data.frame(Gene.ID = c(long_gene, short_gene),
                                      c(tissue_spec_score_long, tissue_spec_score_short))),]
    write.table(unique(tissue_spec_score) , paste0(directory, "tissue_spec_score/", 
                strsplit(gse, "[.]")[[1]][1], "spec_score", Sys.Date(), ".txt"), 
                col.names = F, row.names = F, quote = F)
    
  }
  
}



# Computation -----------------------------------------------------------------

## for human

experiment.list <- list.files("Data/Human/", pattern = "*.RData") # list of all GSE*.RData created with data_extraxtion.R
scoreCalculation(experiment.list = experiment.list, directory = "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/Data/Human/")


## for mouse

experiment.list <- list.files("Data/Mouse/", pattern = "*.RData") # list of all GSE*.RData created with data_extraxtion.R
scoreCalculation(experiment.list = experiment.list, directory = "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/Data/Mouse//")

