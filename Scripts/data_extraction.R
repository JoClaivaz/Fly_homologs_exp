
directory <- "D:/UNIL/Master/Master_Project/Script/Mirjam/"

setwd(directory)

# Load package ------------------------------------------------------------

#source("https://bioconductor.org/biocLite.R")
#biocLite("BgeeDB")
#biocLite("biomaRt")

library(BgeeDB)
library(biomaRt)

# Load data ---------------------------------------------------------------

human <- read.table("Data/Homo_sapiens.GRCh38.pep.all.fa.liso.dom.pairs.txt")
mouse <- read.table("Data/Mus_musculus.GRCm38.pep.all.fa.liso.dom.pairs.txt")
colnames(human) <- c("long_domain", "domain_number_long", "short_domain", "domain_difference")
colnames(mouse) <- c("long_domain", "domain_number_long", "short_domain", "domain_difference")

# Rename protein ID
human$long_domain <- sapply(strsplit(as.character(human$long_domain), "[.]"), "[[", 1)
human$short_domain <- sapply(strsplit(as.character(human$short_domain), "[.]"), "[[", 1)
mouse$long_domain <- sapply(strsplit(as.character(mouse$long_domain), "[.]"), "[[", 1)
mouse$short_domain <- sapply(strsplit(as.character(mouse$short_domain), "[.]"), "[[", 1)


# BiomaRt -----------------------------------------------------------------
# define biomart object
mart_human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mart_mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Using Biomart find corresponding ensemble id for gene of protein with long and short domain

## For human

## gene id and protein id for the longer domain
human_long <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                    filters = "ensembl_peptide_id", values = human$long_domain, mart = mart_human)
colnames(human_long) <- c("Ensembl_Gene_ID_long", "Ensembl_Protein_ID_long")

## gene id and protein id for shorter domain
human_short <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                     filters = "ensembl_peptide_id", values = human$short_domain, mart = mart_human)
colnames(human_short) <- c("Ensembl_Gene_ID_short", "Ensembl_Protein_ID_short")

## link dataset of human --> merge initial data set to human_long
data_human <- merge(human, human_long, by.x = "long_domain", by.y = "Ensembl_Protein_ID_long")
## link the merged data set to human_short
data_human <- merge(data_human, human_short, by.x = "short_domain", by.y = "Ensembl_Protein_ID_short")


## For mouse

## gene id and protein id for the longer domain
mouse_long <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                    filters = "ensembl_peptide_id", values = mouse$long_domain, mart = mart_mouse)
colnames(mouse_long) <- c("Ensembl_Gene_ID_long", "Ensembl_Protein_ID_long")

## gene id and protein id for shorter domain
mouse_short <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                     filters = "ensembl_peptide_id", values = mouse$short_domain, mart = mart_mouse)
colnames(mouse_short) <- c("Ensembl_Gene_ID_short", "Ensembl_Protein_ID_short")

## link dataset of mouse --> merge initial data set to mouse_long
data_mouse <- merge(mouse, mouse_long, by.x = "long_domain", by.y = "Ensembl_Protein_ID_long")
## link the merged data set to mouse_short
data_mouse <- merge(data_mouse, mouse_short, by.x = "short_domain", by.y = "Ensembl_Protein_ID_short")


# Bgee --------------------------------------------------------------------
setwd(paste0(directory, "Data/Bgee/"))

## crete Bgee species object
bgee_human <- Bgee$new(species = "Homo_sapiens", dataType = "rna_seq")
bgee_mouse <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")

## Retrieve annotation for the choosen species and data
annotation_bgee_human <- getAnnotation(bgee_human)
annotation_bgee_mouse <- getAnnotation(bgee_mouse)

### Download the processed RNA-seq data for human and mouse: download all RPKMs and counts
data_bgee_human <- getData(bgee_human)
data_bgee_mouse <- getData(bgee_mouse)


## save data for future analysis -------------------------------------------
save(data_bgee_human, data_bgee_mouse, data_human, data_mouse, mart_human, mart_mouse, 
     file = paste0(directory, "Data/base_data", Sys.Date(), ".RData"))


# Find human-mouse ortholog in the dataset (data_human and data_mouse) --------

## gene id and protein id for the longer domain with the correspondent mouse homologs 
mouse_human_homo_long <- getLDS(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                                filters = "ensembl_peptide_id", values = human$long_domain, mart = mart_human, 
                                attributesL = c("ensembl_gene_id"), 
                                martL = mart_mouse)
colnames(mouse_human_homo_long) <- c("Human_Ensembl_Gene_ID_long", "Human_Ensembl_Protein_ID_long",
                                     "Mouse_Ensembl_Gene_ID_long" )

## gene id and protein id for shorter protein with the correspondent mouse homologs 
mouse_human_homo_short <- getLDS(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                                 filters = "ensembl_peptide_id", values = human$short_domain, mart = mart_human, 
                                 attributesL = c("ensembl_gene_id"), 
                                 martL = mart_mouse)
colnames(mouse_human_homo_short) <- c("Human_Ensembl_Gene_ID_short", "Human_Ensembl_Protein_ID_short",
                                      "Mouse_Ensembl_Gene_ID_short" )

## merge datasets
mouse_human_orthologs <- merge(human[,c("long_domain", "short_domain")], mouse_human_homo_short, 
                               by.x = "short_domain", by.y = "Human_Ensembl_Protein_ID_short")

mouse_human_orthologs <- merge(mouse_human_orthologs, mouse_human_homo_long, by.x = "long_domain", by.y = "Human_Ensembl_Protein_ID_long")

# write mouse human orthologs table:
write.table(mouse_human_orthologs, paste0(directory, "Data/Output/mouse_human_orthologs", Sys.Date(), ".txt"), quote = F, row.names = F)


# Data extraction for each experiment -------------------------------------------------
setwd(directory)
## Load data (if necessary)
load(file = "Data/base_data2017") # biomart, bgee and protein domain informations

## Function to extract, combine and save data coming from bgee and from biomart for each experiment

dataExtraction <- function(data, data.bgee, directory) {
  for (gse in c(1:length(data.bgee))){ 
    
    tissues.list <- list()
    for (tissue in unique(data.bgee[[gse]]$Anatomical.entity.name)){
      bgee <- unique(data.bgee[[gse]][data.bgee[[gse]]$`Anatomical.entity.name` == tissue,][,c("Gene.ID", "RPKM")])
      rpkm <- NULL
      for(gene in unique(bgee$Gene.ID)){
        rpkm_mean <- mean(bgee$RPKM[bgee$Gene.ID==gene])
        rpkm <- c(rpkm, rpkm_mean)
      }
      rpkm_df <- data.frame(Gene.ID = unique(bgee$Gene.ID), RPKM = rpkm)
      ## info for long protein domain
      data_long <- merge(data, rpkm_df, by.x = "Ensembl_Gene_ID_long", by.y = "Gene.ID")
      colnames(data_long)[colnames(data_long) =="RPKM"] <- "RPKM_long"
      ## info for short protein domain
      data_short <- merge(data, rpkm_df, by.x = "Ensembl_Gene_ID_short", by.y = "Gene.ID")
      colnames(data_short)[colnames(data_short) =="RPKM"] <- "RPKM_short"
      data.merged <- merge(data_long, data_short, all = T)
      # normalize rpkm value of paralogs by the long protein domain rpkm
      rpkm_par <- NULL
      for (geneID in unique(data.merged$Ensembl_Gene_ID_long)){
        paral.rpkm.norm <- data.merged$RPKM_short[data.merged$Ensembl_Gene_ID_long==geneID]/data.merged$RPKM_long[data.merged$Ensembl_Gene_ID_long==geneID]
        rpkm_par <- c(rpkm_par, paral.rpkm.norm)
      }
      data.merged$RPKM_short_norm <- rpkm_par
      
      tissues.list[[tissue]] <- data.merged
      
    }
    
    save(tissues.list, file = paste0(directory, unique(data.bgee[[gse]]$Experiment.ID), ".RData"))
    
  }
} 



## for human

dataExtraction(data = data_human, data.bgee = data_bgee_human, directory = "D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Human/")

## for Mouse

dataExtraction(data = data_mouse, data.bgee = data_bgee_mouse, directory = "D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Mouse/")


