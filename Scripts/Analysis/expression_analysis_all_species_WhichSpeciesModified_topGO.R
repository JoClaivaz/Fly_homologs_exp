'''
Joaquim Claivaz
170823

For each gene, which are common to the six drosophilla datasets, the script determined which species have the modification
'''


library(tidyr)

#Loading dataset of common gene of all species
common_genes = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/all_species_genes', 
                          sep = '\t', header = TRUE)

names(common_genes) = c('geneID', 'DROAN', 'DROMO', 'DROPS', 'DROSI', 'DROVI', 'DROYA')
common_genes[,2] = as.character(common_genes[,2])
common_genes[,3] = as.character(common_genes[,3])
common_genes[,4] = as.character(common_genes[,4])
common_genes[,5] = as.character(common_genes[,5])
common_genes[,6] = as.character(common_genes[,6])
common_genes[,7] = as.character(common_genes[,7])

common_genes[,2:7][common_genes[,2:7] != 'modif'] = NA

common_genes$DROAN[common_genes$DROAN == 'modif'] = 'DROAN'
common_genes$DROMO[common_genes$DROMO == 'modif'] = 'DROMO'
common_genes$DROPS[common_genes$DROPS == 'modif'] = 'DROPS'
common_genes$DROSI[common_genes$DROSI == 'modif'] = 'DROSI'
common_genes$DROVI[common_genes$DROVI == 'modif'] = 'DROVI'
common_genes$DROYA[common_genes$DROYA == 'modif'] = 'DROYA'

common_droan = common_genes[,c(1,2)]
common_dromo = common_genes[,c(1,3)]
common_drops = common_genes[,c(1,4)]
common_drosi = common_genes[,c(1,5)]
common_drovi = common_genes[,c(1,6)]
common_droya = common_genes[,c(1,7)]
names(common_droan)[2] = 'status'
names(common_dromo)[2] = 'status'
names(common_drops)[2] = 'status'
names(common_drosi)[2] = 'status'
names(common_drovi)[2] = 'status'
names(common_droya)[2] = 'status'

common_genes_noNA = rbind(common_droan, common_dromo, common_drops, common_drosi, common_drovi, common_droya)
common_genes_noNA = common_genes_noNA[complete.cases(common_genes_noNA),]
common_genes_noNA = aggregate(status ~ geneID, data = common_genes_noNA, paste, collapse = ' ')
table(common_genes_noNA$status)

#2 possibilities
strsplit(common_genes_noNA$status[1], ' ')[[1]]
gregexpr(' ', common_genes_noNA$status[1])[[1]]
#

common_genes_noNA$num_species = NA
for (gene_row in 1:dim(common_genes_noNA)[1]){
  common_genes_noNA$num_species[gene_row] = length(strsplit(common_genes_noNA$status[gene_row], ' ')[[1]]) 
}

table_2 = table(common_genes_noNA$status[common_genes_noNA$num_species == 2])
table_2 = table_2 / sum(table_2)
t_ext = which(table_2 == max(table_2))
t_ext
table_test = table_2
t.value = (mean(table_test) - t_ext) / (sd(table_test) / sqrt(length(table_test))) 
2 * pt(-abs(t.value), df = length(table_test) - 1)
#Groupe1 = DROMO & DROVI

for (gene_row in 1:dim(common_genes_noNA)[1]){
  if (all(grepl('DROMO', common_genes_noNA$status[gene_row]), 
          grepl('DROVI', common_genes_noNA$status[gene_row]))){
    common_genes_noNA$num_species[gene_row] = common_genes_noNA$num_species[gene_row] - 1 
  }
}
table_3 = table(common_genes_noNA$status[common_genes_noNA$num_species == 2])
table_3 = table_3 / sum(table_3)
t_ext = which(table_3 == max(table_3))
t_ext
table_test = table_3
t.value = (mean(table_test) - t_ext) / (sd(table_test) / sqrt(length(table_test))) 
2 * pt(-abs(t.value), df = length(table_test) - 1)
#Groupe2 = Groupe1 + DROAN

for (gene_row in 1:dim(common_genes_noNA)[1]){
  if (all(grepl('DROAN', common_genes_noNA$status[gene_row]),
          grepl('DROMO', common_genes_noNA$status[gene_row]),
          grepl('DROVI', common_genes_noNA$status[gene_row]))){
    common_genes_noNA$num_species[gene_row] = common_genes_noNA$num_species[gene_row] - 1 
  }
}
table_4 = table(common_genes_noNA$status[common_genes_noNA$num_species == 2])
table_4 = table_4 / sum(table_4)
t_ext = which(table_4 == max(table_4))
t_ext 
table_test = table_4
t.value = (mean(table_test) - t_ext) / (sd(table_test) / sqrt(length(table_test))) 
2 * pt(-abs(t.value), df = length(table_test) - 1)
#Groupe3 = Groupe2 + DROPS

for (gene_row in 1:dim(common_genes_noNA)[1]){
  if (all(grepl('DROAN', common_genes_noNA$status[gene_row]),
          grepl('DROMO', common_genes_noNA$status[gene_row]),
          grepl('DROPS', common_genes_noNA$status[gene_row]),
          grepl('DROVI', common_genes_noNA$status[gene_row]))){
    common_genes_noNA$num_species[gene_row] = common_genes_noNA$num_species[gene_row] - 1 
  }
}
table_5 = table(common_genes_noNA$status[common_genes_noNA$num_species == 2])
table_5 = table_5 / sum(table_5)
t_ext = which(table_5 == max(table_5))
t_ext
table_test = table_3
t.value = (mean(table_test) - t_ext) / (sd(table_test) / sqrt(length(table_test))) 
2 * pt(-abs(t.value), df = length(table_test) - 1)
#Groupe4 = Groupe3 + DROYA

for (gene_row in 1:dim(common_genes_noNA)[1]){
  if (all(grepl('DROAN', common_genes_noNA$status[gene_row]),
          grepl('DROMO', common_genes_noNA$status[gene_row]),
          grepl('DROPS', common_genes_noNA$status[gene_row]),
          grepl('DROVI', common_genes_noNA$status[gene_row]),
          grepl('DROYA', common_genes_noNA$status[gene_row]))){
    common_genes_noNA$num_species[gene_row] = common_genes_noNA$num_species[gene_row] - 1 
  }
}
table_6 = table(common_genes_noNA$status[common_genes_noNA$num_species == 2])
table_6 = table_6 / sum(table_6)
t_ext = which(table_6 == max(table_6))
t_ext
table_test = table_3
t.value = (mean(table_test) - t_ext) / (sd(table_test) / sqrt(length(table_test))) 
2 * pt(-abs(t.value), df = length(table_test) - 1)
#Groupe5 = Groupe4 + DROYA

####topGO analysis####
common_genes
common_genes$common_status = NA
for (considered_gene in 1:dim(common_genes)[1]){
  if (!is.na(any(common_genes[considered_gene,] == 'modif'))){
    common_genes$common_status[considered_gene] = 'modif'
  }
}

common_genes_universe = as.vector(common_genes$DROME_geneID)
names(common_genes_universe) = common_genes_universe
common_genes_universe[] = 0
common_genes_modif = as.vector(common_genes$DROME_geneID[!is.na(common_genes$common_status)])
names(common_genes_modif) = common_genes_modif
common_genes_universe[names(common_genes_universe) %in% names(common_genes_modif)] = 1 
common_genes_universe = as.factor(common_genes_universe)

table(common_genes_universe)



###topGO analysis###
#library for topGO analysis
library(topGO)
library(graph)
library(RBGL)
library(Rgraphviz)
#

#biological process
topGO_allspecies_bp = new("topGOdata", 
                     ontology = "BP", 
                     allGenes = common_genes_universe, 
                     geneSel = common_genes_modif,
                     nodeSize = 10, 
                     annot = annFUN.org,
                     mapping="org.Dm.eg.db", 
                     ID = "ensembl")
#summary
topGO_allspecies_bp

#enrichment test
resultFisher_allspecies_bp = runTest(topGO_allspecies_bp, algorithm = "classic", statistic = "fisher")
resultFisher_allspecies_bp

resultKS_allspecies_bp = runTest(topGO_allspecies_bp, algorithm = "classic", statistic = "ks")
resultKS.elim_allspecies_bp = runTest(topGO_allspecies_bp, algorithm = "elim", statistic = "ks")

#analysis
allRes_allspecies_bp = GenTable(topGO_allspecies_bp, classicFisher = resultFisher_allspecies_bp, 
                           classicKS = resultKS_allspecies_bp, elimKS = resultKS.elim_allspecies_bp,
                           orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)

allRes_allspecies_bp
write.csv(allRes_allspecies_bp, file = 'C:/Users/Claivaz/Desktop/allRes_multi')

pValue.classic_allspecies = score(resultKS_allspecies_bp)
pValue.elim_allspecies = score(resultKS.elim_allspecies_bp)[names(pValue.classic_allspecies)]
gstat_allspecies = termStat(topGO_allspecies_bp, names(pValue.classic_allspecies))
gSize_allspecies = gstat_allspecies$Annotated / max(gstat_allspecies$Annotated) * 4

#Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_allspecies = colMap(gstat_allspecies$Significant)
plot(pValue.classic_allspecies, pValue.elim_allspecies, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_allspecies, col = gCol_allspecies)

sel.go_allspecies = names(pValue.classic_allspecies)[pValue.elim_allspecies < pValue.classic_allspecies]
cbind(termStat(topGO_allspecies_bp, sel.go_allspecies),
      elim = pValue.elim_allspecies[sel.go_allspecies],
      classic = pValue.classic_allspecies[sel.go_allspecies])

showSigOfNodes(topGO_allspecies_bp, score(resultKS.elim_allspecies_bp), firstSigNodes = 5, useInfo = 'all')

#Molecular function
topGO_allspecies_mf = new("topGOdata", 
                     ontology = "MF", 
                     allGenes = common_genes_universe, 
                     geneSel = common_genes_modif,
                     nodeSize = 10, 
                     annot = annFUN.org,
                     mapping="org.Dm.eg.db", 
                     ID = "ensembl")
#summary
topGO_allspecies_mf

#enrichment test
resultFisher_allspecies_mf = runTest(topGO_allspecies_mf, algorithm = "classic", statistic = "fisher")
resultFisher_allspecies_mf

resultKS_allspecies_mf = runTest(topGO_allspecies_mf, algorithm = "classic", statistic = "ks")
resultKS.elim_allspecies_mf = runTest(topGO_allspecies_mf, algorithm = "elim", statistic = "ks")

#analysis
allRes_allspecies_mf = GenTable(topGO_allspecies_mf, classicFisher = resultFisher_allspecies_mf, 
                           classicKS = resultKS_allspecies_mf, elimKS = resultKS.elim_allspecies_mf,
                           orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)

allRes_allspecies_mf
write.csv(allRes_allspecies_mf, file = 'C:/Users/Claivaz/Desktop/allRes_multi')