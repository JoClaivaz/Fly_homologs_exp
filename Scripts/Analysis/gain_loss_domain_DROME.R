'''
Joaquim Claivaz
170905

For each complete case (for each species: in control or in one domain modification group, and complete expresssion data)
determine whether an ortholog DROME gain or loss domain
'''
#require library
library(tidyr)

###beginining expression_analysis_all_species.R, to have common_dataset_completecase
###START
droan_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droan_QL', sep = '\t')
dromo_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_dromo_QL', sep = '\t')
drops_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drops_QL', sep = '\t')
drosi_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drosi_QL', sep = '\t')
drovi_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drovi_QL', sep = '\t')
droya_logFC = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droya_QL', sep = '\t')

droan_logFC = separate(droan_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
dromo_logFC = separate(dromo_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
drops_logFC = separate(drops_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
drosi_logFC = separate(drosi_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
drovi_logFC = separate(drovi_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
droya_logFC = separate(droya_logFC, col = 'genes', into = c('sp1', 'sp2'), sep = '_')

droan_logFC$genes = droan_logFC$sp1
for (row_tmp in 1:dim(droan_logFC)[1]){
  if (grepl('DROME', droan_logFC$sp1[row_tmp]) == FALSE){
    droan_logFC$genes[row_tmp] = droan_logFC$sp2[row_tmp]
  }
}
droan_logFC$sp1 = NULL
droan_logFC$sp2 = NULL

dromo_logFC$genes = dromo_logFC$sp1
for (row_tmp in 1:dim(dromo_logFC)[1]){
  if (grepl('DROME', dromo_logFC$sp1[row_tmp]) == FALSE){
    dromo_logFC$genes[row_tmp] = dromo_logFC$sp2[row_tmp]
  }
}
dromo_logFC$sp1 = NULL
dromo_logFC$sp2 = NULL

drops_logFC$genes = drops_logFC$sp1
for (row_tmp in 1:dim(drops_logFC)[1]){
  if (grepl('DROME', drops_logFC$sp1[row_tmp]) == FALSE){
    drops_logFC$genes[row_tmp] = drops_logFC$sp2[row_tmp]
  }
}
drops_logFC$sp1 = NULL
drops_logFC$sp2 = NULL

drosi_logFC$genes = drosi_logFC$sp1
for (row_tmp in 1:dim(drosi_logFC)[1]){
  if (grepl('DROME', drosi_logFC$sp1[row_tmp]) == FALSE){
    drosi_logFC$genes[row_tmp] = drosi_logFC$sp2[row_tmp]
  }
}
drosi_logFC$sp1 = NULL
drosi_logFC$sp2 = NULL

drovi_logFC$genes = drovi_logFC$sp1
for (row_tmp in 1:dim(drovi_logFC)[1]){
  if (grepl('DROME', drovi_logFC$sp1[row_tmp]) == FALSE){
    drovi_logFC$genes[row_tmp] = drovi_logFC$sp2[row_tmp]
  }
}
drovi_logFC$sp1 = NULL
drovi_logFC$sp2 = NULL

droya_logFC$genes = droya_logFC$sp1
for (row_tmp in 1:dim(droya_logFC)[1]){
  if (grepl('DROME', droya_logFC$sp1[row_tmp]) == FALSE){
    droya_logFC$genes[row_tmp] = droya_logFC$sp2[row_tmp]
  }
}
droya_logFC$sp1 = NULL
droya_logFC$sp2 = NULL

droan_logFC$species = 'droan'
dromo_logFC$species = 'dromo'
drops_logFC$species = 'drops'
drosi_logFC$species = 'drosi'
drovi_logFC$species = 'drovi'
droya_logFC$species = 'droya'

species_logFC = rbind(droan_logFC, dromo_logFC, drops_logFC, drosi_logFC, drovi_logFC, droya_logFC)
species_logFC$species = as.factor(species_logFC$species)

#Loading dataset of common gene of all species
common_genes = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/all_species_genes', 
                          sep = '\t', header = TRUE)

#Loading ID converter
ID_converter = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/DROME_id_converter', 
                          sep = '\t')

#Add in common_dataset DROME identifier
common_dataset = data.frame(FB_genes = rep(common_genes$DROME_geneID, 6 * 2), 
                            DROME_genes = rep(NA, 6 * 2 * dim(common_genes)[1]),
                            status = rep(NA, 6 * 2 * dim(common_genes)[1]), 
                            sex = rep(c(rep('male', dim(common_genes)[1]), rep('female',dim(common_genes)[1])), 6), 
                            species = rep(c('droan', 'dromo', 'drops', 'drosi', 'drovi', 'droya'), each = 2 * dim(common_genes)[1]), 
                            logFC = rep(NA, 6 * 2 * dim(common_genes)[1]))

for (com_gene in 1:dim(common_dataset)[1]){
  common_dataset$DROME_genes[com_gene] = as.character(ID_converter$V1[ID_converter$V2 == as.character(common_dataset$FB_genes[com_gene])])
}

#Intersect among common_dataset and species_logFC
for (com_gene in 1:dim(common_dataset)[1]){
  if(length(species_logFC$logFC[(species_logFC$sex == common_dataset$sex[com_gene]) &
                                (species_logFC$species == common_dataset$species[com_gene]) &
                                (species_logFC$genes == common_dataset$DROME_genes[com_gene])]) == 1){
    common_dataset$logFC[com_gene] = species_logFC$logFC[(species_logFC$sex == common_dataset$sex[com_gene]) &
                                                           (species_logFC$species == common_dataset$species[com_gene]) &
                                                           (species_logFC$genes == common_dataset$DROME_genes[com_gene])]
    common_dataset$status[com_gene] = species_logFC$status[(species_logFC$sex == common_dataset$sex[com_gene]) &
                                                             (species_logFC$species == common_dataset$species[com_gene]) &
                                                             (species_logFC$genes == common_dataset$DROME_genes[com_gene])]
  }
}

length(which(common_dataset$status == as.factor(2)))

common_dataset$status2[common_dataset$status == as.factor(2)] = as.character('modif')
common_dataset$status2[common_dataset$status == as.factor(1)] = as.character('control')
common_dataset$status = common_dataset$status2
common_dataset$status2 = NULL

#keep complete case
common_dataset = common_dataset[complete.cases(common_dataset),]

#Verify if whole species have expression data for each DROME gene
str(common_dataset)
common_test = aggregate(logFC ~ DROME_genes, data = common_dataset, paste, collapse = ' ')
common_test = separate(common_test, sep = ' ', col = logFC, into = as.character(c(1:12)))
common_test = common_test[complete.cases(common_test),]
common_test = common_test$DROME_genes

#Keep only gene having complete expression data for whole species
common_dataset_completecase = common_dataset[common_dataset$DROME_genes %in% common_test,]
###END

#Extraction of species2 gene name
common_dataset_completecase_modif = common_dataset_completecase[common_dataset_completecase$status == 'modif',]


droan_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/final_pair_DROME_DROAN_domain_loss')
droan_pair$drome = NA
droan_pair$droan = NA
for (ortho_line in 1:dim(droan_pair)[1]){
  if (grepl('DROME', droan_pair$V1[ortho_line])){
    droan_pair$drome[ortho_line] = as.character(droan_pair$V1[ortho_line])
    droan_pair$droan[ortho_line] = as.character(droan_pair$V3[ortho_line])
  }else if (grepl('DROME', droan_pair$V3[ortho_line])){
    droan_pair$drome[ortho_line] = as.character(droan_pair$V3[ortho_line])
    droan_pair$droan[ortho_line] = as.character(droan_pair$V1[ortho_line])
  }
}
droan_pair$V1 = NULL
droan_pair$V2 = NULL
droan_pair$V3 = NULL
droan_pair$V4 = NULL

dromo_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/final_pair_DROME_DROMO_domain_loss')
dromo_pair$drome = NA
dromo_pair$dromo = NA
for (ortho_line in 1:dim(dromo_pair)[1]){
  if (grepl('DROME', dromo_pair$V1[ortho_line])){
    dromo_pair$drome[ortho_line] = as.character(dromo_pair$V1[ortho_line])
    dromo_pair$dromo[ortho_line] = as.character(dromo_pair$V3[ortho_line])
  }else if (grepl('DROME', dromo_pair$V3[ortho_line])){
    dromo_pair$drome[ortho_line] = as.character(dromo_pair$V3[ortho_line])
    dromo_pair$dromo[ortho_line] = as.character(dromo_pair$V1[ortho_line])
  }
}
dromo_pair$V1 = NULL
dromo_pair$V2 = NULL
dromo_pair$V3 = NULL
dromo_pair$V4 = NULL

drops_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/final_pair_DROME_DROPS_domain_loss')
drops_pair$drome = NA
drops_pair$drops = NA
for (ortho_line in 1:dim(drops_pair)[1]){
  if (grepl('DROME', drops_pair$V1[ortho_line])){
    drops_pair$drome[ortho_line] = as.character(drops_pair$V1[ortho_line])
    drops_pair$drops[ortho_line] = as.character(drops_pair$V3[ortho_line])
  }else if (grepl('DROME', drops_pair$V3[ortho_line])){
    drops_pair$drome[ortho_line] = as.character(drops_pair$V3[ortho_line])
    drops_pair$drops[ortho_line] = as.character(drops_pair$V1[ortho_line])
  }
}
drops_pair$V1 = NULL
drops_pair$V2 = NULL
drops_pair$V3 = NULL
drops_pair$V4 = NULL

drosi_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/final_pair_DROME_DROSI_domain_loss')
drosi_pair$drome = NA
drosi_pair$drosi = NA
for (ortho_line in 1:dim(drosi_pair)[1]){
  if (grepl('DROME', drosi_pair$V1[ortho_line])){
    drosi_pair$drome[ortho_line] = as.character(drosi_pair$V1[ortho_line])
    drosi_pair$drosi[ortho_line] = as.character(drosi_pair$V3[ortho_line])
  }else if (grepl('DROME', drosi_pair$V3[ortho_line])){
    drosi_pair$drome[ortho_line] = as.character(drosi_pair$V3[ortho_line])
    drosi_pair$drosi[ortho_line] = as.character(drosi_pair$V1[ortho_line])
  }
}
drosi_pair$V1 = NULL
drosi_pair$V2 = NULL
drosi_pair$V3 = NULL
drosi_pair$V4 = NULL

drovi_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/final_pair_DROME_DROVI_domain_loss')
drovi_pair$drome = NA
drovi_pair$drovi = NA
for (ortho_line in 1:dim(drovi_pair)[1]){
  if (grepl('DROME', drovi_pair$V1[ortho_line])){
    drovi_pair$drome[ortho_line] = as.character(drovi_pair$V1[ortho_line])
    drovi_pair$drovi[ortho_line] = as.character(drovi_pair$V3[ortho_line])
  }else if (grepl('DROME', drovi_pair$V3[ortho_line])){
    drovi_pair$drome[ortho_line] = as.character(drovi_pair$V3[ortho_line])
    drovi_pair$drovi[ortho_line] = as.character(drovi_pair$V1[ortho_line])
  }
}
drovi_pair$V1 = NULL
drovi_pair$V2 = NULL
drovi_pair$V3 = NULL
drovi_pair$V4 = NULL

droya_pair = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/final_pair_DROME_DROYA_domain_loss')
droya_pair$drome = NA
droya_pair$droya = NA
for (ortho_line in 1:dim(droya_pair)[1]){
  if (grepl('DROME', droya_pair$V1[ortho_line])){
    droya_pair$drome[ortho_line] = as.character(droya_pair$V1[ortho_line])
    droya_pair$droya[ortho_line] = as.character(droya_pair$V3[ortho_line])
  }else if (grepl('DROME', droya_pair$V3[ortho_line])){
    droya_pair$drome[ortho_line] = as.character(droya_pair$V3[ortho_line])
    droya_pair$droya[ortho_line] = as.character(droya_pair$V1[ortho_line])
  }
}
droya_pair$V1 = NULL
droya_pair$V2 = NULL
droya_pair$V3 = NULL
droya_pair$V4 = NULL

for (ortho_line in 1:dim(common_dataset_completecase_modif)[1]){
  if (common_dataset_completecase_modif$species[ortho_line] == 'droan'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = droan_pair$droan[droan_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'dromo'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = dromo_pair$dromo[dromo_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drops'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = drops_pair$drops[drops_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drosi'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = drosi_pair$drosi[drosi_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drovi'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = drovi_pair$drovi[drovi_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'droya'){
    common_dataset_completecase_modif$SPECIE_gene[ortho_line] = droya_pair$droya[droya_pair$drome == common_dataset_completecase_modif$DROME_genes[ortho_line]]
  }
}


#Extraction of domain inferred by pfamscan
common_dataset_completecase_modif$length_dif = NA

droan_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/ortholog_DROME_DROAN_domain')
dromo_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/ortholog_DROME_DROMO_domain')
drops_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/ortholog_DROME_DROPS_domain')
drosi_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/ortholog_DROME_DROSI_domain')
drovi_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/ortholog_DROME_DROVI_domain')
droya_domain = read.table('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/ortholog_DROME_DROYA_domain')

for (ortho_line in 1:dim(common_dataset_completecase_modif)[1]){
  if (common_dataset_completecase_modif$species[ortho_line] == 'droan'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(droan_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(droan_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'dromo'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(dromo_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(dromo_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drops'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(drops_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(drops_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drosi'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(drosi_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(drosi_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'drovi'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(drovi_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(drovi_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }else if (common_dataset_completecase_modif$species[ortho_line] == 'droya'){
    common_dataset_completecase_modif$length_dif[ortho_line] = 
      length(which(droya_domain$V1 == common_dataset_completecase_modif$DROME_genes[ortho_line]))-length(which(droya_domain$V1 == common_dataset_completecase_modif$SPECIE_gene[ortho_line]))
  }
}

common_dataset_completecase_modif$length_dif = as.factor(common_dataset_completecase_modif$length_dif)
table(common_dataset_completecase_modif$length_dif)

#without z-score normalization
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '1'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC in function of gain/loss of domain in DROME', xlab = 'logFC value', cex.main = 0.8, xlim = c(-15, 10))
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '-1'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("modification +1", "modification -1"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

#stats
model_modif = lm(logFC ~ as.factor(length_dif), data = common_dataset_completecase_modif)
qqnorm(residuals(model_modif)); qqline(residuals(model_modif))
anova(model_modif)
kruskal.test(logFC ~ as.factor(length_dif), data = common_dataset_completecase_modif)

#with z-score normalization
common_dataset_completecase_modif_zscore = common_dataset_completecase_modif
mean_droan = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droan'])
sd_droan = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droan'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droan'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droan'] 
                                                                                                   - mean_droan) / sd_droan
mean_dromo = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'dromo'])
sd_dromo = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'dromo'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'dromo'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'dromo'] 
                                                                                                   - mean_dromo) / sd_dromo
mean_drops = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drops'])
sd_drops = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drops'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drops'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drops'] 
                                                                                                   - mean_drops) / sd_drops
mean_drosi = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drosi'])
sd_drosi = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drosi'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drosi'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drosi'] 
                                                                                                   - mean_drosi) / sd_drosi
mean_drovi = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drovi'])
sd_drovi = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drovi'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drovi'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'drovi'] 
                                                                                                   - mean_drovi) / sd_drovi
mean_droya = mean(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droya'])
sd_droya = sd(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droya'])
common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droya'] = (common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$species == 'droya']
                                                                                                   - mean_droya) / sd_droya

hist(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$length_dif == '1'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC in function of gain/loss of domain in DROME (z-score normalization)', xlab = 'logFC value', cex.main = 0.8, xlim = c(-8, 8))
hist(common_dataset_completecase_modif_zscore$logFC[common_dataset_completecase_modif_zscore$length_dif == '-1'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("modification +1", "modification -1"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

#stats
model_modif_zscore = lm(logFC ~ as.factor(length_dif), data = common_dataset_completecase_modif_zscore)
qqnorm(residuals(model_modif_zscore)); qqline(residuals(model_modif_zscore))
anova(model_modif_zscore)
kruskal.test(logFC ~ as.factor(length_dif), data = common_dataset_completecase_modif_zscore)

#save the data
save.image(file = 'C:/Users/Claivaz/Desktop/gain_loss_domain_DROME_data')