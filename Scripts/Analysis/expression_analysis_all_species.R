'''
Joaquim Claivaz
170725

Expression analysis (logFC pairwise, reference DROME, QL analysis) whole species
'''
####Data organization####
#Loading and organization of logFC (QL) data for all species
library(tidyr)

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
                                      (species_logFC$genes == common_dataset$DROME_genes[com_gene])]) != 0){
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
length(common_test) == length(unique(common_dataset_completecase$DROME_genes))

####Data analysis: full data & complete case####
###common_dataset analysis
##Data visualization
#whole data logFC distribution
d = density(common_dataset$logFC[!is.na(common_dataset$logFC)])
plot(d, main = 'logFC distribution across species', ylab = 'density', xlab = 'logFC value')
hist(common_dataset$logFC[!is.na(common_dataset$logFC)], breaks = 1000, freq = FALSE, add = TRUE)

#logFC distribution in function of status
hist(common_dataset$logFC[common_dataset$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = 'logFC distribution across species', ylab = 'density', xlab = 'logFC value')
hist(common_dataset$logFC[common_dataset$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topright', c('control', 'modif'), fill=c(rgb(0, 0 , 1, 0.5),rgb(1, 0 , 0, 0.5)), cex=0.8, horiz=T)

#boxplot logFC ~ species
boxplot(logFC ~ species, data = common_dataset, main = 'logFC in function of the species', ylab = 'logFC value')

#ANOVA condition
model_species_status_sex = lm(common_dataset$logFC ~ common_dataset$sex * as.factor(common_dataset$status) / common_dataset$species)
summary(model_species_status_sex)

qqnorm(residuals(model_species_status_sex)); qqline(residuals(model_species_status_sex))
#non normal dostribution
anova(model_species_status_sex)

###common_dataset_completecase analysis
##Data visualization
#whole data logFC distribution
d = density(common_dataset_completecase$logFC[!is.na(common_dataset_completecase$logFC)])
plot(d, main = 'logFC distribution across species', ylab = 'density', xlab = 'logFC value')
hist(common_dataset_completecase$logFC[!is.na(common_dataset_completecase$logFC)], breaks = 1000, freq = FALSE, add = TRUE)

#logFC distribution in function of status
hist(common_dataset_completecase$logFC[common_dataset_completecase$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = 'logFC distribution across species', ylab = 'density', xlab = 'logFC value', ylim = c(0, 0.7))
hist(common_dataset_completecase$logFC[common_dataset_completecase$status == 'control'], breaks = 200, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topright', c('control', 'modif'), fill=c(rgb(0, 0 , 1, 0.5),rgb(1, 0 , 0, 0.5)), cex=0.8, horiz=T)

#boxplot logFC ~ species
boxplot(logFC ~ status + species, data = common_dataset_completecase, 
        main = 'logFC in function of the species and the modification', ylab = 'logFC value', cex.axis = 0.7, 
        las = 2, cex = 0.4)

#ANOVA condition
model_species_status_sex = lm(common_dataset_completecase$logFC ~ common_dataset_completecase$sex * as.factor(common_dataset_completecase$status) * common_dataset_completecase$species)
summary(model_species_status_sex)

qqnorm(residuals(model_species_status_sex)); qqline(residuals(model_species_status_sex))
#non normal dostribution
anova(model_species_status_sex)

kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase)
kruskal.test(logFC ~ as.factor(sex), data = common_dataset_completecase)

kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase[common_dataset_completecase$sex == 'male',])
kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase[common_dataset_completecase$sex == 'female',])

###Zscore normalization specific to species (and avoid specie factor in analysis, but keeping sex)
##Normalization step specific to species
common_dataset_completecase_zscore = common_dataset_completecase

mean_droan = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droan'])
sd_droan = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droan'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droan'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droan'] 
                                                                                     - mean_droan) / sd_droan

mean_dromo = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'dromo'])
sd_dromo = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'dromo'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'dromo'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'dromo'] 
                                                                                     - mean_dromo) / sd_dromo

mean_drops = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drops'])
sd_drops = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drops'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drops'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drops'] 
                                                                                     - mean_drops) / sd_drops

mean_drosi = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drosi'])
sd_drosi = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drosi'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drosi'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drosi'] 
                                                                                     - mean_drosi) / sd_drosi

mean_drovi = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drovi'])
sd_drovi = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drovi'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drovi'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'drovi'] 
                                                                                     - mean_drovi) / sd_drovi

mean_droya = mean(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droya'])
sd_droya = sd(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droya'])
common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droya'] = (common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$species == 'droya'] 
                                                                                     - mean_droya) / sd_droya

#Distribution logFC data normalized                                                                                                                                                                         - mean_droya) / sd_droya
d = density(common_dataset_completecase_zscore$logFC[!is.na(common_dataset_completecase_zscore$logFC)])
plot(d, main = 'logFC distribution across species', ylab = 'density', xlab = 'logFC value')
hist(common_dataset_completecase_zscore$logFC[!is.na(common_dataset_completecase_zscore$logFC)], breaks = 1000, freq = FALSE, add = TRUE)

#logFC distribution in function of status
hist(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5),
     main = 'logFC normalized distribution', ylab = 'density', xlab = 'logFC value', ylim = c(0, 0.7))
hist(common_dataset_completecase_zscore$logFC[common_dataset_completecase_zscore$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topright', c('control', 'modif'), fill=c(rgb(0, 0 , 1, 0.5), rgb(1, 0 , 0, 0.5)), cex=0.8, horiz=T)

#boxplot logFC ~ species
boxplot(logFC ~ status + species, data = common_dataset_completecase_zscore, 
        main = 'logFC normalized in function of the species and the modification', ylab = 'logFC value', cex.axis = 0.7, 
        las = 2, cex = 0.4)

#ANOVA conditio
model_species_status_sex = lm(common_dataset_completecase_zscore$logFC ~ as.factor(common_dataset_completecase_zscore$sex) + as.factor(common_dataset_completecase_zscore$status) + as.factor(common_dataset_completecase_zscore$species))
summary(model_species_status_sex)

qqnorm(residuals(model_species_status_sex)); qqline(residuals(model_species_status_sex))
#non normal dostribution
anova(model_species_status_sex)

kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase_zscore)
kruskal.test(logFC ~ sex, data = common_dataset_completecase_zscore)

kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase_zscore[common_dataset_completecase_zscore$sex == 'male',])
kruskal.test(logFC ~ as.factor(status), data = common_dataset_completecase_zscore[common_dataset_completecase_zscore$sex == 'female',])
