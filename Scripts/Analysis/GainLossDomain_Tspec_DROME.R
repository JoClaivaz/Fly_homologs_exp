'''
Joaquim Claivaz
170926

Consider only ortholog 1:1
For each complete case (for each species: in control or in one domain modification group and complete expresssion data)
determine whether an ortholog DROME gain or loss domain and see effect on expression
Determine also Tspec
'''
#require library
library(tidyr)
#

#FUN
data_organization_GL = function(logFC_table){
  logFC_data = read.table(logFC_table, sep = '\t')
  logFC_data = separate(logFC_data, col = 'genes', into = c('sp1', 'sp2'), sep = '_')
  logFC_data$genes = logFC_data$sp1
  for (row_tmp in 1:dim(logFC_data)[1]){
    if (grepl('DROME', logFC_data$sp1[row_tmp]) == FALSE){
      logFC_data$genes[row_tmp] = logFC_data$sp2[row_tmp]
    }
  }
  logFC_data$sp1 = NULL
  logFC_data$sp2 = NULL
  logFC_data$species = strsplit(strsplit(logFC_table, 'logFC_drome_')[[1]][2], '_QL')[[1]][1]
  
  return(logFC_data)
}

extraction_name_specie2 = function(specie2_table){
  specie2_pair = read.table(specie2_table)
  specie2_pair$drome = NA
  specie2_name = tolower(strsplit(strsplit(specie2_table, 'pair_DROME_')[[1]][2], '_domain')[[1]][1])
  specie2_pair$specie2_name = NA
  
  for (ortho_line in 1:dim(specie2_pair)[1]){
    if (grepl('DROME', specie2_pair$V1[ortho_line])){
      specie2_pair$drome[ortho_line] = as.character(specie2_pair$V1[ortho_line])
      specie2_pair$specie2_name[ortho_line] = as.character(specie2_pair$V3[ortho_line])
    }else if (grepl('DROME', specie2_pair$V3[ortho_line])){
      specie2_pair$drome[ortho_line] = as.character(specie2_pair$V3[ortho_line])
      specie2_pair$specie2_name[ortho_line] = as.character(specie2_pair$V1[ortho_line])
    }
  }
  specie2_pair$V1 = NULL
  specie2_pair$V2 = NULL
  specie2_pair$V3 = NULL
  specie2_pair$V4 = NULL
  
  return(specie2_pair)
}

log_transformation_rpkm = function(data_frame){
  for (i in 2:dim(data_frame)[2]){
    data_frame[,i] = log2(data_frame[,i] + 0.000001)
    data_frame[,i][data_frame[,i] < 1] = 0
  }
  return(data_frame)
}

#tau calculation from Nadeza work
fTau <- function(x){
  if(all(!is.na(x)))  {
    if(min(x, na.rm=TRUE) >= 0)    {
      if(max(x)!=0)      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
    } 
  } else {
    res <- NA
  } 
  return(res)
}
#

####Gain/loss domain effect on expression data####
#Data organization
droan_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droan_QL')
dromo_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_dromo_QL')
drops_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drops_QL')
drosi_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drosi_QL')
drovi_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drovi_QL')
droya_logFC = data_organization_GL('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droya_QL')

species_logFC = rbind(droan_logFC, dromo_logFC, drops_logFC, drosi_logFC, drovi_logFC, droya_logFC)
species_logFC$species = as.factor(species_logFC$species)
#

#Loading dataset of common gene of all species
common_genes = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/all_species_genes', 
                          sep = '\t', header = TRUE)
#

#Loading ID converter
ID_converter = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/DROME_id_converter', 
                          sep = '\t')
#

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
#

#Intersect among common_dataset and species_logFC (consider only ortholog 1:1 between pairwise species)
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

common_dataset$status2[common_dataset$status == as.factor(2)] = as.character('modif')
common_dataset$status2[common_dataset$status == as.factor(1)] = as.character('control')
common_dataset$status = common_dataset$status2
common_dataset$status2 = NULL
#

#keep complete case
common_dataset = common_dataset[complete.cases(common_dataset),]
#

#Verify if whole species have expression data for each DROME gene
common_test = aggregate(logFC ~ DROME_genes, data = common_dataset, paste, collapse = ' ')
common_test = separate(common_test, sep = ' ', col = logFC, into = as.character(c(1:12)))
common_test = common_test[complete.cases(common_test),]
common_test = common_test$DROME_genes
#

#Keep only gene having complete expression data for whole species
common_dataset_completecase = common_dataset[common_dataset$DROME_genes %in% common_test,]
#

#Extraction of species2 gene name
common_dataset_completecase_modif = common_dataset_completecase[common_dataset_completecase$status == 'modif',]

droan_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/final_pair_DROME_DROAN_domain_loss')
droan_pair$droan = droan_pair$specie2_name
droan_pair$specie2_name = NULL
dromo_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/final_pair_DROME_DROMO_domain_loss')
dromo_pair$dromo = dromo_pair$specie2_name
dromo_pair$specie2_name = NULL
drops_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/final_pair_DROME_DROPS_domain_loss')
drops_pair$drops = drops_pair$specie2_name
drops_pair$specie2_name = NULL
drosi_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/final_pair_DROME_DROSI_domain_loss')
drosi_pair$drosi = drosi_pair$specie2_name
drosi_pair$specie2_name = NULL
drovi_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/final_pair_DROME_DROVI_domain_loss')
drovi_pair$drovi = drovi_pair$specie2_name
drovi_pair$specie2_name = NULL
droya_pair = extraction_name_specie2('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/final_pair_DROME_DROYA_domain_loss')
droya_pair$droya = droya_pair$specie2_name
droya_pair$specie2_name = NULL

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
#

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

#Analysis
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '1'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC in function of gain/loss of domain in DROME', xlab = 'logFC value', cex.main = 0.8, xlim = c(-15, 10))
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '-1'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("modification +1", "modification -1"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
model_modif = lm(logFC ~ as.factor(length_dif) + sex, data = common_dataset_completecase_modif)
anova(model_modif)
#no effect of sex, length highly significative

####Tspecificity analysis####
#store gene list with at least one modification in pairwise comparison
modified_gene = unique(common_dataset_completecase$FB_genes[common_dataset_completecase$status == 'modif'])
#

#load expression dataset from flybase of DROME
flybase_expression = read.table('D:/UNIL/Master/Master_Project/Data/flybase/gene_rpkm_report_fb_2017_04_formated',
                                sep = '\t', header = TRUE)
#

#filter out not considered genes
flybase_expression = flybase_expression[flybase_expression$FBgn %in% common_dataset_completecase$FB_genes,]
#

#keep only parent library name: modENCODE_mRNA-Seq_tissues
flybase_expression = flybase_expression[flybase_expression$Parent_library_name == 'modENCODE_mRNA-Seq_tissues',]
#





#see different state available (29 on 124) and filter out
unique(flybase_expression$RNASource_name)
#chosen one: development time 4d, without VirF
keep_state = grepl('4d', unique(as.character(flybase_expression$RNASource_name)))
flybase_expression = flybase_expression[flybase_expression$RNASource_name %in% 
                                          unique(as.character(flybase_expression$RNASource_name))[keep_state],]
notkeep_state = grepl('VirF', unique(as.character(flybase_expression$RNASource_name)))
keep_state = grepl(FALSE, notkeep_state)
flybase_expression = flybase_expression[flybase_expression$RNASource_name %in% 
                                          unique(as.character(flybase_expression$RNASource_name))[keep_state],]
flybase_expression$RNASource_name = as.character(flybase_expression$RNASource_name)
flybase_expression$RNASource_name[flybase_expression$RNASource_name == 'mE_mRNA_A_4d_dig_sys'] = 'mE_mRNA_A_UniS_4d_digsys'
flybase_expression$RNASource_name[flybase_expression$RNASource_name == 'mE_mRNA_A_4d_carcass'] = 'mE_mRNA_A_UniS_4d_carcass'
flybase_expression$RNASource_name[flybase_expression$RNASource_name == 'mE_mRNA_A_MateM_4d_acc_gland'] = 'mE_mRNA_A_MateM_4d_accgland'

#data organization
flybase_expression = separate(flybase_expression, col = 'RNASource_name', sep = '_', into = c('x1', 'x2', 'x3', 'Sex', 'x4', 'tissue'))
flybase_expression$Sex = as.factor(flybase_expression$Sex)
flybase_expression$tissue = as.factor(flybase_expression$tissue)
str(flybase_expression)

flybase_expression = aggregate(RPKM_value ~ FBgn + tissue + Sex, data = flybase_expression, paste, collapse = ' ')
flybase_expression_unisex = flybase_expression[flybase_expression$Sex == 'UniS',]
flybase_expression_male = flybase_expression[flybase_expression$Sex == 'MateM',]
flybase_expression_female = flybase_expression[flybase_expression$Sex == 'MateF',]

flybase_expression_male = rbind(flybase_expression_male, flybase_expression_unisex)
flybase_expression_female = rbind(flybase_expression_female, flybase_expression_unisex)
flybase_expression_male$Sex = NULL
flybase_expression_female$Sex = NULL
str(flybase_expression_female)

flybase_expression_male = aggregate(RPKM_value ~ FBgn, data = flybase_expression_male, paste, collapse = '_')
flybase_expression_male = separate(flybase_expression_male, col = 'RPKM_value', sep = '_', into = c('accgland', 'head', 'testis', 'carcass', 'digsys'))
flybase_expression_female = aggregate(RPKM_value ~ FBgn, data = flybase_expression_female, paste, collapse = '_')
flybase_expression_female = separate(flybase_expression_female, col = 'RPKM_value', sep = '_', into = c('head', 'ovary', 'carcass', 'digsys'))
flybase_expression_female[,2] = as.numeric(flybase_expression_female[,2])
flybase_expression_female[,3] = as.numeric(flybase_expression_female[,3])
flybase_expression_female[,4] = as.numeric(flybase_expression_female[,4])
flybase_expression_female[,5] = as.numeric(flybase_expression_female[,5])
flybase_expression_male[,2] = as.numeric(flybase_expression_male[,2])
flybase_expression_male[,3] = as.numeric(flybase_expression_male[,3])
flybase_expression_male[,4] = as.numeric(flybase_expression_male[,4])
flybase_expression_male[,5] = as.numeric(flybase_expression_male[,5])
flybase_expression_male[,6] = as.numeric(flybase_expression_male[,6])
str(flybase_expression_male)
str(flybase_expression_female)

#Tspec calculation
flybase_expression_female = log_transformation_rpkm(flybase_expression_female)
flybase_expression_female$tspec = NA
for (gene_row in 1:dim(flybase_expression_female)[1]){
  flybase_expression_female$tspec[gene_row] = fTau(flybase_expression_female[gene_row, 3:dim(flybase_expression_female)[2]-1])
}

flybase_expression_male = log_transformation_rpkm(flybase_expression_male)
flybase_expression_male$tspec = NA
for (gene_row in 1:dim(flybase_expression_male)[1]){
  flybase_expression_male$tspec[gene_row] = fTau(flybase_expression_male[gene_row, 3:dim(flybase_expression_male)[2]-1])
}

#status inference of tspec: ubiquitous or specific (threshold 0.8)
flybase_expression_female$status = 'ubiquitous'
flybase_expression_male$status = 'ubiquitous'
flybase_expression_female$status[flybase_expression_female$tspec > 0.8] = 'specific'
flybase_expression_male$status[flybase_expression_male$tspec > 0.8] = 'specific'
flybase_expression_female$status = as.factor(flybase_expression_female$status)
flybase_expression_male$status = as.factor(flybase_expression_male$status)

#for tissue specific gene infere which tissue specificity
flybase_expression_female$specificity = 'ubiquitous'
for (gene_row in 1:dim(flybase_expression_female)[1]){
  if (flybase_expression_female$status[gene_row] == 'specific'){
    flybase_expression_female$specificity[gene_row] = names(flybase_expression_female[gene_row,2:5])[which(flybase_expression_female[gene_row,2:5] == max(flybase_expression_female[gene_row,2:5]))] 
  }
}

flybase_expression_male$specificity = 'ubiquitous'
for (gene_row in 1:dim(flybase_expression_male)[1]){
  if (flybase_expression_male$status[gene_row] == 'specific'){
    flybase_expression_male$specificity[gene_row] = names(flybase_expression_male[gene_row,2:6])[which(flybase_expression_male[gene_row,2:6] == max(flybase_expression_male[gene_row,2:6]))] 
  }
}

###analysis
#dataset modification and control
flybase_expression_male_modified = flybase_expression_male[flybase_expression_male$FBgn %in% modified_gene,]
flybase_expression_male_control = flybase_expression_male[!(flybase_expression_male$FBgn %in% modified_gene),]
flybase_expression_female_modified = flybase_expression_female[flybase_expression_female$FBgn %in% modified_gene,]
flybase_expression_female_control = flybase_expression_female[!(flybase_expression_female$FBgn %in% modified_gene),]

hist(flybase_expression_male_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in fly male in function of the status', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_male_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topleft', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(flybase_expression_female_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in fly female in function of the status', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_female_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('topleft', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

fb_female_c = table(flybase_expression_female_control$specificity)
fb_female_m = table(flybase_expression_female_modified$specificity)
fb_male_c = table(flybase_expression_male_control$specificity)
fb_male_m = table(flybase_expression_male_modified$specificity)

col = palette()
considered_table = fb_male_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in control male genes", col = col[1:6])

considered_table = fb_male_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in modified male genes", col = col[1:6])

considered_table = fb_female_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in control female genes", col = col[2:6])

considered_table = fb_female_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in modified female genes", col = col[2:6])

#test anova application
model_test = lm(c(flybase_expression_female_control$tspec, flybase_expression_female_modified$tspec) ~ 
                  as.factor(c(rep(1, length(flybase_expression_female_control$tspec)),
                              rep(0, length(flybase_expression_female_modified$tspec)))))
qqnorm(residuals(model_test)); qqline(residuals(model_test))

kruskal.test(c(flybase_expression_female_control$tspec, flybase_expression_female_modified$tspec) ~ 
               as.factor(c(rep(1, length(flybase_expression_female_control$tspec)),
                           rep(0, length(flybase_expression_female_modified$tspec)))))

kruskal.test(c(flybase_expression_male_control$tspec, flybase_expression_male_modified$tspec) ~ 
               as.factor(c(rep(1, length(flybase_expression_male_control$tspec)),
                           rep(0, length(flybase_expression_male_modified$tspec)))))

#table proportion
fb_female_m_p = fb_female_m / sum(fb_female_m)
fb_female_c_p = fb_female_c / sum(fb_female_c)
fb_male_m_p = fb_male_m / sum(fb_male_m)
fb_male_c_p = fb_male_c / sum(fb_male_c)

#chi-square test: comparison between control and modification 
n_female = sum(fb_female_c +  fb_female_m)
n_female_c = sum(fb_female_c)
n_female_m = sum(fb_female_m)
p_female = n_female_m / (n_female_c + n_female_m)
s_female = fb_female_c +  fb_female_m
exp_female_m = s_female * p_female
exp_female_c = s_female * (1 - p_female)

chi_female = sum((fb_female_m - exp_female_m)^2 / exp_female_m) + sum((fb_female_c - exp_female_c)^2 / exp_female_c)
df_female = length(fb_female_c) - 1
pchisq(chi_female, df = df_female, lower.tail = F)

n_male = sum(fb_male_c +  fb_male_m)
n_male_c = sum(fb_male_c)
n_male_m = sum(fb_male_m)
p_male = n_male_m / (n_male_c + n_male_m)
s_male = fb_male_c +  fb_male_m
exp_male_m = s_male * p_male
exp_male_c = s_male * (1 - p_male)

chi_male = sum((fb_male_m - exp_male_m)^2 / exp_male_m) + sum((fb_male_c - exp_male_c)^2 / exp_male_c)
df_male = length(fb_male_c) - 1
pchisq(chi_male, df = df_male, lower.tail = F)

#in female, modification occured differentially in function of the tissue, some tissues are more subject to domain modification
#not in male