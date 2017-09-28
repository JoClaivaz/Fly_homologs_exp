'''
Joaquim Claivaz
170926

For each complete case (for each species: in control or in one domain modification group and complete expresssion data)
determine whether an ortholog DROME gain or loss domain and see effect on expression
Determine also Tspec
Consider only ortholog 1:1
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

sort_RNASource_name_keep = function(dataset_flybase, regexp_pattern){
  keep_state = grepl(regexp_pattern, unique(as.character(dataset_flybase$RNASource_name)))
  dataset_flybase = dataset_flybase[dataset_flybase$RNASource_name %in% 
                                      unique(as.character(dataset_flybase$RNASource_name))[keep_state],]
  return(dataset_flybase)
}

sort_RNASource_name_notkeep = function(dataset_flybase, regexp_pattern){
  notkeep_state = grepl(regexp_pattern, unique(as.character(dataset_flybase$RNASource_name)))
  keep_state = grepl(FALSE, notkeep_state)
  dataset_flybase = dataset_flybase[dataset_flybase$RNASource_name %in% 
                                      unique(as.character(dataset_flybase$RNASource_name))[keep_state],]
  return(dataset_flybase)
}

flybase_expression_organization_bysex = function(fly_exp_cond, considered_sex, with_UniS = TRUE){
  fly_exp_cond = separate(fly_exp_cond, col = 'RNASource_name', sep = '_', into = c('x1', 'x2', 'x3', 'Sex', 'x4', 'tissue'))
  fly_exp_cond$Sex = as.factor(fly_exp_cond$Sex)
  fly_exp_cond$tissue = as.factor(fly_exp_cond$tissue)
  fly_exp_cond = aggregate(RPKM_value ~ FBgn + tissue + Sex, data = fly_exp_cond, paste, collapse = ' ')
  fly_exp_cond_tmp = fly_exp_cond[fly_exp_cond$Sex == considered_sex,]
  
  if(with_UniS == TRUE){
    fly_exp_cond_unisex = fly_exp_cond[fly_exp_cond$Sex == 'UniS',]
    fly_exp_cond = rbind(fly_exp_cond_tmp, fly_exp_cond_unisex)
  }else{
    fly_exp_cond = fly_exp_cond_tmp
  }
  
  tissue_list = as.vector(unique(fly_exp_cond$tissue))
  fly_exp_cond$Sex = NULL
  fly_exp_cond = aggregate(RPKM_value ~ FBgn, data = fly_exp_cond, paste, collapse = '_')
  fly_exp_cond = separate(fly_exp_cond, col = 'RPKM_value', sep = '_', into = tissue_list)
  
  for(col_num in 2:dim(fly_exp_cond)[2]){
    fly_exp_cond[,col_num] = as.numeric(fly_exp_cond[,col_num])
  }
  
  return(fly_exp_cond)
}

flybase_expression_organization_l3 = function(fly_exp_cond){
  fly_exp_cond = separate(fly_exp_cond, col = 'RNASource_name', sep = '_', into = c('x1', 'x2', 'x3', 'x4', 'tissue'))
  fly_exp_cond$tissue = as.factor(fly_exp_cond$tissue)
  fly_exp_cond = aggregate(RPKM_value ~ FBgn + tissue, data = fly_exp_cond, paste, collapse = ' ')
  
  tissue_list = as.vector(unique(fly_exp_cond$tissue))
  fly_exp_cond = aggregate(RPKM_value ~ FBgn, data = fly_exp_cond, paste, collapse = '_')
  fly_exp_cond = separate(fly_exp_cond, col = 'RPKM_value', sep = '_', into = tissue_list)
  
  for(col_num in 2:dim(fly_exp_cond)[2]){
    fly_exp_cond[,col_num] = as.numeric(fly_exp_cond[,col_num])
  }
  
  return(fly_exp_cond)
}

flybase_expression_organization_A_bytissue = function(fly_exp_cond, considered_tissue){
  fly_exp_cond = separate(fly_exp_cond, col = 'RNASource_name', sep = '_', into = c('x1', 'x2', 'x3', 'time', 'tissue'))
  fly_exp_cond$tissue = as.factor(fly_exp_cond$tissue)
  fly_exp_cond = aggregate(RPKM_value ~ FBgn + tissue + time, data = fly_exp_cond, paste, collapse = ' ')
  fly_exp_cond = fly_exp_cond[fly_exp_cond$tissue == considered_tissue,]
  
  time_list = as.vector(unique(fly_exp_cond$time))
  fly_exp_cond = aggregate(RPKM_value ~ FBgn, data = fly_exp_cond, paste, collapse = '_')
  fly_exp_cond = separate(fly_exp_cond, col = 'RPKM_value', sep = '_', into = time_list)
  
  for(col_num in 2:dim(fly_exp_cond)[2]){
    fly_exp_cond[,col_num] = as.numeric(fly_exp_cond[,col_num])
  }
  
  return(fly_exp_cond)
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

Tspec_inference = function(flybase_expression_data){
  flybase_expression_data = log_transformation_rpkm(flybase_expression_data)
  flybase_expression_data$tspec = NA
  for (gene_row in 1:dim(flybase_expression_data)[1]){
    flybase_expression_data$tspec[gene_row] = fTau(flybase_expression_data[gene_row, 3:dim(flybase_expression_data)[2]-1])
  }
  
  return(flybase_expression_data)
}

specificity_inference = function(flybase_expression_data, threshold_specificity = 0.8){
  column_tmp = dim(flybase_expression_data)[2]-1
  
  #status inference of tspec: ubiquitous or specific (threshold 0.8)
  flybase_expression_data$status = 'ubiquitous'
  flybase_expression_data$status[flybase_expression_data$tspec > threshold_specificity] = 'specific'
  flybase_expression_data$status = as.factor(flybase_expression_data$status)
  
  #for tissue specific gene infere which tissue specificity
  flybase_expression_data$specificity = 'ubiquitous'
  
  for (gene_row in 1:dim(flybase_expression_data)[1]){
    if (flybase_expression_data$status[gene_row] == 'specific'){
      flybase_expression_data$specificity[gene_row] = names(flybase_expression_data[gene_row,2:column_tmp])[which(flybase_expression_data[gene_row,2:column_tmp] == max(flybase_expression_data[gene_row,2:column_tmp]))] 
    }
  }
  
  return(flybase_expression_data) 
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

#Intersect among common_dataset and species_logFC (remove duplicated genes (in-paralogs) present in species2)
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
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '1' & common_dataset_completecase_modif$sex == 'female'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC in function of gain/loss of domain in DROME female', xlab = 'logFC value', cex.main = 0.8, xlim = c(-15, 10))
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '-1' & common_dataset_completecase_modif$sex == 'female'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("modification +1", "modification -1"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '1' & common_dataset_completecase_modif$sex == 'male'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC in function of gain/loss of domain in DROME male', xlab = 'logFC value', cex.main = 0.8, xlim = c(-15, 10))
hist(common_dataset_completecase_modif$logFC[common_dataset_completecase_modif$length_dif == '-1' & common_dataset_completecase_modif$sex == 'male'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("modification +1", "modification -1"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

model_modif = lm(logFC ~ as.factor(length_dif) + as.factor(sex), data = common_dataset_completecase_modif)
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
#unique(flybase_expression$RNASource_name)
#

####1.chosen state for tissue spec in function of sex: development time 4d, with VirF (virgine female ?!)####
#sort dataset and filter out not considered states
flybase_expression_4d = sort_RNASource_name_keep(dataset_flybase = flybase_expression, 
                                                 regexp_pattern = '4d')
  
flybase_expression_4d$RNASource_name = as.character(flybase_expression_4d$RNASource_name)
flybase_expression_4d$RNASource_name[flybase_expression_4d$RNASource_name == 'mE_mRNA_A_4d_dig_sys'] = 'mE_mRNA_A_UniS_4d_digsys'
flybase_expression_4d$RNASource_name[flybase_expression_4d$RNASource_name == 'mE_mRNA_A_4d_carcass'] = 'mE_mRNA_A_UniS_4d_carcass'
flybase_expression_4d$RNASource_name[flybase_expression_4d$RNASource_name == 'mE_mRNA_A_MateM_4d_acc_gland'] = 'mE_mRNA_A_MateM_4d_accgland'
#

#data organization
flybase_expression_4d_male = flybase_expression_organization_bysex(fly_exp_cond = flybase_expression_4d, 
                                                                  considered_sex = 'MateM', with_UniS = TRUE)
flybase_expression_4d_female = flybase_expression_organization_bysex(fly_exp_cond = flybase_expression_4d, 
                                                                   considered_sex = 'MateF', with_UniS = TRUE)
flybase_expression_4d_virf = flybase_expression_organization_bysex(fly_exp_cond = flybase_expression_4d, 
                                                                   considered_sex = 'Virf', with_UniS = TRUE)
#

#Tspec calculation
flybase_expression_4d_female = Tspec_inference(flybase_expression_4d_female)
flybase_expression_4d_male = Tspec_inference(flybase_expression_4d_male)
flybase_expression_4d_virf = Tspec_inference(flybase_expression_4d_virf)
#

#Inference of ubiquitous and specific gene and determine which it's the most expressed tissue for specific genes
flybase_expression_4d_female = specificity_inference(flybase_expression_data = flybase_expression_4d_female,
                                                     threshold_specificity = 0.8)
flybase_expression_4d_male = specificity_inference(flybase_expression_data = flybase_expression_4d_male,
                                                     threshold_specificity = 0.8)
flybase_expression_4d_virf = specificity_inference(flybase_expression_data = flybase_expression_4d_virf,
                                                     threshold_specificity = 0.8)
#

###analysis
#dataset modification and control
flybase_expression_4d_male_modified = flybase_expression_4d_male[flybase_expression_4d_male$FBgn %in% modified_gene,]
flybase_expression_4d_male_control = flybase_expression_4d_male[!(flybase_expression_4d_male$FBgn %in% modified_gene),]
flybase_expression_4d_female_modified = flybase_expression_4d_female[flybase_expression_4d_female$FBgn %in% modified_gene,]
flybase_expression_4d_female_control = flybase_expression_4d_female[!(flybase_expression_4d_female$FBgn %in% modified_gene),]
flybase_expression_4d_virf_modified = flybase_expression_4d_virf[flybase_expression_4d_virf$FBgn %in% modified_gene,]
flybase_expression_4d_virf_control = flybase_expression_4d_virf[!(flybase_expression_4d_virf$FBgn %in% modified_gene),]

hist(flybase_expression_4d_male_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in fly male in function of the status', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_4d_male_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(flybase_expression_4d_female_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in fly mate female in function of the status', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_4d_female_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(flybase_expression_4d_virf_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in fly virgin female in function of the status', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_4d_virf_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

fb_female_c = table(flybase_expression_4d_female_control$specificity)
fb_female_m = table(flybase_expression_4d_female_modified$specificity)
fb_male_c = table(flybase_expression_4d_male_control$specificity)
fb_male_m = table(flybase_expression_4d_male_modified$specificity)
fb_virf_c = table(flybase_expression_4d_virf_control$specificity)
fb_virf_m = table(flybase_expression_4d_virf_modified$specificity)

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
    main = "Specificity in control mate female genes", col = col[2:6])

considered_table = fb_female_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in modified mate female genes", col = col[2:6])

considered_table = fb_virf_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in control virgin female genes", col = col[2:6])

considered_table = fb_virf_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in modified virgin female genes", col = col[2:6])

#test anova application
anova_df_female = rbind(flybase_expression_4d_female_control, flybase_expression_4d_female_modified)
anova_df_female$modif_status = c(rep(0, length(flybase_expression_4d_female_control$tspec)),
                                 rep(1, length(flybase_expression_4d_female_modified$tspec)))
anova_df_female$modif_status = as.factor(anova_df_female$modif_status)

anova_df_male = rbind(flybase_expression_4d_male_control, flybase_expression_4d_male_modified)
anova_df_male$modif_status = c(rep(0, length(flybase_expression_4d_male_control$tspec)),
                                 rep(1, length(flybase_expression_4d_male_modified$tspec)))
anova_df_male$modif_status = as.factor(anova_df_male$modif_status)

anova_df_virf = rbind(flybase_expression_4d_virf_control, flybase_expression_4d_virf_modified)
anova_df_virf$modif_status = c(rep(0, length(flybase_expression_4d_virf_control$tspec)),
                                 rep(1, length(flybase_expression_4d_virf_modified$tspec)))
anova_df_virf$modif_status = as.factor(anova_df_virf$modif_status)

model_test_female = aov(anova_df_female$tspec ~ anova_df_female$modif_status)
qqnorm(residuals(model_test_female)); qqline(residuals(model_test_female))
summary(model_test_female)
kruskal.test(anova_df_female$tspec ~ anova_df_female$modif_status)

model_test_male = aov(anova_df_male$tspec ~ anova_df_male$modif_status)
qqnorm(residuals(model_test_male)); qqline(residuals(model_test_male))
summary(model_test_male)
kruskal.test(anova_df_male$tspec ~ anova_df_male$modif_status)

model_test_virf = aov(anova_df_virf$tspec ~ anova_df_virf$modif_status)
qqnorm(residuals(model_test_virf)); qqline(residuals(model_test_virf))
summary(model_test_virf)
kruskal.test(anova_df_virf$tspec ~ anova_df_virf$modif_status)

#table proportion
fb_female_m_p = fb_female_m / sum(fb_female_m)
fb_female_c_p = fb_female_c / sum(fb_female_c)
fb_male_m_p = fb_male_m / sum(fb_male_m)
fb_male_c_p = fb_male_c / sum(fb_male_c)
fb_virf_m_p = fb_virf_m / sum(fb_virf_m)
fb_virf_c_p = fb_virf_c / sum(fb_virf_c)

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

n_virf = sum(fb_virf_c +  fb_virf_m)
n_virf_c = sum(fb_virf_c)
n_virf_m = sum(fb_virf_m)
p_virf = n_virf_m / (n_virf_c + n_virf_m)
s_virf = fb_virf_c +  fb_virf_m
exp_virf_m = s_virf * p_virf
exp_virf_c = s_virf * (1 - p_virf)

chi_virf = sum((fb_virf_m - exp_virf_m)^2 / exp_virf_m) + sum((fb_virf_c - exp_virf_c)^2 / exp_virf_c)
df_virf = length(fb_virf_c) - 1
pchisq(chi_virf, df = df_virf, lower.tail = F)

#in female, modification occured differentially in function of the tissue, some tissues are more subject to domain modification
#not in male
#slignthy in virF

####2.chosen state for tissue spec: development time L3####
#sort dataset and filter out not considered states
flybase_expression_l3 = sort_RNASource_name_keep(dataset_flybase = flybase_expression, 
                                                 regexp_pattern = 'L3')

flybase_expression_l3$RNASource_name = as.character(flybase_expression_l3$RNASource_name)
unique(flybase_expression_l3$RNASource_name)
flybase_expression_l3$RNASource_name[flybase_expression_l3$RNASource_name == 'mE_mRNA_L3_CNS'] = 'mE_mRNA_L3_No_CNS'
flybase_expression_l3$RNASource_name[flybase_expression_l3$RNASource_name == 'mE_mRNA_L3_Wand_imag_disc'] = 'mE_mRNA_L3_Wand_imagdisc'
flybase_expression_l3$RNASource_name[flybase_expression_l3$RNASource_name == 'mE_mRNA_L3_Wand_dig_sys'] = 'mE_mRNA_L3_Wand_digsys'
#

#data organization
flybase_expression_l3 = flybase_expression_organization_l3(fly_exp_cond = flybase_expression_l3)
#

#Tspec calculation
flybase_expression_l3 = Tspec_inference(flybase_expression_l3)
#

#Inference of ubiquitous and specific gene and determine which it's the most expressed tissue for specific genes
flybase_expression_l3 = specificity_inference(flybase_expression_data = flybase_expression_l3,
                                                     threshold_specificity = 0.8)
#

###analysis
#dataset modification and control
flybase_expression_l3_modified = flybase_expression_l3[flybase_expression_l3$FBgn %in% modified_gene,]
flybase_expression_l3_control = flybase_expression_l3[!(flybase_expression_l3$FBgn %in% modified_gene),]

hist(flybase_expression_l3_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in L3 fly', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_l3_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

fb_l3_c = table(flybase_expression_l3_control$specificity)
fb_l3_m_tmp = table(flybase_expression_l3_modified$specificity)
fb_l3_m = fb_l3_c
for (list_tmp in 1:dim(fb_l3_m)){
  if (names(fb_l3_m[list_tmp]) %in% names(fb_l3_m_tmp)){
    fb_l3_m[list_tmp] = fb_l3_m_tmp[names(fb_l3_m[list_tmp])]
  }else{
    fb_l3_m[list_tmp] = 0
  }
}

col = palette()
considered_table = fb_l3_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in control L3 genes", col = col)

considered_table = fb_l3_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in modified L3 genes", col = col)

#test anova application
anova_df_l3 = rbind(flybase_expression_l3_control, flybase_expression_l3_modified)
anova_df_l3$modif_status = c(rep(0, length(flybase_expression_l3_control$tspec)),
                                 rep(1, length(flybase_expression_l3_modified$tspec)))
anova_df_l3$modif_status = as.factor(anova_df_l3$modif_status)

model_test_l3 = aov(anova_df_l3$tspec ~ anova_df_l3$modif_status)
qqnorm(residuals(model_test_l3)); qqline(residuals(model_test_l3))
summary(model_test_l3)
kruskal.test(anova_df_l3$tspec ~ anova_df_l3$modif_status)

#table proportion
fb_l3_c_p = fb_l3_c / sum(fb_l3_c)
fb_l3_m_p = fb_l3_m / sum(fb_l3_m)

#chi-square test: comparison between control and modification 
n_l3 = sum(fb_l3_c +  fb_l3_m)
n_l3_c = sum(fb_l3_c)
n_l3_m = sum(fb_l3_m)
p_l3 = n_l3_m / (n_l3_c + n_l3_m)
s_l3 = fb_l3_c +  fb_l3_m
exp_l3_m = s_l3 * p_l3
exp_l3_c = s_l3 * (1 - p_l3)

chi_l3 = sum((fb_l3_m - exp_l3_m)^2 / exp_l3_m) + sum((fb_l3_c - exp_l3_c)^2 / exp_l3_c)
df_l3 = length(fb_l3_c) - 1
pchisq(chi_l3, df = df_l3, lower.tail = F)

#no difference in tissue proportion between modification groups

####3.chosen state for tissue spec: development time A (1d, 4d, 20d) in carcass and dig_sys####
#sort dataset and filter out not considered states
flybase_expression_A = sort_RNASource_name_keep(dataset_flybase = flybase_expression, 
                                                 regexp_pattern = '_A_')
flybase_expression_A = sort_RNASource_name_notkeep(dataset_flybase = flybase_expression_A, 
                                                regexp_pattern = 'Mate')
flybase_expression_A = sort_RNASource_name_notkeep(dataset_flybase = flybase_expression_A, 
                                                regexp_pattern = 'VirF')


flybase_expression_A$RNASource_name = as.character(flybase_expression_A$RNASource_name)
unique(flybase_expression_A$RNASource_name)
flybase_expression_A$RNASource_name[flybase_expression_A$RNASource_name == 'mE_mRNA_A_1d_dig_sys'] = 'mE_mRNA_A_1d_digsys'
flybase_expression_A$RNASource_name[flybase_expression_A$RNASource_name == 'mE_mRNA_A_4d_dig_sys'] = 'mE_mRNA_A_4d_digsys'
flybase_expression_A$RNASource_name[flybase_expression_A$RNASource_name == 'mE_mRNA_A_20d_dig_sys'] = 'mE_mRNA_A_20d_digsys'
#

#data organization
flybase_expression_A_carcass = flybase_expression_organization_A_bytissue(fly_exp_cond = flybase_expression_A, considered_tissue = 'carcass')
flybase_expression_A_digsys = flybase_expression_organization_A_bytissue(fly_exp_cond = flybase_expression_A, considered_tissue = 'digsys')
#

#Tspec calculation
flybase_expression_A_carcass = Tspec_inference(flybase_expression_A_carcass)
flybase_expression_A_digsys = Tspec_inference(flybase_expression_A_digsys)
#

#Inference of ubiquitous and specific gene and determine which it's the most expressed tissue for specific genes
flybase_expression_A_carcass = specificity_inference(flybase_expression_data = flybase_expression_A_carcass,
                                              threshold_specificity = 0.8)
flybase_expression_A_digsys = specificity_inference(flybase_expression_data = flybase_expression_A_digsys,
                                                     threshold_specificity = 0.8)
#

###analysis
#dataset modification and control
flybase_expression_A_carcass_modified = flybase_expression_A_carcass[flybase_expression_A_carcass$FBgn %in% modified_gene,]
flybase_expression_A_carcass_control = flybase_expression_A_carcass[!(flybase_expression_A_carcass$FBgn %in% modified_gene),]
flybase_expression_A_digsys_modified = flybase_expression_A_digsys[flybase_expression_A_digsys$FBgn %in% modified_gene,]
flybase_expression_A_digsys_control = flybase_expression_A_digsys[!(flybase_expression_A_digsys$FBgn %in% modified_gene),]

hist(flybase_expression_A_carcass_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in carcass tissue in fly', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_A_carcass_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(flybase_expression_A_digsys_modified$tspec, breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of Tspec in digsys tissue in fly', xlab = 'Tspec value', cex.main = 0.9)
hist(flybase_expression_A_digsys_control$tspec, breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('top', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

fb_A_c_c = table(flybase_expression_A_carcass_control$specificity)
fb_A_c_m_tmp = table(flybase_expression_A_carcass_modified$specificity)
fb_A_d_c = table(flybase_expression_A_digsys_control$specificity)
fb_A_d_m_tmp = table(flybase_expression_A_digsys_modified$specificity)
fb_A_c_m = fb_A_c_c 
for (list_tmp in 1:dim(fb_A_c_m)){
  if (names(fb_A_c_m[list_tmp]) %in% names(fb_A_c_m_tmp)){
    fb_A_c_m[list_tmp] = fb_A_c_m_tmp[names(fb_A_c_m[list_tmp])]
  }else{
    fb_A_c_m[list_tmp] = 0
  }
}
fb_A_d_m = fb_A_d_c
for (list_tmp in 1:dim(fb_A_d_m)){
  if (names(fb_A_d_m[list_tmp]) %in% names(fb_A_d_m_tmp)){
    fb_A_d_m[list_tmp] = fb_A_d_m_tmp[names(fb_A_d_m[list_tmp])]
  }else{
    fb_A_d_m[list_tmp] = 0
  }
}

col = palette()
considered_table = fb_A_c_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in carcass control genes", col = col)

considered_table = fb_A_c_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in carcass modified genes", col = col)

col = palette()
considered_table = fb_A_d_c
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in digsys control genes", col = col)

considered_table = fb_A_d_m
lbls = paste(names(considered_table), sep= "")
pie(considered_table, labels = lbls, cex = 0.5, 
    main = "Specificity in digsys modified genes", col = col)

#test anova application
anova_df_carcass = rbind(flybase_expression_A_carcass_control, flybase_expression_A_carcass_modified)
anova_df_carcass$modif_status = c(rep(0, length(flybase_expression_A_carcass_control$tspec)),
                             rep(1, length(flybase_expression_A_carcass_modified$tspec)))
anova_df_carcass$modif_status = as.factor(anova_df_carcass$modif_status)

model_test_carcass = aov(anova_df_carcass$tspec ~ anova_df_carcass$modif_status)
qqnorm(residuals(model_test_carcass)); qqline(residuals(model_test_carcass))
summary(model_test_carcass)
kruskal.test(anova_df_carcass$tspec ~ anova_df_carcass$modif_status)

anova_df_digsys = rbind(flybase_expression_A_digsys_control, flybase_expression_A_digsys_modified)
anova_df_digsys$modif_status = c(rep(0, length(flybase_expression_A_digsys_control$tspec)),
                                  rep(1, length(flybase_expression_A_digsys_modified$tspec)))
anova_df_digsys$modif_status = as.factor(anova_df_digsys$modif_status)

model_test_digsys = aov(anova_df_digsys$tspec ~ anova_df_digsys$modif_status)
qqnorm(residuals(model_test_digsys)); qqline(residuals(model_test_digsys))
summary(model_test_digsys)
kruskal.test(anova_df_digsys$tspec ~ anova_df_digsys$modif_status)

#table proportion
fb_A_c_c_p = fb_A_c_c / sum(fb_A_c_c)
fb_A_c_m_p = fb_A_c_m / sum(fb_A_c_m)
fb_A_d_c_p = fb_A_d_c / sum(fb_A_d_c)
fb_A_d_m_p = fb_A_d_m / sum(fb_A_d_m)

#chi-square test: comparison between control and modification 
n_A_c_c = sum(fb_A_c_c +  fb_A_c_m)
n_A_c_c = sum(fb_A_c_c)
n_A_c_m = sum(fb_A_c_m)
p_A_c = n_A_c_m / (n_A_c_c + n_A_c_m)
s_A_c = fb_A_c_c +  fb_A_c_m
exp_A_c_m = s_A_c * p_A_c
exp_A_c_c = s_A_c * (1 - p_A_c)

chi_A_c = sum((fb_A_c_m - exp_A_c_m)^2 / exp_A_c_m) + sum((fb_A_c_c - exp_A_c_c)^2 / exp_A_c_c)
df_A_c = length(fb_A_c_c) - 1
pchisq(chi_A_c, df = df_A_c, lower.tail = F)

n_A_d_c = sum(fb_A_d_c +  fb_A_d_m)
n_A_d_c = sum(fb_A_d_c)
n_A_d_m = sum(fb_A_d_m)
p_A_d = n_A_d_m / (n_A_d_c + n_A_d_m)
s_A_d = fb_A_d_c +  fb_A_d_m
exp_A_d_m = s_A_d * p_A_d
exp_A_d_c = s_A_d * (1 - p_A_d)

chi_A_d = sum((fb_A_d_m - exp_A_d_m)^2 / exp_A_d_m) + sum((fb_A_d_c - exp_A_d_c)^2 / exp_A_d_c)
df_A_d = length(fb_A_d_c) - 1
pchisq(chi_A_d, df = df_A_d, lower.tail = F)

#no difference in time point proportion between modification groups