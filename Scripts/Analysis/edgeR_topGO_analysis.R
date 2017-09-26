'''
Joaquim Claivaz
170925

EdgeR and topGO analysis, pairwise species comparisons
'''
####library needed####
##edgeR part
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages('tidyr')
library(edgeR)
library(tidyr)
##topGO part
#biocLite(c("topGO", "org.Dm.eg.db", "graph", "RBGL", "Rgraphviz"))
library(topGO)
library(graph)
library(RBGL)
library(Rgraphviz)
#

####Functions####
data_organization = function(modif_group_file, control_group_file){
  
  #store Upper case name of the second species
  specie2_UC = strsplit(modif_group_file, 'data_expression_DROME_')[[1]][2]
  
  #load expression datasets and rename the columns
  modif_group = read.table(modif_group_file, sep = '\t')
  control_group = read.table(control_group_file, sep = '\t')
  colnames(modif_group) = c('specie', 'experiment', 'geneID_fbgn', 'anatomical_entity_name', 'stage_name', 'sex', 'read')
  colnames(control_group) = c('specie', 'experiment', 'geneID_fbgn', 'anatomical_entity_name', 'stage_name', 'sex', 'read')
  
  #add status factor
  modif_group$status = 'modification'
  control_group$status = 'no_modification'
  
  #check the common level/state amongst the two datasets and reduce the dataframe
  stage_common = intersect(modif_group$stage_name[modif_group$specie == 'DROME'],
                           modif_group$stage_name[modif_group$specie == specie2_UC])
  modif_group_reduce = modif_group[modif_group$stage_name == stage_common,]
  control_group_reduce = control_group[control_group$stage_name == stage_common,]
  
  anat_common = intersect(modif_group_reduce$anatomical_entity_name[modif_group_reduce$specie == 'DROME'],
                          modif_group_reduce$anatomical_entity_name[modif_group_reduce$specie == specie2_UC])
  modif_group_reduce = modif_group_reduce[modif_group_reduce$anatomical_entity_name == anat_common,]
  control_group_reduce = control_group_reduce[control_group_reduce$anatomical_entity_name == anat_common,]
  
  reduced_data =(rbind(modif_group_reduce, control_group_reduce))
  
  #group replicates
  reduced_data = aggregate(`read` ~ specie + experiment + geneID_fbgn + anatomical_entity_name +
                             stage_name + sex + status, data = reduced_data, paste, collapse = ' ')
  
  #test the separate function, some species have 4 considered libraries while others have 2
  options(warn = 2)
  separate_try = try(separate(data = reduced_data, col = read, sep = ' ', into = c('read1', 'read2')))
  options(warn = 1)
  if(class(separate_try) == 'try-error'){
    reduced_data = separate(data = reduced_data, col = read, sep = ' ', into = c('read1', 'read2', 'read3', 'read4'))
    reduced_data$read3 = as.integer(reduced_data$read3)
    reduced_data$read4 = as.integer(reduced_data$read4)
  }else{
    reduced_data = separate(data = reduced_data, col = read, sep = ' ', into = c('read1', 'read2'))
  }
  
  #unique identifier per conditions
  reduced_data = unite(reduced_data, 'sex_specie_status', c(sex, specie, status))
  reduced_data$experiment = NULL
  reduced_data$anatomical_entity_name = NULL
  reduced_data$stage_name = NULL
  
  reduced_data$read1 = as.integer(reduced_data$read1)
  reduced_data$read2 = as.integer(reduced_data$read2)
  
  return(reduced_data)
}

expression_data_considering_ortholog_4read = function(reduced_data, ortholog_pair_file){
  
  ortholog_pair = read.table(ortholog_pair_file, sep = '\t')
  
  read1_drome_female = c()
  read2_drome_female = c()
  read1_specie2_female = c()
  read2_specie2_female = c()
  read1_drome_male = c()
  read2_drome_male = c()
  read1_specie2_male = c()
  read2_specie2_male = c()
  
  for (line_pair in 1:dim(ortholog_pair)[1]){
    if (grepl('DROME', ortholog_pair[line_pair,5])){
      read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
    } else {
      read1_specie2_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read2_specie2_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read1_specie2_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read2_specie2_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
    }
    if (grepl('DROME', ortholog_pair[line_pair,6])){
      read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
    } else {
      read1_specie2_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read2_specie2_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read1_specie2_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read2_specie2_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
    }
  }
  
  pair_name = unite(ortholog_pair, 'ortho_pair', c(V5, V6))
  data_output = data.frame(pair_name$ortho_pair, read1_drome_female, read2_drome_female, read1_drome_male, read2_drome_male,
                           read1_specie2_female, read2_specie2_female, read1_specie2_male, read2_specie2_male)
  names(data_output)[1] = 'ortho_pair'
  
  return(data_output)
}

expression_data_considering_ortholog_6read = function(reduced_data, ortholog_pair_file){
  
  ortholog_pair = read.table(ortholog_pair_file, sep = '\t')
  
  read1_drome_female = c()
  read2_drome_female = c()
  read1_specie2_female = c()
  read2_specie2_female = c()
  read3_specie2_female = c()
  read4_specie2_female = c()
  read1_drome_male = c()
  read2_drome_male = c()
  read1_specie2_male = c()
  read2_specie2_male = c()
  read3_specie2_male = c()
  read4_specie2_male = c()
  
  for (line_pair in 1:dim(ortholog_pair)[1]){
    if (grepl('DROME', ortholog_pair[line_pair,5])){
      read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
    } else {
      read1_specie2_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read2_specie2_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read3_specie2_female[line_pair] = reduced_data$read3[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read4_specie2_female[line_pair] = reduced_data$read4[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][1]
      read1_specie2_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read2_specie2_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read3_specie2_male[line_pair] = reduced_data$read3[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
      read4_specie2_male[line_pair] = reduced_data$read4[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,1]][2]
    }
    if (grepl('DROME', ortholog_pair[line_pair,6])){
      read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
    } else {
      read1_specie2_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read2_specie2_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read3_specie2_female[line_pair] = reduced_data$read3[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read4_specie2_female[line_pair] = reduced_data$read4[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][1]
      read1_specie2_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read2_specie2_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read3_specie2_male[line_pair] = reduced_data$read3[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
      read4_specie2_male[line_pair] = reduced_data$read4[as.character(reduced_data$geneID_fbgn) == ortholog_pair[line_pair,3]][2]
    }
  }
  
  pair_name = unite(ortholog_pair, 'ortho_pair', c(V5, V6))
  data_output = data.frame(pair_name$ortho_pair, read1_drome_female, read2_drome_female, read1_drome_male, read2_drome_male,
                           read1_specie2_female, read2_specie2_female, read3_specie2_female, read4_specie2_female, 
                           read1_specie2_male, read2_specie2_male, read3_specie2_male, read4_specie2_male)
  names(data_output)[1] = 'ortho_pair'
  
  return(data_output)
}

expression_data_to_dataframe = function(reduced_data, ortholog_pair_file_modif, ortholog_pair_file_control){
  if(dim(reduced_data)[2] == 4){
    data_modification = expression_data_considering_ortholog_4read(reduced_data, ortholog_pair_file_modif)
    data_control = expression_data_considering_ortholog_4read(reduced_data, ortholog_pair_file_control)
  }else if(dim(reduced_data)[2] == 6){
    data_modification = expression_data_considering_ortholog_6read(reduced_data, ortholog_pair_file_modif)
    data_control = expression_data_considering_ortholog_6read(reduced_data, ortholog_pair_file_control)
  }
  
  #in control some NA expression value, need to remove it for some species
  data_control = data_control[complete.cases(data_control),]
  #
  
  data_modif_control = rbind(data_control, data_modification)
  
  return(list(data_modif_control, data_control$ortho_pair, data_modification$ortho_pair))
}

gene_length_normalization = function(pair_length_file, edgeR_dataset){
  pair_length_data = read.table(file = pair_length_file, sep = '\t')
  pair_length = edgeR_dataset$genes
  pair_length = separate(pair_length, genes, c('s1', 's2'))
  s1_length = c()
  s2_length = c()
  for (pair in 1 : dim(pair_length)[1]){
    s1_length[pair] = pair_length_data$V2[pair_length_data$V1 == pair_length$s1[pair]]
    s2_length[pair] = pair_length_data$V2[pair_length_data$V1 == pair_length$s2[pair]]
  }
  pair_length = data.frame(pair_length, s1_length, s2_length)
  normalization_length = pair_length$s1_length / pair_length$s2_length
  
  #normalization of count in function of length of ortholog pair (specie1$count/(specie1Length/specie2Length))
  edgeR_dataset$counts[,1] = edgeR_dataset$counts[,1] / normalization_length
  edgeR_dataset$counts[,2] = edgeR_dataset$counts[,2] / normalization_length
  edgeR_dataset$counts[,3] = edgeR_dataset$counts[,3] / normalization_length
  edgeR_dataset$counts[,4] = edgeR_dataset$counts[,4] / normalization_length
  
  return(edgeR_dataset)
}

dispersion_estimation = function(edgeR_dataset){
  if(dim(edgeR_dataset$counts)[2] == 8){
    group = rep(1:4, each = 2)
    design = model.matrix(~0 + group, data = edgeR_dataset$samples)
  }else if(dim(edgeR_dataset$counts)[2] == 12){
    group = c(rep(1:2, each = 2), rep(3:4, each = 4))
    design = model.matrix(~0 + group, data = edgeR_dataset$samples)
  }
  edgeR_dataset = estimateDisp(edgeR_dataset, design)
  
  return(edgeR_dataset)
}

QL_fitting = function(edgeR_dataset){
  if(dim(edgeR_dataset$counts)[2] == 8){
    group = rep(1:4, each = 2)
    design = model.matrix(~0 + group, data = edgeR_dataset$samples)
    considered_contrasts = makeContrasts(MvsF_sp1 = group2 - group1, MvsF_sp2 = group4 -group3,
                                         sp1VSsp2_M = group4 - group2, sp1VSsp2_F = group3 - group1, levels = design)
  }else if(dim(edgeR_dataset$counts)[2] == 12){
    group = c(rep(1:2, each = 2), rep(3:4, each = 4))
    design = model.matrix(~0 + group, data = edgeR_dataset$samples)
    considered_contrasts = makeContrasts(MvsF_sp1 = group2 - group1, MvsF_sp2 = group4 -group3,
                                         sp1VSsp2_M = group4 - group2, sp1VSsp2_F = group3 - group1, levels = design)
  }
  fit_QL = glmQLFit(edgeR_dataset, design)
  
  return(list(fit_QL, considered_contrasts))
}

logFC_into_dataframe = function(QL_sp1VSsp2_female, QL_sp1VSsp2_male, data_pair_modif, data_pair_control){
  female_modif = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif],
                            logFC = QL_sp1VSsp2_female$table$logFC[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif], 
                            status = 'modif', sex = 'female') 
  male_modif = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif],
                          logFC = QL_sp1VSsp2_male$table$logFC[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif], 
                          status = 'modif', sex = 'male') 
  female_control = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control],
                              logFC = QL_sp1VSsp2_female$table$logFC[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control], 
                              status = 'control', sex = 'female') 
  male_control = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control],
                            logFC = QL_sp1VSsp2_male$table$logFC[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control], 
                            status = 'control', sex = 'male')
  
  logFC_sp1VSsp2_QL = rbind( female_modif, male_modif, female_control, male_control)
  
  return(logFC_sp1VSsp2_QL)
}

pval_into_dataframe = function(QL_sp1VSsp2_female, QL_sp1VSsp2_male, data_pair_modif, data_pair_control){
  female_modif = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif],
                            logFC = QL_sp1VSsp2_female$table$PValue[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif], 
                            status = 'modif', sex = 'female') 
  male_modif = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif],
                          logFC = QL_sp1VSsp2_male$table$PValue[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif], 
                          status = 'modif', sex = 'male') 
  female_control = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control],
                              logFC = QL_sp1VSsp2_female$table$PValue[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control], 
                              status = 'control', sex = 'female') 
  male_control = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control],
                            logFC = QL_sp1VSsp2_male$table$PValue[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control], 
                            status = 'control', sex = 'male')
  
  pval_sp1VSsp2_QL = rbind( female_modif, male_modif, female_control, male_control)
  
  return(pval_sp1VSsp2_QL)
}

topGO_data_organization = function(gene_list, considered_gene, gene_ID){
  gene_OMA = separate(data = data.frame(genes = as.character(gene_list$genes[[1]][gene_list$genes[[1]] %in% considered_gene])),
                      col = 'genes',into = c('sp1', 'sp2'), sep = '_')
  gene_fbgn = c()
  for (pair in 1:dim(gene_OMA)[1]){
    if (grepl('DROME', gene_OMA$sp1[pair])){
      gene_fbgn[pair] = as.character(gene_ID$V2[gene_ID$V1 == gene_OMA$sp1[pair]])
    }else{
      gene_fbgn[pair] = as.character(gene_ID$V2[gene_ID$V1 == gene_OMA$sp2[pair]])
    }
  }
  return(gene_fbgn)
}

topGO_analysis = function(genes_list_modif, genes_list_control, GO_enrichment){
  all_genes = c(rep(0, length(genes_list_control)),rep(1, length(genes_list_modif)))
  names(all_genes) = c(genes_list_control, genes_list_modif)
  genes_modif = all_genes[all_genes == 1]
  all_genes = as.factor(all_genes)
  genes_modif = as.factor(genes_modif)
  
  topGO_GO_enrichment = new("topGOdata", 
                            ontology = GO_enrichment, 
                            allGenes = all_genes, 
                            geneSel = genes_modif,
                            nodeSize = 10, 
                            annot = annFUN.org,
                            mapping="org.Dm.eg.db", 
                            ID = "ensembl")
  
  resultFisher = runTest(topGO_GO_enrichment, algorithm = "classic", statistic = "fisher")
  resultKS = runTest(topGO_GO_enrichment, algorithm = "classic", statistic = "ks")
  resultKS.elim = runTest(topGO_GO_enrichment, algorithm = "elim", statistic = "ks")
  
  allRes = GenTable(topGO_GO_enrichment, classicFisher = resultFisher, 
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)
  
  return(list(topGO_GO_enrichment, allRes, resultKS.elim))
}
#

####Data Organization####
reduced_data_droan = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN_nomodif')
reduced_data_dromo = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROMO',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROMO_nomodif')
reduced_data_drops = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROPS',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROPS_nomodif')
reduced_data_drosi = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROSI',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROSI_nomodif')
reduced_data_drovi = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROVI',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROVI_nomodif')
reduced_data_droya = data_organization(modif_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROYA',
                                       control_group_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROYA_nomodif')
#

####Expression data into ortholog's dataframe and store control/modification pair genes####
list_droan = expression_data_to_dataframe(reduced_data = reduced_data_droan,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROAN_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROAN_domain_nomodif_fbgn')
data_droan = list_droan[[1]]
gene_modif_droan = list_droan[[3]]
gene_control_droan = list_droan[[2]]

list_dromo = expression_data_to_dataframe(reduced_data = reduced_data_dromo,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROMO_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROMO_domain_nomodif_fbgn')
data_dromo = list_dromo[[1]]
gene_modif_dromo = list_dromo[[3]]
gene_control_dromo = list_dromo[[2]]

list_drops = expression_data_to_dataframe(reduced_data = reduced_data_drops,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROPS_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROPS_domain_nomodif_fbgn')
data_drops = list_drops[[1]]
gene_modif_drops = list_drops[[3]]
gene_control_drops = list_drops[[2]]

list_drosi = expression_data_to_dataframe(reduced_data = reduced_data_drosi,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROSI_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROSI_domain_nomodif_fbgn')
data_drosi = list_drosi[[1]]
gene_modif_drosi = list_drosi[[3]]
gene_control_drosi = list_drosi[[2]]

list_drovi = expression_data_to_dataframe(reduced_data = reduced_data_drovi,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROVI_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROVI_domain_nomodif_fbgn')
data_drovi = list_drovi[[1]]
gene_modif_drovi = list_drovi[[3]]
gene_control_drovi = list_drovi[[2]]

list_droya = expression_data_to_dataframe(reduced_data = reduced_data_droya,
                                          ortholog_pair_file_modif = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROYA_domain_loss_fbgn', 
                                          ortholog_pair_file_control = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROYA_domain_nomodif_fbgn')
data_droya = list_droya[[1]]
gene_modif_droya = list_droya[[3]]
gene_control_droya = list_droya[[2]]
#

####edgeR analysis####
edgeR_droan = DGEList(counts = data_droan[,2:9], group = rep(1:4, each = 2), genes = data_droan[,1])
edgeR_dromo = DGEList(counts = data_dromo[,2:9], group = rep(1:4, each = 2), genes = data_dromo[,1])
edgeR_drops = DGEList(counts = data_drops[,2:9], group = rep(1:4, each = 2), genes = data_drops[,1])
edgeR_drovi = DGEList(counts = data_drovi[,2:9], group = rep(1:4, each = 2), genes = data_drovi[,1])
edgeR_droya = DGEList(counts = data_droya[,2:9], group = rep(1:4, each = 2), genes = data_droya[,1])
edgeR_drosi = DGEList(counts = data_drosi[,2:13], group = c(rep(1:2, each = 2),rep(3:4, each = 4)), genes = data_drosi[,1])

#filter out not or low expressed gene
keep = rowSums(cpm(edgeR_droan) > 1) >= 2
edgeR_droan = edgeR_droan[keep, , keep.lib.sizes=FALSE]
keep = rowSums(cpm(edgeR_dromo) > 1) >= 2
edgeR_dromo = edgeR_dromo[keep, , keep.lib.sizes=FALSE]
keep = rowSums(cpm(edgeR_drops) > 1) >= 2
edgeR_drops = edgeR_drops[keep, , keep.lib.sizes=FALSE]
keep = rowSums(cpm(edgeR_drosi) > 1) >= 2
edgeR_drosi = edgeR_drosi[keep, , keep.lib.sizes=FALSE]
keep = rowSums(cpm(edgeR_drovi) > 1) >= 2
edgeR_drovi = edgeR_drovi[keep, , keep.lib.sizes=FALSE]
keep = rowSums(cpm(edgeR_droya) > 1) >= 2
edgeR_droya = edgeR_droya[keep, , keep.lib.sizes=FALSE]
#

#TMM normalization
edgeR_droan = calcNormFactors(edgeR_droan)
edgeR_dromo = calcNormFactors(edgeR_dromo)
edgeR_drops = calcNormFactors(edgeR_drops)
edgeR_drosi = calcNormFactors(edgeR_drosi)
edgeR_drovi = calcNormFactors(edgeR_drovi)
edgeR_droya = calcNormFactors(edgeR_droya)
#

#Gene length normalization
edgeR_droan = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROAN_protein_length.txt',
                                        edgeR_droan)
edgeR_dromo = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROMO_protein_length.txt',
                                        edgeR_dromo)
edgeR_drops = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROPS_protein_length.txt',
                                        edgeR_drops)
edgeR_drosi = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROSI_protein_length.txt',
                                        edgeR_drosi)
edgeR_drovi = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROVI_protein_length.txt',
                                        edgeR_drovi)
edgeR_droya = gene_length_normalization(pair_length_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROYA_protein_length.txt',
                                        edgeR_droya)
#

#dispersion estimation using glm (multi factors and non normal distribution) and plots
edgeR_droan = dispersion_estimation(edgeR_dataset = edgeR_droan)
plotMDS(edgeR_droan, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_droan)

edgeR_dromo = dispersion_estimation(edgeR_dataset = edgeR_dromo)
plotMDS(edgeR_dromo, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_dromo)

edgeR_drops = dispersion_estimation(edgeR_dataset = edgeR_drops)
plotMDS(edgeR_drops, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_drops)

edgeR_drosi = dispersion_estimation(edgeR_dataset = edgeR_drosi)
plotMDS(edgeR_drosi, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_drosi)

edgeR_drovi = dispersion_estimation(edgeR_dataset = edgeR_drovi)
plotMDS(edgeR_drovi, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_drovi)

edgeR_droya = dispersion_estimation(edgeR_dataset = edgeR_droya)
plotMDS(edgeR_droya, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_droya)
#

#QL fitting (better for low number of replicates than LRT)
list_droan = QL_fitting(edgeR_dataset = edgeR_droan)
droan_QL_fit = list_droan[[1]]
droan_contrasts = list_droan[[2]]

list_dromo = QL_fitting(edgeR_dataset = edgeR_dromo)
dromo_QL_fit = list_dromo[[1]]
dromo_contrasts = list_dromo[[2]]

list_drops = QL_fitting(edgeR_dataset = edgeR_drops)
drops_QL_fit = list_drops[[1]]
drops_contrasts = list_drops[[2]]

list_drosi = QL_fitting(edgeR_dataset = edgeR_drosi)
drosi_QL_fit = list_drosi[[1]]
drosi_contrasts = list_drosi[[2]]

list_drovi = QL_fitting(edgeR_dataset = edgeR_drovi)
drovi_QL_fit = list_drovi[[1]]
drovi_contrasts = list_drovi[[2]]

list_droya = QL_fitting(edgeR_dataset = edgeR_droya)
droya_QL_fit = list_droya[[1]]
droya_contrasts = list_droya[[2]]
#

#test DGE significance and print most significant DGE
droan_male_QL = glmQLFTest(droan_QL_fit, contrast = droan_contrasts[,'sp1VSsp2_M']) 
topTags(droan_male_QL)
droan_female_QL = glmQLFTest(droan_QL_fit, contrast = droan_contrasts[,'sp1VSsp2_F']) 
topTags(droan_female_QL)

dromo_male_QL = glmQLFTest(dromo_QL_fit, contrast = dromo_contrasts[,'sp1VSsp2_M']) 
topTags(dromo_male_QL)
dromo_female_QL = glmQLFTest(dromo_QL_fit, contrast = dromo_contrasts[,'sp1VSsp2_F']) 
topTags(dromo_female_QL)

drops_male_QL = glmQLFTest(drops_QL_fit, contrast = drops_contrasts[,'sp1VSsp2_M']) 
topTags(drops_male_QL)
drops_female_QL = glmQLFTest(drops_QL_fit, contrast = drops_contrasts[,'sp1VSsp2_F']) 
topTags(drops_female_QL)

drosi_male_QL = glmQLFTest(drosi_QL_fit, contrast = drosi_contrasts[,'sp1VSsp2_M']) 
topTags(drosi_male_QL)
drosi_female_QL = glmQLFTest(drosi_QL_fit, contrast = drosi_contrasts[,'sp1VSsp2_F']) 
topTags(drosi_female_QL)

drovi_male_QL = glmQLFTest(drovi_QL_fit, contrast = drovi_contrasts[,'sp1VSsp2_M']) 
topTags(drovi_male_QL)
drovi_female_QL = glmQLFTest(drovi_QL_fit, contrast = drovi_contrasts[,'sp1VSsp2_F']) 
topTags(drovi_female_QL)

droya_male_QL = glmQLFTest(droya_QL_fit, contrast = droya_contrasts[,'sp1VSsp2_M']) 
topTags(droya_male_QL)
droya_female_QL = glmQLFTest(droya_QL_fit, contrast = droya_contrasts[,'sp1VSsp2_F']) 
topTags(droya_female_QL)
#

#store logFC into dataframe and save dataset
logFC_drome_droan_QL = logFC_into_dataframe(QL_sp1VSsp2_female = droan_female_QL, 
                                            QL_sp1VSsp2_male = droan_male_QL,
                                            data_pair_modif = gene_modif_droan,
                                            data_pair_control = gene_control_droan)
write.table(logFC_drome_droan_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droan_QL', sep = '\t')
logFC_drome_dromo_QL = logFC_into_dataframe(QL_sp1VSsp2_female = dromo_female_QL, 
                                            QL_sp1VSsp2_male = dromo_male_QL,
                                            data_pair_modif = gene_modif_dromo,
                                            data_pair_control = gene_control_dromo)
write.table(logFC_drome_dromo_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_dromo_QL', sep = '\t')
logFC_drome_drops_QL = logFC_into_dataframe(QL_sp1VSsp2_female = drops_female_QL, 
                                            QL_sp1VSsp2_male = drops_male_QL,
                                            data_pair_modif = gene_modif_drops,
                                            data_pair_control = gene_control_drops)
write.table(logFC_drome_drops_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drops_QL', sep = '\t')
logFC_drome_drosi_QL = logFC_into_dataframe(QL_sp1VSsp2_female = drosi_female_QL, 
                                            QL_sp1VSsp2_male = drosi_male_QL,
                                            data_pair_modif = gene_modif_drosi,
                                            data_pair_control = gene_control_drosi)
write.table(logFC_drome_drosi_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drosi_QL', sep = '\t')
logFC_drome_drovi_QL = logFC_into_dataframe(QL_sp1VSsp2_female = drovi_female_QL, 
                                            QL_sp1VSsp2_male = drovi_male_QL,
                                            data_pair_modif = gene_modif_drovi,
                                            data_pair_control = gene_control_drovi)
write.table(logFC_drome_drovi_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_drovi_QL', sep = '\t')
logFC_drome_droya_QL = logFC_into_dataframe(QL_sp1VSsp2_female = droya_female_QL, 
                                            QL_sp1VSsp2_male = droya_male_QL,
                                            data_pair_modif = gene_modif_droya,
                                            data_pair_control = gene_control_droya)
write.table(logFC_drome_droya_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droya_QL', sep = '\t')
#

#store pvalue into dataframe and save dataset
pval_drome_droan_QL = pval_into_dataframe(QL_sp1VSsp2_female = droan_female_QL, 
                                            QL_sp1VSsp2_male = droan_male_QL,
                                            data_pair_modif = gene_modif_droan,
                                            data_pair_control = gene_control_droan)
write.table(pval_drome_droan_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_droan_QL', sep = '\t')
pval_drome_dromo_QL = pval_into_dataframe(QL_sp1VSsp2_female = dromo_female_QL, 
                                            QL_sp1VSsp2_male = dromo_male_QL,
                                            data_pair_modif = gene_modif_dromo,
                                            data_pair_control = gene_control_dromo)
write.table(pval_drome_dromo_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_dromo_QL', sep = '\t')
pval_drome_drops_QL = pval_into_dataframe(QL_sp1VSsp2_female = drops_female_QL, 
                                            QL_sp1VSsp2_male = drops_male_QL,
                                            data_pair_modif = gene_modif_drops,
                                            data_pair_control = gene_control_drops)
write.table(pval_drome_drops_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_drops_QL', sep = '\t')
pval_drome_drosi_QL = pval_into_dataframe(QL_sp1VSsp2_female = drosi_female_QL, 
                                            QL_sp1VSsp2_male = drosi_male_QL,
                                            data_pair_modif = gene_modif_drosi,
                                            data_pair_control = gene_control_drosi)
write.table(pval_drome_drosi_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_drosi_QL', sep = '\t')
pval_drome_drovi_QL = pval_into_dataframe(QL_sp1VSsp2_female = drovi_female_QL, 
                                            QL_sp1VSsp2_male = drovi_male_QL,
                                            data_pair_modif = gene_modif_drovi,
                                            data_pair_control = gene_control_drovi)
write.table(pval_drome_drovi_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_drovi_QL', sep = '\t')
pval_drome_droya_QL = pval_into_dataframe(QL_sp1VSsp2_female = droya_female_QL, 
                                            QL_sp1VSsp2_male = droya_male_QL,
                                            data_pair_modif = gene_modif_droya,
                                            data_pair_control = gene_control_droya)
write.table(pval_drome_droya_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_droya_QL', sep = '\t')
#

#boxplot of logFC in function of the sex and the status (control or domain modification) 
boxplot(droan_female_QL$table$logFC[droan_female_QL$genes[[1]] %in% gene_modif_droan], 
        droan_female_QL$table$logFC[droan_female_QL$genes[[1]] %in% gene_control_droan],
        droan_male_QL$table$logFC[droan_male_QL$genes[[1]] %in% gene_modif_droan],
        droan_male_QL$table$logFC[droan_male_QL$genes[[1]] %in% gene_control_droan],
        main = 'Distribution of differential expression between DROME and DROAN (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
boxplot(dromo_female_QL$table$logFC[dromo_female_QL$genes[[1]] %in% gene_modif_dromo], 
        dromo_female_QL$table$logFC[dromo_female_QL$genes[[1]] %in% gene_control_dromo],
        dromo_male_QL$table$logFC[dromo_male_QL$genes[[1]] %in% gene_modif_dromo],
        dromo_male_QL$table$logFC[dromo_male_QL$genes[[1]] %in% gene_control_dromo],
        main = 'Distribution of differential expression between DROME and DROMO (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
boxplot(drops_female_QL$table$logFC[drops_female_QL$genes[[1]] %in% gene_modif_drops], 
        drops_female_QL$table$logFC[drops_female_QL$genes[[1]] %in% gene_control_drops],
        drops_male_QL$table$logFC[drops_male_QL$genes[[1]] %in% gene_modif_drops],
        drops_male_QL$table$logFC[drops_male_QL$genes[[1]] %in% gene_control_drops],
        main = 'Distribution of differential expression between DROME and DROPS (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
boxplot(drosi_female_QL$table$logFC[drosi_female_QL$genes[[1]] %in% gene_modif_drosi], 
        drosi_female_QL$table$logFC[drosi_female_QL$genes[[1]] %in% gene_control_drosi],
        drosi_male_QL$table$logFC[drosi_male_QL$genes[[1]] %in% gene_modif_drosi],
        drosi_male_QL$table$logFC[drosi_male_QL$genes[[1]] %in% gene_control_drosi],
        main = 'Distribution of differential expression between DROME and DROSI (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
boxplot(drovi_female_QL$table$logFC[drovi_female_QL$genes[[1]] %in% gene_modif_drovi], 
        drovi_female_QL$table$logFC[drovi_female_QL$genes[[1]] %in% gene_control_drovi],
        drovi_male_QL$table$logFC[drovi_male_QL$genes[[1]] %in% gene_modif_drovi],
        drovi_male_QL$table$logFC[drovi_male_QL$genes[[1]] %in% gene_control_drovi],
        main = 'Distribution of differential expression between DROME and DROVI (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
boxplot(droya_female_QL$table$logFC[droya_female_QL$genes[[1]] %in% gene_modif_droya], 
        droya_female_QL$table$logFC[droya_female_QL$genes[[1]] %in% gene_control_droya],
        droya_male_QL$table$logFC[droya_male_QL$genes[[1]] %in% gene_modif_droya],
        droya_male_QL$table$logFC[droya_male_QL$genes[[1]] %in% gene_control_droya],
        main = 'Distribution of differential expression between DROME and DROYA (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)
#

#distribution of logFC (histogram)
hist(logFC_drome_droan_QL$logFC[logFC_drome_droan_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROAN ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_droan_QL$logFC[logFC_drome_droan_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(logFC_drome_dromo_QL$logFC[logFC_drome_dromo_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROMO ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_dromo_QL$logFC[logFC_drome_dromo_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(logFC_drome_drops_QL$logFC[logFC_drome_drops_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROPS ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_drops_QL$logFC[logFC_drome_drops_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(logFC_drome_drosi_QL$logFC[logFC_drome_drosi_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROSI ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_drosi_QL$logFC[logFC_drome_drosi_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(logFC_drome_drovi_QL$logFC[logFC_drome_drovi_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROVI ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_drovi_QL$logFC[logFC_drome_drovi_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

hist(logFC_drome_droya_QL$logFC[logFC_drome_droya_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROYA ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_droya_QL$logFC[logFC_drome_droya_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)
#

#ANOVA analysis: effect of status (domain modification or control) on logFC value, no interraction effects
droan_model_QL = aov(logFC ~ status + sex, data = logFC_drome_droan_QL)
summary(droan_model_QL)
dromo_model_QL = aov(logFC ~ status + sex, data = logFC_drome_dromo_QL)
summary(dromo_model_QL)
drops_model_QL = aov(logFC ~ status + sex, data = logFC_drome_drops_QL)
summary(drops_model_QL)
drosi_model_QL = aov(logFC ~ status + sex, data = logFC_drome_drosi_QL)
summary(drosi_model_QL)
drovi_model_QL = aov(logFC ~ status + sex, data = logFC_drome_drovi_QL)
summary(drovi_model_QL)
droya_model_QL = aov(logFC ~ status + sex, data = logFC_drome_droya_QL)
summary(droya_model_QL)
#

####topGO analysis####
#Load gene ID table converter (if not available run: 'DROME_id_extraction.py' and 'DROME_ID_conversion.py')
gene_ID = read.table(file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/DROME_id_converter')
#

#Data organization for topGO (correct fbgn identifier instead of OMA identifier)
droan_modif_topGO = topGO_data_organization(gene_list = droan_male_QL, 
                                            considered_gene = gene_modif_droan,
                                            gene_ID = gene_ID)
droan_control_topGO = topGO_data_organization(gene_list = droan_male_QL, 
                                            considered_gene = gene_control_droan,
                                            gene_ID = gene_ID)
dromo_modif_topGO = topGO_data_organization(gene_list = dromo_male_QL, 
                                            considered_gene = gene_modif_dromo,
                                            gene_ID = gene_ID)
dromo_control_topGO = topGO_data_organization(gene_list = dromo_male_QL, 
                                              considered_gene = gene_control_dromo,
                                              gene_ID = gene_ID)
drops_modif_topGO = topGO_data_organization(gene_list = drops_male_QL, 
                                            considered_gene = gene_modif_drops,
                                            gene_ID = gene_ID)
drops_control_topGO = topGO_data_organization(gene_list = drops_male_QL, 
                                              considered_gene = gene_control_drops,
                                              gene_ID = gene_ID)
drosi_modif_topGO = topGO_data_organization(gene_list = drosi_male_QL, 
                                            considered_gene = gene_modif_drosi,
                                            gene_ID = gene_ID)
drosi_control_topGO = topGO_data_organization(gene_list = drosi_male_QL, 
                                              considered_gene = gene_control_drosi,
                                              gene_ID = gene_ID)
drovi_modif_topGO = topGO_data_organization(gene_list = drovi_male_QL, 
                                            considered_gene = gene_modif_drovi,
                                            gene_ID = gene_ID)
drovi_control_topGO = topGO_data_organization(gene_list = drovi_male_QL, 
                                              considered_gene = gene_control_drovi,
                                              gene_ID = gene_ID)
droya_modif_topGO = topGO_data_organization(gene_list = droya_male_QL, 
                                            considered_gene = gene_modif_droya,
                                            gene_ID = gene_ID)
droya_control_topGO = topGO_data_organization(gene_list = droya_male_QL, 
                                              considered_gene = gene_control_droya,
                                              gene_ID = gene_ID)
#

#topGO analysis
droan_topGO_BP = topGO_analysis(genes_list_modif = droan_modif_topGO,
                                genes_list_control = droan_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(droan_topGO_BP[[1]], score(droan_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(droan_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/droan_topGO_BP')

droan_topGO_MF = topGO_analysis(genes_list_modif = droan_modif_topGO,
                                genes_list_control = droan_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(droan_topGO_MF[[1]], score(droan_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(droan_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/droan_topGO_MF')

dromo_topGO_BP = topGO_analysis(genes_list_modif = dromo_modif_topGO,
                                genes_list_control = dromo_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(dromo_topGO_BP[[1]], score(dromo_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(dromo_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/dromo_topGO_BP')

dromo_topGO_MF = topGO_analysis(genes_list_modif = dromo_modif_topGO,
                                genes_list_control = dromo_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(dromo_topGO_MF[[1]], score(dromo_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(dromo_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/dromo_topGO_MF')

drops_topGO_BP = topGO_analysis(genes_list_modif = drops_modif_topGO,
                                genes_list_control = drops_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(drops_topGO_BP[[1]], score(drops_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drops_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/drops_topGO_BP')

drops_topGO_MF = topGO_analysis(genes_list_modif = drops_modif_topGO,
                                genes_list_control = drops_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(drops_topGO_MF[[1]], score(drops_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drops_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/drops_topGO_MF')

drosi_topGO_BP = topGO_analysis(genes_list_modif = drosi_modif_topGO,
                                genes_list_control = drosi_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(drosi_topGO_BP[[1]], score(drosi_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drosi_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/drosi_topGO_BP')

drosi_topGO_MF = topGO_analysis(genes_list_modif = drosi_modif_topGO,
                                genes_list_control = drosi_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(drosi_topGO_MF[[1]], score(drosi_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drosi_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/drosi_topGO_MF')

drovi_topGO_BP = topGO_analysis(genes_list_modif = drovi_modif_topGO,
                                genes_list_control = drovi_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(drovi_topGO_BP[[1]], score(drovi_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drovi_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/drovi_topGO_BP')

drovi_topGO_MF = topGO_analysis(genes_list_modif = drovi_modif_topGO,
                                genes_list_control = drovi_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(drovi_topGO_MF[[1]], score(drovi_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(drovi_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/drovi_topGO_MF')

droya_topGO_BP = topGO_analysis(genes_list_modif = droya_modif_topGO,
                                genes_list_control = droya_control_topGO,
                                GO_enrichment = 'BP')
showSigOfNodes(droya_topGO_BP[[1]], score(droya_topGO_BP[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(droya_topGO_BP[[2]], file = 'C:/Users/Claivaz/Desktop/droya_topGO_BP')

droya_topGO_MF = topGO_analysis(genes_list_modif = droya_modif_topGO,
                                genes_list_control = droya_control_topGO,
                                GO_enrichment = 'MF')
showSigOfNodes(droya_topGO_MF[[1]], score(droya_topGO_MF[[3]]), firstSigNodes = 5, useInfo = 'all')
#write.csv(droya_topGO_MF[[2]], file = 'C:/Users/Claivaz/Desktop/droya_topGO_MF')
