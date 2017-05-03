###26 april 2017, Joaquim Claivaz

#This script allows the dataframe creation containing Tspecificity index, 
#the group of specificity, the length and the position of the loss between short and long paralog.

#Inputs: RData files from the function data_extraction.R, threshold for ubiquitous and tissue specific creation, and the destination folder
#Output: dataframe in format text containing: ID_short, ID_long, length_longdom, modif_dl, length_dl, tspec_short, tspec_long, max_exp_tissues_short, max_exp_tissues_long, specificity_short, specificity_long and specificity_shift

#Rscript "D:/UNIL/Master/Master_Project/Script/Current/DataframeCreation_TSpec_GroupShift.R --rd D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Human/GSE30352.RData --od D:/UNIL/Master/Master_Project/Script/try_170426 --tu 0.8 --ts 0.8


####Option for terminal running with different arguments####

#!/usr/bin/Rscript
# module add R/3.4.0

arg <- commandArgs(trailingOnly=TRUE)
library(optparse)

option_list <- list(                   
  make_option("--rd",  type="character",   default="D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Human/GSE30352.RData",help=("Complete path of RData file from data_extraction script (%default)"),metavar="RData_file"),
  make_option("--od",  type="character",   default="D:/UNIL/Master/Master_Project/Script/try_170426",help=("Complete path of the output files (%default)"),metavar="output_directory"),
  make_option("--tu", type="double", default=0.8, help=("Tspec values under this threshold will be considered as ubiquitous gene (%default)"), metavar="threshold_ubiquity"),
  make_option("--ts", type="double", default=0.8, help=("Tspec values above this threshold will be considered as tissue sepcific gene  (%default)"), metavar="threshold_specificity")
)

parser_object <- OptionParser(usage = "Usage: %prog [Options]", option_list=option_list, description="...")
opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE), positional_arguments=TRUE)

RData_file  = opt$options$rd
output_directory  = opt$options$rd
threshold_ubiquity = opt$options$tu
threshold_specificity = opt$options$ts

######FUNCTIONs######
####extract expression data from'data_extraction.R' output to new dataframe####
expression_data_to_dataframe = function(raw_data, col_rpkm, col_gene){
  data_frame = data.frame(raw_data[[1]][col_gene])
  for (tissues in names(raw_data)){
    data_frame = cbind(data_frame, raw_data[[tissues]][col_rpkm])
  }
  colnames(data_frame) = c(names(raw_data[[1]])[col_gene], names(raw_data))
  return(data_frame)
}

####log transformation of all RPKM present in dataframe from expression_data_to dataframe and normalization of all RPKM<1 equal to 0####
log_transformation_rpkm = function(data_frame){
  for (i in 2:dim(data_frame)[2]){
    data_frame[,i] = log10(data_frame[,i] + 0.000001)
    data_frame[,i][data_frame[,i] < 1] = 0
  }
  return(data_frame)
}

####tau calculation from Nadeza works####
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

####Tissue specificity index calculations####
tspec_calculation = function(data_frame){
  tau = data.frame(geneID = data_frame[,1], tspec = rep(NA, length(data_frame[,1])))
  for (i in 1:dim(data_frame)[1]){
    tmp_rpkm = NULL
    for (j in 2:dim(data_frame)[2]){
      tmp_rpkm = c(tmp_rpkm, data_frame[i,j])
    }
    tau[i,2] = fTau(tmp_rpkm)
  }
  return(tau)
}

####extraction of number of domain involved in modification, the position of the modification and the length of the longest paralog####
extraction_modificationdomain_lengthlongparalog = function(raw_data, col_length_paralog, col_modification){
  data_frame = data.frame(length_long_para = raw_data[[1]][col_length_paralog], num_domain_modif = rep(NA, length(raw_data[[1]][col_length_paralog])), position_modif = rep(NA,length(raw_data[[1]][col_length_paralog])))
  for(i in 1:dim(data_frame)[1]){
    if (grepl("\\|", raw_data[[1]][col_modification][i,])!=T){
      if (length(strsplit(as.character(raw_data[[1]][col_modification][i,]), "-")[[1]]) > 2){
        data_frame[i,3] = "f&b"
        data_frame[i,2] = 2
      }else{
        data_frame[i,3] = strsplit(as.character(raw_data[[1]][col_modification][i,]), "-")[[1]][1]
        data_frame[i,2] = strsplit(as.character(raw_data[[1]][col_modification][i,]), "-")[[1]][2]
      }
    }
  }
  return(data_frame)
}

####determination of tissue with maximal expression####
tissue_maximal_expression = function(raw_data){
  max_expression_tissue = NULL
  for (i in 1:dim(raw_data)[1]){
    tissue_exp_tmp = -1
    tissue = NA
    for (j in 2:dim(raw_data)[2]){
      if (all(is.na(raw_data[i,j]) == F, tissue_exp_tmp < raw_data[i,j])){
        tissue_exp_tmp = raw_data[i,j]
        tissue = colnames(raw_data)[j]
      }
    }
    max_expression_tissue[i] = tissue
  }
  return(max_expression_tissue)
}

########START running script########
##Input parameters
#RData_file="D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Human/GSE30352.RData"
#output_directory  = "D:/UNIL/Master/Master_Project/Script/try_170426"
#threshold_specificity=0.8
#threshold_ubiquity=0.8
##
load(RData_file)

domain_modification_length = extraction_modificationdomain_lengthlongparalog(tissues.list, 4, 5)

GeneExpression_short = expression_data_to_dataframe(tissues.list, 8, 6)
GeneExpression_long = expression_data_to_dataframe(tissues.list, 7, 1)

GeneExpression_short_transformed = log_transformation_rpkm(GeneExpression_short)
GeneExpression_long_transformed = log_transformation_rpkm(GeneExpression_long)

tspec_short_tmp = tspec_calculation(GeneExpression_short_transformed[!duplicated(GeneExpression_short_transformed),])
tspec_long_tmp = tspec_calculation(GeneExpression_long_transformed[!duplicated(GeneExpression_long_transformed),])

maximal_tissue_expression_short_tmp = tissue_maximal_expression(GeneExpression_short_transformed[!duplicated(GeneExpression_short_transformed),])
maximal_tissue_expression_long_tmp = tissue_maximal_expression(GeneExpression_long_transformed[!duplicated(GeneExpression_long_transformed),])

data_frame = data.frame(Ensembl_gene_ID_short=GeneExpression_short_transformed[,1], Ensembl_gene_ID_long=GeneExpression_long_transformed[,1], num_domain_long=domain_modification_length[,1], 
                      num_domain_modif=domain_modification_length[,2], position_modif=domain_modification_length[,3],
                      tspec_short = rep(NA, dim(GeneExpression_short_transformed)[1]), tspec_long = rep(NA, dim(GeneExpression_short_transformed)[1]),
                      maximal_tissue_expression_short = rep(NA, dim(GeneExpression_short_transformed)[1]), maximal_tissue_expression_long = rep(NA, dim(GeneExpression_short_transformed)[1]))


#insert in dataframe specific long and short tspec, maximal_tissue_expression according to the gene ID
for (i in 1:dim(tspec_short_tmp)[1]){
  col_index_short = grep(tspec_short_tmp[i,1], data_frame[,1])
  for (j in col_index_short){
    data_frame$tspec_short[j] = tspec_short_tmp[i,2]
    data_frame$maximal_tissue_expression_short[j] = maximal_tissue_expression_short_tmp[i]
  }
}

for (i in 1:dim(tspec_long_tmp)[1]){
  col_index_long = grep(tspec_long_tmp[i,1], data_frame[,2])
  for (j in col_index_long){
    data_frame$tspec_long[j] = tspec_long_tmp[i,2]
    data_frame$maximal_tissue_expression_long[j] = maximal_tissue_expression_long_tmp[i]
  }
}

#remove raws containing NA
data_frame = data_frame[complete.cases(data_frame),]
data_frame = cbind(data_frame, specificity_short = rep(NA, dim(data_frame)[1]), specificity_long = rep(NA, dim(data_frame)[1]), specificity_shift = rep(NA, dim(data_frame)[1]))

#determine the specificity of each paralog genes in function of tissue specificity index
for (i in 1:dim(data_frame)[1]){
  if (data_frame[i,]$tspec_long >= threshold_specificity){
    data_frame$specificity_long[i] = 'specific'
  }else if (data_frame[i,]$tspec_long < threshold_ubiquity){
    data_frame$specificity_long[i] = 'ubiquitous'
  }else {
    data_frame$specificity_long[i] = 'undefined'
  }
  
  if (data_frame[i,]$tspec_short >= threshold_specificity){
    data_frame$specificity_short[i] = 'specific'
  }else if (data_frame[i,]$tspec_short<threshold_ubiquity){
    data_frame$specificity_short[i] = 'ubiquitous'
  }else {
    data_frame$specificity_short[i] = 'undefined'
  }
}

#determine the group of shift specificity
for (i in 1:dim(data_frame)[1]){
  if(data_frame$specificity_long[i] == 'specific'){
    if (data_frame$specificity_short[i] == 'specific'){
      if (data_frame[i,]$maximal_tissue_expression_short == data_frame[i,]$maximal_tissue_expression_long){
        data_frame$specificity_shift[i] = 'no_shift'
      }else{
        data_frame$specificity_shift[i] = 'specificity_shift'
      }
    }else if (data_frame$specificity_short[i] == 'ubiquitous'){
      data_frame$specificity_shift[i] = 'specific_to_ubiquitous'
    }else{
      data_frame$specificity_shift[i] = 'undefined'
    }
  }else if (data_frame$specificity_long[i] == 'ubiquitous'){
    if(data_frame$specificity_short[i] == 'ubiquitous'){
      data_frame$specificity_shift[i] = 'no_shift'
    }else if (data_frame$specificity_short[i] == 'specific'){
      data_frame$specificity_shift[i] = 'ubiquitous_to_specific'
    }else{
      data_frame$specificity_shift[i] = 'undefined'
    }
  }else{
    data_frame$specificity_shift[i] = 'undefined'
  }
}

#set tspec as numeric, and delete first column X
data_frame$tspec_short = as.numeric(data_frame$tspec_short)
data_frame$tspec_long = as.numeric(data_frame$tspec_long)

#write dataframe in text format in the ouput folder
write.csv(data_frame, file = paste0(output_directory, '/', strsplit(basename(RData_file), '.RData')[1], '_df_out.txt'))
########END running script########

###Test
# load('D:/UNIL/Master/Master_Project/Script/Mirjam/Data/Human/GSE30352.RData')
# 
# data_test = data.frame(tissues.list[['testis']]$domain_difference)
# for (tissue in names(tissues.list)){
#   data_test = cbind(data_test, tissues.list[[tissue]]$RPKM_long, tissues.list[[tissue]]$RPKM_short)
# }
# data_test = data_test[complete.cases(data_test),]
# data_test = data_test[grepl("\\|",data_test$tissues.list...testis....domain_difference) !=T ,]
# dim(data_frame)[1] == dim(data_test)[1]
###