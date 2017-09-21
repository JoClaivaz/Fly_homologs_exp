'''
Joaquim Claivaz
170907

Tspec analysis on modified and control group of genes in multispecies analysis
'''
#library needed
library(tidyr)

###Fun
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


#Need of common_dataset_completecase generated in gain_loss_DROME.R => load it
load(file = 'C:/Users/Claivaz/Desktop/gain_loss_domain_DROME_data')
str(common_dataset_completecase)

#store gene list with at least one modification in pairwise comparison
modified_gene = unique(common_dataset_completecase$FB_genes[common_dataset_completecase$status == 'modif'])

#load expression dataset from flybase of DROME
flybase_expression = read.table('D:/UNIL/Master/Master_Project/Data/flybase/gene_rpkm_report_fb_2017_04_formated',
                                sep = '\t', header = TRUE)

#filter out not considered genes
flybase_expression = flybase_expression[flybase_expression$FBgn %in% common_dataset_completecase$FB_genes,]
str(flybase_expression)

#keep only parent library name: modENCODE_mRNA-Seq_tissues
flybase_expression = flybase_expression[flybase_expression$Parent_library_name == 'modENCODE_mRNA-Seq_tissues',]

#see different state available (29 on 124) and filter out
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