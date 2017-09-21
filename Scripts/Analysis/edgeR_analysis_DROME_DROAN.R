'''
Joaquim Claivaz
170830

EdgeR analysis pairwise species
'''
####library####
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages('tidyr')
library(edgeR)
library(tidyr)
#

####DROMEvsDROAN####
#Data organization# START
modif_group = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN', sep = '\t')
control_group = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN_nomodif', sep = '\t')
colnames(modif_group) = c('specie', 'experiment', 'geneID_fbgn', 'anatomical_entity_name', 'stage_name', 'sex', 'read')
colnames(control_group) = c('specie', 'experiment', 'geneID_fbgn', 'anatomical_entity_name', 'stage_name', 'sex', 'read')

modif_group$status = 'modification'
control_group$status = 'no_modification'

#check the common level amongst the two datasets and reduce the dataframe
anat_common = intersect(modif_group$anatomical_entity_name[modif_group$specie == 'DROME'], modif_group$anatomical_entity_name[modif_group$specie == 'DROAN'])
modif_group_reduce = modif_group[modif_group$anatomical_entity_name == anat_common,]
control_group_reduce = control_group[control_group$anatomical_entity_name == anat_common,]

reduced_data = rbind(modif_group_reduce, control_group_reduce)

#group replicates
reduced_data = aggregate(`read` ~ specie + experiment + geneID_fbgn + anatomical_entity_name +
                               stage_name + sex + status, data = reduced_data, paste, collapse = ' ')


reduced_data = separate(data = reduced_data, col = read, sep = ' ', into = c('read1', 'read2'))

#unique identifier per conditions
reduced_data = unite(reduced_data, 'sex_specie_status', c(sex, specie, status))
reduced_data$experiment = NULL
reduced_data$anatomical_entity_name = NULL
reduced_data$stage_name = NULL

str(reduced_data)
reduced_data$read1 = as.integer(reduced_data$read1)
reduced_data$read2 = as.integer(reduced_data$read2)
str(reduced_data)
#Data organization# END

#group the different pair with modification
pair_modification = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROAN_domain_loss_fbgn', sep = '\t')

levels(as.factor(reduced_data$sex_specie_status))

read1_drome_female = c()
read2_drome_female = c()
read1_droan_female = c()
read2_droan_female = c()
read1_drome_male = c()
read2_drome_male = c()
read1_droan_male = c()
read2_droan_male = c()


for (line_pair in 1:dim(pair_modification)[1]){
  if (grepl('DROME', pair_modification[line_pair,5])){
    read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][1]
    read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][1]
    read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][2]
    read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][2]
  } else {
    read1_droan_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][1]
    read2_droan_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][1]
    read1_droan_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][2]
    read2_droan_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,1]][2]
  }
  if (grepl('DROME', pair_modification[line_pair,6])){
    read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][1]
    read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][1]
    read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][2]
    read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][2]
  } else {
    read1_droan_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][1]
    read2_droan_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][1]
    read1_droan_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][2]
    read2_droan_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_modification[line_pair,3]][2]
  }
}

pair_modif_name = unite(pair_modification, 'ortho_pair', c(V5, V6))
data_modification = data.frame(pair_modif_name$ortho_pair, read1_drome_female, read2_drome_female, read1_drome_male, read2_drome_male, read1_droan_female, read2_droan_female, read1_droan_male, read2_droan_male)

#group the different pair control
pair_control = read.table('D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROAN_domain_nomodif_fbgn', sep = '\t')

read1_drome_female = c()
read2_drome_female = c()
read1_droan_female = c()
read2_droan_female = c()
read1_drome_male = c()
read2_drome_male = c()
read1_droan_male = c()
read2_droan_male = c()


for (line_pair in 1:dim(pair_control)[1]){
  if (grepl('DROME', pair_control[line_pair,5])){
    read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][1]
    read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][1]
    read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][2]
    read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][2]
  } else {
    read1_droan_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][1]
    read2_droan_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][1]
    read1_droan_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][2]
    read2_droan_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,1]][2]
  }
  if (grepl('DROME', pair_control[line_pair,6])){
    read1_drome_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][1]
    read2_drome_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][1]
    read1_drome_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][2]
    read2_drome_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][2]
  } else {
    read1_droan_female[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][1]
    read2_droan_female[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][1]
    read1_droan_male[line_pair] = reduced_data$read1[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][2]
    read2_droan_male[line_pair] = reduced_data$read2[as.character(reduced_data$geneID_fbgn) == pair_control[line_pair,3]][2]
  }
}

pair_control_name = unite(pair_control, 'ortho_pair', c(V5, V6))
data_control = data.frame(pair_control_name$ortho_pair, read1_drome_female, read2_drome_female, read1_drome_male, read2_drome_male, read1_droan_female, read2_droan_female, read1_droan_male, read2_droan_male)

str(data_control)

names(data_modification)[1] = 'ortho_pair'
names(data_control)[1] = 'ortho_pair'

#in control some NA expression value, need to remove it for some species
data_control = data_control[complete.cases(data_control),]
#

data_modif_control = rbind(data_control, data_modification)
####edgeR####
######DGE QL DROMEvsDROAN / design with two factors: strain, sex (and replicates)######
species = c(rep('DROME', 4), rep('DROAN', 4))
sex = c(rep('female', 2), rep('male', 2), rep('female', 2), rep('male', 2))
replicates = c(rep(c('exp1', 'exp2'), 4))

edgeR_complete = DGEList(counts = data_modif_control[,2:9], group = rep(1:4, each = 2), genes = data_modif_control[,1])

#filter out not or low expressed gene
keep = rowSums(cpm(edgeR_complete) > 1) >= 2
edgeR_complete = edgeR_complete[keep, , keep.lib.sizes=FALSE]

#TMM normalization
edgeR_complete = calcNormFactors(edgeR_complete)

###Consider gene length difference###START###
#https://support.bioconductor.org/p/79379/
#https://www.biostars.org/p/99310/

#Extract gene length for each pair
pair_length_data = read.table(file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/DROME_DROAN_protein_length.txt', sep = '\t')

pair_length = edgeR_complete$genes
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
edgeR_complete$counts[,1] = edgeR_complete$counts[,1] / normalization_length
edgeR_complete$counts[,2] = edgeR_complete$counts[,2] / normalization_length
edgeR_complete$counts[,3] = edgeR_complete$counts[,3] / normalization_length
edgeR_complete$counts[,4] = edgeR_complete$counts[,4] / normalization_length
###END###

#dispersion estimation using glm (multi factors and non normal distribution)
group = rep(1:4, each = 2)
design = model.matrix(~0 + group, data = edgeR_complete$samples)

edgeR_complete = estimateDisp(edgeR_complete, design)

plotMDS(edgeR_complete, xlim = c(-5, 8), cex = 0.7, col = rep(palette(),each = 2))
plotBCV(edgeR_complete)

#if only modification group take into account
y = DGEList(counts = data_modification[,2:9], group = rep(1:4, each = 2), genes = data_modification[,1])
keep = rowSums(cpm(y) > 1) >= 2
y = y[keep, , keep.lib.sizes=FALSE]
y = calcNormFactors(y)
pair_length = y$genes
pair_length = separate(pair_length, genes, c('s1', 's2'))
s1_length = c()
s2_length = c()
for (pair in 1 : dim(pair_length)[1]){
  s1_length[pair] = pair_length_data$V2[pair_length_data$V1 == pair_length$s1[pair]]
  s2_length[pair] = pair_length_data$V2[pair_length_data$V1 == pair_length$s2[pair]]
}
pair_length = data.frame(pair_length, s1_length, s2_length)
normalization_length = pair_length$s1_length / pair_length$s2_length
y$counts[,1] = y$counts[,1] / normalization_length
y$counts[,2] = y$counts[,2] / normalization_length
y$counts[,3] = y$counts[,3] / normalization_length
y$counts[,4] = y$counts[,4] / normalization_length
y = estimateDisp(y, design)
plotBCV(y)
#

####QL (better for low number of replicates, for LRT see old script)
###DROMEvsDROAN
fit_complete_QL = glmQLFit(edgeR_complete, design)

considered_contrasts = makeContrasts(MvsF_sp1 = group2 - group1, MvsF_sp2 = group4 -group3,
                                     sp1VSsp2_M = group4 - group2, sp1VSsp2_F = group3 - group1, levels = design) 

#maleVSfemale DROME
QL_mVSf_sp1 = glmQLFTest(fit_complete_QL, contrast = considered_contrasts[,'MvsF_sp1']) 
topTags(QL_mVSf_sp1)
#maleVSfemale DROAN
QL_mVSf_sp2 = glmQLFTest(fit_complete_QL, contrast = considered_contrasts[,'MvsF_sp2']) 
topTags(QL_mVSf_sp2)
#DROMEvsDROAN male
QL_sp1VSsp2_male = glmQLFTest(fit_complete_QL, contrast = considered_contrasts[,'sp1VSsp2_M']) 
topTags(QL_sp1VSsp2_male)
#DROMEvsDROAN female
QL_sp1VSsp2_female = glmQLFTest(fit_complete_QL, contrast = considered_contrasts[,'sp1VSsp2_F']) 
topTags(QL_sp1VSsp2_female)

#store pair control and modif
data_pair_control = data_control$ortho_pair
data_pair_modif = data_modification$ortho_pair

####number of barplot_edgeR.R
#1st plot
length(data_pair_control)
length(data_pair_modif)
#2nd plot
length(QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif])
length(QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control])

##plot
boxplot(QL_sp1VSsp2_female$table$logFC[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif], 
        QL_sp1VSsp2_female$table$logFC[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control],
        QL_sp1VSsp2_male$table$logFC[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif],
        QL_sp1VSsp2_male$table$logFC[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control],
        main = 'Distribution of differential expression between DROME and DROAN (QL)', 
        ylab = 'LogFC', names = c('female\nmodification', 'female\ncontrol', 'male\nmodification', 'male\ncontrol'),
        cex.main = 0.8, cex.axis = 0.8)

####ANOVA or equivalent non parametric####
###QL data
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

logFC_drome_droan_QL = rbind( female_modif, male_modif, female_control, male_control)
str(logFC_drome_droan_QL)

#distribution
hist(logFC_drome_droan_QL$logFC[logFC_drome_droan_QL$status == 'modif'], breaks = 100, freq = F, col = rgb(1, 0 , 0, 0.5), main = 'Distribution of logFC between DROME and DROAN ortholog (QL)', xlab = 'logFC value', cex.main = 0.9)
hist(logFC_drome_droan_QL$logFC[logFC_drome_droan_QL$status == 'control'], breaks = 100, freq = F, col = rgb(0, 0 , 1, 0.5), add = T)
legend('right', c("Domain modification", "No modification"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), cex = 0.8, horiz = F)

boxplot(logFC ~ status * sex, data = logFC_drome_droan_QL)

#condition anova
model_QL = aov(logFC ~ status * sex, data = logFC_drome_droan_QL)

qqnorm(residuals(model_QL)); qqline(residuals(model_QL))
plot(logFC_drome_droan_QL$logFC ~ logFC_drome_droan_QL$status * logFC_drome_droan_QL$sex)
#no normal distribution, ANOVA cant be applied

kruskal.test(logFC ~ status, data = logFC_drome_droan_QL)
kruskal.test(logFC ~ sex, data = logFC_drome_droan_QL)

#in function of sex
kruskal.test(logFC ~ status, data = logFC_drome_droan_QL[logFC_drome_droan_QL$sex == 'male',])
kruskal.test(logFC ~ status, data = logFC_drome_droan_QL[logFC_drome_droan_QL$sex == 'female',])

###Write table
write.table(logFC_drome_droan_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/logFC_drome_droan_QL', sep = '\t')

###p-value datas
female_modif_k = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif],
                          pval = QL_sp1VSsp2_female$table$PValue[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_modif], 
                          status = 'modif', sex = 'female') 
male_modif_k = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif],
                          pval = QL_sp1VSsp2_male$table$PValue[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif], 
                        status = 'modif', sex = 'male') 
female_control_k = data.frame(genes = QL_sp1VSsp2_female$genes[[1]][QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control],
                              pval = QL_sp1VSsp2_female$table$PValue[QL_sp1VSsp2_female$genes[[1]] %in% data_pair_control], 
                            status = 'control', sex = 'female') 
male_control_k = data.frame(genes = QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control],
                            pval = QL_sp1VSsp2_male$table$PValue[QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control], 
                          status = 'control', sex = 'male')

pvalue_drome_droan_QL = rbind( female_modif_k, male_modif_k, female_control_k, male_control_k)
write.table(pvalue_drome_droan_QL, 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/pvalue/pvalue_drome_droan_QL', sep = '\t')
####END p-value data

####Gene onthology####
###Load gene ID table converter
gene_ID = read.table(file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/DROME_id_converter')

###install right database
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Dm.eg.db")

###modification
genes_list = separate(data = data.frame(genes = as.character(QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_modif])),
                      col = 'genes',into = c('sp1', 'sp2'), sep = '_') 
genes_list_drome_modif = c()
for (pair in 1:dim(genes_list)[1]){
  if (grepl('DROME', genes_list$sp1[pair])){
    genes_list_drome_modif[pair] = as.character(gene_ID$V2[gene_ID$V1 == genes_list$sp1[pair]])
  }else{
    genes_list_drome_modif[pair] = as.character(gene_ID$V2[gene_ID$V1 == genes_list$sp2[pair]])
  }
}

###control
genes_list = separate(data = data.frame(genes = as.character(QL_sp1VSsp2_male$genes[[1]][QL_sp1VSsp2_male$genes[[1]] %in% data_pair_control])),
                      col = 'genes',into = c('sp1', 'sp2'), sep = '_') 
genes_list_drome_control = c()
for (pair in 1:dim(genes_list)[1]){
  if (grepl('DROME', genes_list$sp1[pair])){
    genes_list_drome_control[pair] = as.character(gene_ID$V2[gene_ID$V1 == genes_list$sp1[pair]])
  }else{
    genes_list_drome_control[pair] = as.character(gene_ID$V2[gene_ID$V1 == genes_list$sp2[pair]])
  }
}


###Gene ontology enrichment modif group
##package installation
#source("https://bioconductor.org/biocLite.R")
#biocLite("topGO")
#biocLite(c("graph", "RBGL", "Rgraphviz"))
library(topGO)
library(graph)
library(RBGL)
library(Rgraphviz)

all_genes = c(genes_list_drome_control, genes_list_drome_modif)
all_genes = as.vector(all_genes)
names(all_genes) = c(genes_list_drome_control, genes_list_drome_modif)
all_genes = c(rep(0, length(genes_list_drome_control)),rep(1, length(genes_list_drome_modif)))
names(all_genes) = c(genes_list_drome_control, genes_list_drome_modif)
genes_modif = all_genes[all_genes == 1]

all_genes = as.factor(all_genes)
genes_modif = as.factor(genes_modif)

#biological process
topGO_droan_bp = new("topGOdata", 
                   ontology = "BP", 
                   allGenes = all_genes, 
                   geneSel = genes_modif,
                   nodeSize = 10, 
                   annot = annFUN.org,
                   mapping="org.Dm.eg.db", 
                   ID = "ensembl")
#summary
topGO_droan_bp

#enrichment test
resultFisher_droan_bp = runTest(topGO_droan_bp, algorithm = "classic", statistic = "fisher")
resultFisher_droan_bp

resultKS_droan_bp = runTest(topGO_droan_bp, algorithm = "classic", statistic = "ks")
resultKS.elim_droan_bp = runTest(topGO_droan_bp, algorithm = "elim", statistic = "ks")

#analysis
allRes_droan_bp = GenTable(topGO_droan_bp, classicFisher = resultFisher_droan_bp, 
                   classicKS = resultKS_droan_bp, elimKS = resultKS.elim_droan_bp,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)

allRes_droan_bp
write.csv(allRes_droan_bp, file = 'C:/Users/Claivaz/Desktop/allRes_droan')

pValue.classic_droan = score(resultKS_droan_bp)
pValue.elim_droan = score(resultKS.elim_droan_bp)[names(pValue.classic_droan)]
gstat_droan = termStat(topGO_droan_bp, names(pValue.classic_droan))
gSize_droan = gstat_droan$Annotated / max(gstat_droan$Annotated) * 4

#Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_droan = colMap(gstat_droan$Significant)
plot(pValue.classic_droan, pValue.elim_droan, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_droan, col = gCol_droan)

sel.go_droan = names(pValue.classic_droan)[pValue.elim_droan < pValue.classic_droan]
cbind(termStat(topGO_droan_bp, sel.go_droan),
      elim = pValue.elim_droan[sel.go_droan],
      classic = pValue.classic_droan[sel.go_droan])

showSigOfNodes(topGO_droan_bp, score(resultKS.elim_droan_bp), firstSigNodes = 5, useInfo = 'all')

#Molecular function
topGO_droan_mf = new("topGOdata", 
                     ontology = "MF", 
                     allGenes = all_genes, 
                     geneSel = genes_modif,
                     nodeSize = 10, 
                     annot = annFUN.org,
                     mapping="org.Dm.eg.db", 
                     ID = "ensembl")
#summary
topGO_droan_mf

#enrichment test
resultFisher_droan_mf = runTest(topGO_droan_mf, algorithm = "classic", statistic = "fisher")
resultFisher_droan_mf

resultKS_droan_mf = runTest(topGO_droan_mf, algorithm = "classic", statistic = "ks")
resultKS.elim_droan_mf = runTest(topGO_droan_mf, algorithm = "elim", statistic = "ks")

#analysis
allRes_droan_mf = GenTable(topGO_droan_mf, classicFisher = resultFisher_droan_mf, 
                           classicKS = resultKS_droan_mf, elimKS = resultKS.elim_droan_mf,
                           orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)

allRes_droan_mf
write.csv(allRes_droan_mf, file = 'C:/Users/Claivaz/Desktop/allRes_droan')
