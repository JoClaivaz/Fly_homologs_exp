# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:59:55 2017

@author: Claivaz
Extraction of gene expression data in function of the presence of the gene identifier in the pair file
"""
def data_expression_extraction_final_pair(folder_specie_exp, specie_1, specie_2, final_pair_file, data_expression_file):
    
    import os
    
    #Extract the considered gene identifiers
    final_pair = open(final_pair_file, 'r')
    gene_ID_sp1 = []
    gene_ID_sp2 = []
    for pair_line in final_pair:
        if pair_line.split('\t')[0] not in gene_ID_sp1 and specie_1 in pair_line.split('\t')[4]:
            gene_ID_sp1.append(pair_line.split('\t')[0])
        elif pair_line.split('\t')[0] not in gene_ID_sp2 and specie_2 in pair_line.split('\t')[4]:
            gene_ID_sp2.append(pair_line.split('\t')[0])
        
        if pair_line.split('\t')[2] not in gene_ID_sp1 and specie_1 in pair_line.split('\t')[5]:
            gene_ID_sp1.append(pair_line.split('\t')[2])
        elif pair_line.split('\t')[2] not in gene_ID_sp2 and specie_2 in pair_line.split('\t')[5]:
            gene_ID_sp2.append(pair_line.split('\t')[2])
   
    final_pair.close()
    #
    
    #Extract information if data present in gene_ID
    expression_files_sp1 = os.listdir(folder_specie_exp + '/' + specie_1)
    expression_files_sp2 = os.listdir(folder_specie_exp + '/' + specie_2)
    
    data_expression = open(data_expression_file, 'w')
    data_expression.write('#specie\texperiment\tgeneID_fbgn\tanatomical_entity_name\tstage_name\tsex\tread\n')
    
    for expression_file in expression_files_sp1:
        file_tmp = open(folder_specie_exp + '/' + specie_1 + '/' + expression_file, 'r')
        for file_line in file_tmp:
            if file_line.split('\t')[3] in gene_ID_sp1 :
                data_expression.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(specie_1,
                                                                       file_line.split('\t')[0],
                                                                       file_line.split('\t')[3],
                                                                       file_line.split('\t')[5],
                                                                       file_line.split('\t')[7],
                                                                       file_line.split('\t')[8],
                                                                       file_line.split('\t')[10]))
        file_tmp.close()
    
    for expression_file in expression_files_sp2:
        file_tmp = open(folder_specie_exp + '/' + specie_2 + '/' + expression_file, 'r')
        for file_line in file_tmp:
            if file_line.split('\t')[3] in gene_ID_sp2 :
                data_expression.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(specie_2,
                                                                       file_line.split('\t')[0],
                                                                       file_line.split('\t')[3],
                                                                       file_line.split('\t')[5],
                                                                       file_line.split('\t')[7],
                                                                       file_line.split('\t')[8],
                                                                       file_line.split('\t')[10]))
        file_tmp.close()
    data_expression.close()

##DROMEvsDROAN
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROAN',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROAN_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROAN',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROAN_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN_nomodif')

##DROMEvsDROMO
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROMO',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROMO_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROMO')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROMO',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROMO_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROMO_nomodif')

##DROMEvsDROPS
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROPS',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROPS_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROPS')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROPS',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROPS_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROPS_nomodif')

##DROMEvsDROSI
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROSI',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROSI_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROSI')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROSI',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROSI_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROSI_nomodif')

##DROMEvsDROVI
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROVI',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROVI_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROVI')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROVI',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROVI_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROVI_nomodif')

##DROMEvsDROYA
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROYA',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROYA_domain_loss_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROYA')
#nomodif/control
data_expression_extraction_final_pair(folder_specie_exp = 'D:/UNIL/Master/Master_Project/Data/Bgee',
                                      specie_1 = 'DROME',
                                      specie_2 = 'DROYA',
                                      final_pair_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROYA_domain_nomodif_fbgn',
                                      data_expression_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROYA_nomodif')
              