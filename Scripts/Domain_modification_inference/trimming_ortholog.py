# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:28:30 2017

@author: Claivaz
Script extraction of pair not considered for analysis: different of 1:1 ortholog
"""

def pair_no_1_1_ortholog(ortholog_oma_file, output_file):
    ortholog_oma = open('D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROAN', 'r')
    output_oma = open(output_file, 'w')
    
    for ortholog_pair in ortholog_oma:
        if '1:1' not in ortholog_pair:
            pair_split = ortholog_pair.split('\t')
            if 'DROME' in pair_split[0]:
                output_oma.write('%s\t%s\n' % (pair_split[0], pair_split[1]))
            else:
                output_oma.write('%s\t%s\n' % (pair_split[1], pair_split[0]))
    
    ortholog_oma.close()
    output_oma.close()
            
            

#DROMEvsDROAN
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROAN',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROAN_not_consider')

#DROMEvsDROMO
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROMO',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROMO_not_consider')

#DROMEvsDROPS
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROPS',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROPS_not_consider')

#DROMEvsDROSI
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROSI',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROSI_not_consider')

#DROMEvsDROVI
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROVI',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROVI_not_consider')

#DROMEvsDROYA
pair_no_1_1_ortholog(ortholog_oma_file = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROYA',
                     output_file = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROYA_not_consider')