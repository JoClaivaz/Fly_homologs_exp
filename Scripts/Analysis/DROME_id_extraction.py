# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170724
Extraction of unique DROME gene ID (DROMEXXXXX) present in the different datasets.
Then conversion from DROMEXXXX in FBgnXXXXX 
"""

import os
path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/'
expression_file = os.listdir(path_file)
DROME_identifier = []

for considered_file in expression_file:
    if '_QL' in considered_file:
        open_file = open(path_file + considered_file, 'r')
        #skip the headers
        open_file.readline()
        
        for line_file in open_file:
            considered_pair = line_file.split('\t')[1].replace('"','').split('_')
            if 'DROME' in considered_pair[0] and considered_pair[0] not in DROME_identifier:
                DROME_identifier.append(considered_pair[0])
            elif 'DROME' in considered_pair[1] and considered_pair[1] not in DROME_identifier:
                DROME_identifier.append(considered_pair[1])
        open_file.close()
        
#START New step to ensure to have all gene
nomodif_droan_gene = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROAN_domain_nomodif_fbgn', 'r')
modif_droan_gene = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROAN_domain_loss_fbgn', 'r')

for considered_pair in nomodif_droan_gene:
    if 'DROME' in considered_pair.split('\t')[4]:
        if considered_pair.split('\t')[4] not in DROME_identifier:
            DROME_identifier.append(considered_pair.split('\t')[4])
            
for considered_pair in nomodif_droan_gene:
    if 'DROME' in considered_pair.split('\t')[5]:
        if considered_pair.split('\t')[5] not in DROME_identifier:
            DROME_identifier.append(considered_pair.split('\t')[5])


#END new step

DROME_name = open('D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROME.fa', 'r')

FBgn_identifier = []
DROME_identifier_ordered = []

for DROME_line in DROME_name:
    if '>' in DROME_line:
        DROME_gene = DROME_line.split('|')[0].replace('>', '').replace(' ', '')
        if DROME_gene in DROME_identifier:
            DROME_identifier_ordered.append(DROME_gene)
            DROME_split = DROME_line.replace(' ', '').replace('\n', '').split('|')
            for DROME_slice in DROME_split:
                if 'FBgn' in DROME_slice:
                    FBgn_identifier.append(DROME_slice)
                    break
DROME_name.close()

DROME_FBgn_names = open(path_file + 'DROME_FBgn_names', 'w')

for index_tmp in range(len(DROME_identifier_ordered)):
    DROME_FBgn_names.write(DROME_identifier_ordered[index_tmp] + '\t' + 
                           FBgn_identifier[index_tmp] + '\n')
DROME_FBgn_names.close()

DROME_only_FBgn_names = open(path_file + 'DROME_only_FBgn_names', 'w')
DROME_only_FBgn_names.write('\n'.join(FBgn_identifier))
DROME_only_FBgn_names.close()