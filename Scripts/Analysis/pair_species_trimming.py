# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170727
Trimming pair gene present in each species pair dataset as no modification or simple modification
If complex modification in one specie, gene is removed
"""

#FUN
def pair_specie_to_dict(x, path_file):
    
    dict_fly = {}
    modif_file = open(path_file + 'final_pair_DROME_%s_domain_loss_fbgn' % (x), 'r')
    nomodif_file =  open(path_file + 'ortholog_DROME_%s_domain_nomodif_fbgn' % (x), 'r')
    
    for nomodif_line in nomodif_file:
        if 'DROME' in nomodif_line.replace('\n', '').split('\t')[4]:
            dict_fly[nomodif_line.replace('\n', '').split('\t')[0]] = 'no_modif'
        elif 'DROME' in nomodif_line.replace('\n', '').split('\t')[5]:
            dict_fly[nomodif_line.replace('\n', '').split('\t')[2]] = 'no_modif'
    
    for modif_line in modif_file:
        if 'DROME' in modif_line.replace('\n', '').split('\t')[4]:
            dict_fly[modif_line.replace('\n', '').split('\t')[0]] = 'modif'
        elif 'DROME' in modif_line.replace('\n', '').split('\t')[5]:
            dict_fly[modif_line.replace('\n', '').split('\t')[2]] = 'modif'
    
    modif_file.close()
    nomodif_file.close()
    
    return dict_fly



#Store in dictionary the considered pair gene for each species
DROAN_dict = pair_specie_to_dict(x = 'DROAN', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

DROMO_dict = pair_specie_to_dict(x = 'DROMO', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

DROPS_dict = pair_specie_to_dict(x = 'DROPS', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

DROSI_dict = pair_specie_to_dict(x = 'DROSI', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

DROVI_dict = pair_specie_to_dict(x = 'DROVI', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

DROYA_dict = pair_specie_to_dict(x = 'DROYA', 
                                 path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/')

#Filter genes: keep ones present in each dictionary
considered_genes = []

for tested_gene in DROAN_dict.keys():
    try:
        DROMO_dict[tested_gene]
        DROPS_dict[tested_gene]
        DROSI_dict[tested_gene]
        DROVI_dict[tested_gene]
        DROYA_dict[tested_gene]
        considered_genes.append(tested_gene)
    except:
        pass

#Write output file for analysis
output_file = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/all_species_genes', 'w')

output_file.write('DROME_geneID\tDROAN_status\tDROMO_status\tDROPS_status\tDROSI_status\tDROVI_status\tDROYA_status\n')

for tested_gene in considered_genes:
    output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (tested_gene, DROAN_dict[tested_gene], DROMO_dict[tested_gene],
                                                        DROPS_dict[tested_gene], DROSI_dict[tested_gene],
                                                        DROVI_dict[tested_gene], DROYA_dict[tested_gene]))
output_file.close()