# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:26:55 2017

@author: Claivaz
Extract Fbgn (flybase identifier) corresponding to OMA identifier and write again the 
final_pair_specie1_specie2_domain_loss_fbgn (new file with right fbgn name)
"""
def extract_fbgn_OMAname(final_pair_file_in, final_pair_file_out, ortholog_OMA):
    
    pair_in = open(final_pair_file_in, 'r')
    final_pair_file = open(final_pair_file_out, 'w')
    
    for pair_line in pair_in:
        OMA_file = open(ortholog_OMA, 'r')
        
        for OMA_line in OMA_file:
            if pair_line.replace('\t', ' ').split(' ')[0] in OMA_line:
                OMA_split = OMA_line.split(' | ')
                for OMA_string in OMA_split:
                    if 'FBgn' in OMA_string:
                        spec1_tmp_fbgn = OMA_string.replace('\n', '')
                        spec1_tmp_OMA = pair_line.replace('\t', ' ').split(' ')[0] 
                        break
            if pair_line.replace('\t', ' ').split(' ')[2] in OMA_line:
                OMA_split = OMA_line.split(' | ')
                for OMA_string in OMA_split:
                    if 'FBgn' in OMA_string:
                        spec2_tmp_fbgn = OMA_string.replace('\n', '')
                        spec2_tmp_OMA = pair_line.replace('\t', ' ').split(' ')[2]
                        break
                    
        final_pair_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (spec1_tmp_fbgn, pair_line.replace('\t', ' ').split(' ')[1].replace('\n', ''), 
                                                            spec2_tmp_fbgn, pair_line.replace('\t', ' ').split(' ')[3].replace('\n', ''),
                                                            spec1_tmp_OMA, spec2_tmp_OMA))
        
    OMA_file.close()
    pair_in.close()
    final_pair_file.close()
    
###DROME_DROAN
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/final_pair_DROME_DROAN_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROAN_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROAN')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/ortholog_DROME_DROAN_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROAN_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROAN')



###DROME_DROMO
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/final_pair_DROME_DROMO_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROMO_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROMO')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/ortholog_DROME_DROMO_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROMO_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROMO')




###DROME_DROPS
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/final_pair_DROME_DROPS_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROPS_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROPS')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/ortholog_DROME_DROPS_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROPS_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROPS')




###DROME_DROSI
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/final_pair_DROME_DROSI_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROSI_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROSI')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/ortholog_DROME_DROSI_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROSI_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROSI')




###DROME_DROVI
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/final_pair_DROME_DROVI_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROVI_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROVI')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/ortholog_DROME_DROVI_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROVI_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROVI')




###DROME_DROYA
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/final_pair_DROME_DROYA_domain_loss',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/final_pair_DROME_DROYA_domain_loss_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROYA')
#nomodif/control
extract_fbgn_OMAname(final_pair_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/ortholog_DROME_DROYA_domain_nomodif',
                     final_pair_file_out = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/ortholog_DROME_DROYA_domain_nomodif_fbgn',
                     ortholog_OMA = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROYA')