# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170520
Extraction of pfamscan data, in function of pairwise_OMA
"""
def extract_all_ortholog_from_pfamscan(pairwise_oma_in, pfamscan_in, all_ortolog_out):
    
    pairwise_oma = open(pairwise_oma_in, 'r')
    domain_out = open(all_ortolog_out, 'w')
    
    pair_number = 0
    for pair in pairwise_oma:
        if ':' in pair.split('\t')[2] :
            pair_number += 1
            pfamscan_out = open(pfamscan_in, 'r')
            for domain_line in pfamscan_out:
                
                if pair.split('\t')[0] in domain_line:
                    domain_line = domain_line.replace(' ', '\t')
                    while '\t\t' in domain_line:
                        domain_line = domain_line.replace('\t\t', '\t')
                    domain_out.write(str(pair_number) + '\t' + domain_line.split('\t')[0] + '\t' +
                                     domain_line.split('\t')[5] + '\t' + domain_line.split('\t')[8] + '\t' +
                                     domain_line.split('\t')[9] + '\t' + domain_line.split('\t')[10] + '\n')
                
                elif pair.split('\t')[1] in domain_line:
                    domain_line = domain_line.replace(' ', '\t')
                    while '\t\t' in domain_line:
                        domain_line = domain_line.replace('\t\t', '\t')
                    domain_out.write(str(pair_number) + '\t' + domain_line.split('\t')[0] + '\t' +
                                     domain_line.split('\t')[5] + '\t' + domain_line.split('\t')[8] + '\t' +
                                     domain_line.split('\t')[9] + '\t' + domain_line.split('\t')[10] + '\n')
            pfamscan_out.close()
    pairwise_oma.close()
    domain_out.close()            
    
#Run the function / WINDOWS
#DROME vs DROYA
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROYA/pairwise_ortholog_DROME_DROYA', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROYA_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROYA_domain_parsed')

#DROME vs DROPS
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROPS/pairwise_ortholog_DROME_DROPS', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROPS_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROPS_domain_parsed')

#DROME vs DROAN
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROAN', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROAN_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROAN_domain_parsed')

#DROME vs DROMO
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROMO', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROMO_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROMO_domain_parsed')

#DROME vs DROSI
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROSI', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROSI_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROSI_domain_parsed')

#DROME vs DROVI
extract_all_ortholog_from_pfamscan(pairwise_oma_in = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROVI', 
                                      pfamscan_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/ortholog_DROME_DROVI_domain',
                                      all_ortolog_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROVI_domain_parsed')