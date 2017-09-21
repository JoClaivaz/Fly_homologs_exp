# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170522
Extraction of the unique ortholog pair, considering that the loss it's not part of repetition domains and inferred by domainDIFF
"""
def extract_ortholog_domain_loss_domainDiff (domain_loss_in, domain_modification_in, output_file):
    
    domain_loss = open(domain_loss_in, 'r')
    
    pair_out = open(output_file, 'w')
    
    
    for domain_loss_line in domain_loss:
        
        try:
            domain_loss_line.split(' ')[2]
            domain_modification = open(domain_modification_in, 'r')
            for domain_modification_line in domain_modification:
                
                if domain_loss_line.replace('\n', '').split(' ')[2] in domain_modification_line and domain_loss_line.split(' ')[1] in domain_modification_line:
                    pair_out.write(domain_modification_line)
                    domain_modification_in.close()
        except:
            pass
    
    domain_loss.close()
    pair_out.close()
        
#Run the function / WINDOWS
#DROME vs DROPS
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/DROME_DROPS_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/ortholog_DROME_DROPS_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/final_pair_DROME_DROPS_domain_loss')

#DROME vs DROYA
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/DROME_DROYA_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/ortholog_DROME_DROYA_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/final_pair_DROME_DROYA_domain_loss')

#DROME vs DROAN
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/DROME_DROAN_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/ortholog_DROME_DROAN_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/final_pair_DROME_DROAN_domain_loss')

#DROME vs DROMO
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/DROME_DROMO_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/ortholog_DROME_DROMO_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/final_pair_DROME_DROMO_domain_loss')

#DROME vs DROSI
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/DROME_DROSI_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/ortholog_DROME_DROSI_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/final_pair_DROME_DROSI_domain_loss')

#DROME vs DROVI
extract_ortholog_domain_loss_domainDiff(domain_loss_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/DROME_DROVI_domain_loss', 
                                        domain_modification_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/ortholog_DROME_DROVI_domain_modifications', 
                                        output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/final_pair_DROME_DROVI_domain_loss')
