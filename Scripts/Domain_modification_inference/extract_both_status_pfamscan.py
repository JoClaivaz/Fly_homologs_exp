# -*- coding: utf-8 -*-
"""
Created on Tue May 23 20:12:16 2017

@author: Claivaz
Extract domain pfamscan of a pair which is present in DROME_DROYA_status as both (meaning both are present in pfamscan at least for one domain)
"""
def extract_both_status_pfamscan(status_file, pfamscan_parsed_in, output_file):  
    
    status_in = open(status_file, 'r')
    pair_number = 0
    pair_both_present = open(output_file, 'w')

    for status_line in status_in:
        if 'both' in status_line:
            status_line_tmp = status_line.replace('"','')
            pair_number += 1
            pfamscan_in = open(pfamscan_parsed_in, 'r')
            for pfamscan_line in pfamscan_in:
                if status_line_tmp.split(',')[1] in pfamscan_line:
                    pair_both_present.write(str(pair_number) + '\t' + pfamscan_line.split('\t')[1] + 
                                            '\t' + pfamscan_line.split('\t')[2] + 
                                            '\t' + pfamscan_line.split('\t')[3] + 
                                            '\t' + pfamscan_line.split('\t')[4] + 
                                            '\t' + pfamscan_line.split('\t')[5])
                elif status_line_tmp.split(',')[3] in pfamscan_line:
                    pair_both_present.write(str(pair_number) + '\t' + pfamscan_line.split('\t')[1] + 
                                            '\t' + pfamscan_line.split('\t')[2] + 
                                            '\t' + pfamscan_line.split('\t')[3] + 
                                            '\t' + pfamscan_line.split('\t')[4] + 
                                            '\t' + pfamscan_line.split('\t')[5])
            pfamscan_in.close()
    pair_both_present.close()
    status_in.close()

#Run FUN in Windows
#DROMEvsDROYA
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROYA_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROYA_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROYA')

#DROMEvsDROPS
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROPS_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROPS_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROPS')

#DROMEvsDROAN
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROAN_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROAN_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROAN')

#DROMEvsDROMO
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROMO_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROMO_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROMO')

#DROMEvsDROSI
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROSI_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROSI_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROSI')
#DROMEvsDROVI
extract_both_status_pfamscan(status_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/DROME_DROVI_status', 
                             pfamscan_parsed_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_ortholog_DROME_DROVI_domain_parsed', 
                             output_file = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_both_orthologPair_DROME_DROVI')

