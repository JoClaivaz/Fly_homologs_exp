# -*- coding: utf-8 -*-
'''
Joaquim Claivaz
170512
Extraction of orthologs pairs between specie 1 and specie 2
without domain modification (obtained from output of pfamscan), with 1 domain modification 
and with complex modification (more than 1 domain differences or same length but different composition).
'''

def extract_ortholog_modification(PairOrtho_in, Domain_file_in, Domain_file_out):
    
    PairOrtho = open(PairOrtho_in, 'r')
    Domain_file_nomodif = open(Domain_file_out + '_domain_nomodif', 'w')
    Domain_file_modif_complex = open(Domain_file_out + '_domain_complex_modif', 'w')
    Domain_file_modif_1 = open(Domain_file_out + '_1_domain_modif', 'w')
    
    for pair_ortho in PairOrtho:
        spec1_domain = []
        spec2_domain = []
        Domain_file = open(Domain_file_in, 'r')
        for ortho_domain_line in Domain_file:
            if pair_ortho.split('\t')[0] in ortho_domain_line:
                ortho_domain_line = ortho_domain_line.replace(' ', '\t')
                while '\t\t' in ortho_domain_line:
                    ortho_domain_line = ortho_domain_line.replace('\t\t', '\t')
                spec1_domain.append(ortho_domain_line.split('\t')[6])
            elif pair_ortho.split('\t')[1] in ortho_domain_line:
                ortho_domain_line = ortho_domain_line.replace(' ', '\t')
                while '\t\t' in ortho_domain_line:
                    ortho_domain_line = ortho_domain_line.replace('\t\t', '\t')
                spec2_domain.append(ortho_domain_line.split('\t')[6])
        if spec1_domain == spec2_domain and len(spec1_domain) > 0:
            Domain_file_nomodif.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec1_domain)) +
                                      '\t' + pair_ortho.split('\t')[1] + '\tnomodif\n')
        elif len(spec1_domain) == len(spec2_domain) and len(spec1_domain) > 0:
            Domain_file_modif_complex.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec1_domain)) +
                                                '\t' + pair_ortho.split('\t')[1] + '\tcomplex_modif\n')
        elif len(spec1_domain) > 0 and len(spec2_domain) > 0:
            if len(spec1_domain)+1 == len(spec2_domain) or len(spec2_domain)+1 == len(spec1_domain):
                Domain_file_modif_1.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec1_domain)) +
                                              '\t' + pair_ortho.split('\t')[1] + '\tmodif\n')
            else:
                Domain_file_modif_complex.write(pair_ortho.split('\t')[0] + '\t' + str(len(spec1_domain)) +
                                                '\t' + pair_ortho.split('\t')[1] + '\tcomplex_modif\n')
        
                
    
    PairOrtho.close()
    Domain_file_nomodif.close()
    Domain_file_modif_complex.close()
    Domain_file_modif_1.close()

#Run the function / WINDOWS
#DROME vs DROPS
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROPS',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/ortholog_DROME_DROPS_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROPS/ortholog_DROME_DROPS')

#DROME vs DROYA
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROYA',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/ortholog_DROME_DROYA_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROYA/ortholog_DROME_DROYA')

#DROME vs DROAN
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROAN',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/ortholog_DROME_DROAN_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROAN/ortholog_DROME_DROAN')

#DROME vs DROMO
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROMO',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/ortholog_DROME_DROMO_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROMO/ortholog_DROME_DROMO')

#DROME vs DROSI
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROSI',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/ortholog_DROME_DROSI_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROSI/ortholog_DROME_DROSI')

#DROME vs DROVI
extract_ortholog_modification(PairOrtho_in = 'D:/UNIL/Master/Master_Project/Data/OMA/pairwise_ortholog_DROME_DROVI',
                              Domain_file_in = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/ortholog_DROME_DROVI_domain', 
                              Domain_file_out = 'D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/all_DROSO/DROME_DROVI/ortholog_DROME_DROVI')