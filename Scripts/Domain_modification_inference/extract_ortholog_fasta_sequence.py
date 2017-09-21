# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170512
Extraction of fasta sequence amongst orthologs pairs between specie 1 and specie 2 
(obtained from OMA), in function of the pairwise_ortholog file, only considering 
subset of one parameter (e.g. '1:1', unique ortholog, ':' all pair).
"""

###Function
def extract_ortholog_fasta_sequence (Specie1, Specie2, PairOrtho, output, parameter):

    pair_ortholog = open(PairOrtho, 'r')
    unique_ortholog = open(output, 'w')
    write = False
    spec1_tmp = []
    spec2_tmp = []
    
    for pair in pair_ortholog:
        if parameter in pair.split('\t')[2] :
            SPEC1 = open(Specie1, 'r')
            SPEC2 = open(Specie2, 'r')
            
            if pair.split('\t')[0] not in spec1_tmp:
                for seq_spec1 in SPEC1:
                    if write == True and '>' in seq_spec1:
                        write = False
                        SPEC1.close()
                        break
                    elif write == True:
                        unique_ortholog.write(seq_spec1)
                    elif pair.split('\t')[0] in seq_spec1:
                        write = True
                        unique_ortholog.write(seq_spec1)
                        spec1_tmp.append(pair.split('\t')[0])
                
            if pair.split('\t')[1] not in spec2_tmp:        
                for seq_spec2 in SPEC2:
                    if write == True and '>' in seq_spec2:
                        write = False                
                        SPEC2.close()  
                        break
                    elif write == True:
                        unique_ortholog.write(seq_spec2)
                    elif pair.split('\t')[1] in seq_spec2:
                        write = True
                        unique_ortholog.write(seq_spec2)
                        spec2_tmp.append(pair.split('\t')[1])
                                
    
       
    pair_ortholog.close()
    unique_ortholog.close()
    

#Run the function / WINDOWS
#DROME vs DROYA, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROYA/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROYA/OMA.2.1.1/DB/DROYA.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROYA/pairwise_ortholog_DROME_DROYA', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROYA', 
                                parameter = ':')

#DROME vs DROPS, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROPS/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROPS/OMA.2.1.1/DB/DROPS.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROPS/pairwise_ortholog_DROME_DROPS', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROPS', 
                                parameter = ':')     

#DROME vs DROAN, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROAN.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROAN', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROAN', 
                                parameter = ':')     

#DROME vs DROMO, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROMO.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROMO', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROMO', 
                                parameter = ':')    

#DROME vs DROSI, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROSI.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROSI', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROSI', 
                                parameter = ':')     
 
#DROME vs DROVI, all orthologs
extract_ortholog_fasta_sequence(Specie1 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/OMA.2.1.1/DB/DROVI.fa', 
                                PairOrtho = 'D:/UNIL/Master/Master_Project/OMA/DROME_DROAN_DROMO_DROSI_DROVI/pairwise_ortholog_DROME_DROVI', 
                                output = 'D:/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROVI', 
                                parameter = ':')     




#Run the function / UBUNTU
#DROME vs DROYA, all orthologs
extract_ortholog_fasta_sequence(Specie1 = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROYA/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROYA/OMA.2.1.1/DB/DROYA.fa', 
                                PairOrtho = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROYA/pairwise_ortholog_DROME_DROYA', 
                                output = '/media/jclaivaz/Data/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROYA', 
                                parameter = ':')

#DROME vs DROPS, all orthologs
extract_ortholog_fasta_sequence(Specie1 = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROPS/OMA.2.1.1/DB/DROME.fa',
                                Specie2 = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROPS/OMA.2.1.1/DB/DROPS.fa', 
                                PairOrtho = '/media/jclaivaz/Data/UNIL/Master/Master_Project/OMA/DROME_DROPS/pairwise_ortholog_DROME_DROPS', 
                                output = '/media/jclaivaz/Data/UNIL/Master/Master_Project/Data/OMA/ortholog_DROME_DROPS', 
                                parameter = ':')     
