# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:32:20 2017

@author: Claivaz
This script formates the expression data of DROME from flybase, allowing to be read in R
"""

flybase_file = open('D:/UNIL/Master/Master_Project/Data/flybase/gene_rpkm_report_fb_2017_04.tsv', 'r')
flybase_file_formate = open('D:/UNIL/Master/Master_Project/Data/flybase/gene_rpkm_report_fb_2017_04_formated', 'w')

flybase_file_formate.write('##Release_ID\tFBgn#\tGeneSymbol\tParent_library_FBlc#\tParent_library_name\tRNASource_FBlc#\tRNASource_name\tRPKM_value\tBin_value\tUnique_exon_base_count\tTotal_exon_base_count\tCount_used\n'.replace('#', ''))

for fb_line in flybase_file:
    if '#' not in fb_line and 'WARNING' not in fb_line:
        flybase_file_formate.write(fb_line)
        
flybase_file.close()
flybase_file_formate.close()
    




