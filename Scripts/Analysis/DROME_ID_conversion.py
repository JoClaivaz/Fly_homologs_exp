# -*- coding: utf-8 -*-
"""
Joaquim Claivaz
170724
Create file containing the different identifier for DROME genes: DROMEXXX FBgnXXX CGXXX
"""

path_file = 'D:/UNIL/Master/Master_Project/Data/expression_analysis/all_droso/'
FlyBase_IDs = open(path_file + 'FlyBase_IDs.txt', 'r')
DROME_FBgn_names = open(path_file + 'DROME_FBgn_names', 'r')
dict_flybase = {}

#skip the headers
FlyBase_IDs.readline()

for line_flybase in FlyBase_IDs:
    line_tmp = line_flybase.replace('\n', '').replace("'", '').split('\t')
    dict_flybase[line_tmp[0]] = line_tmp[3]
FlyBase_IDs.close()
    
DROME_id_converter = open(path_file + 'DROME_id_converter', 'w')

for DROME_line in DROME_FBgn_names:
    line_split = DROME_line.replace('\n', '').split('\t')
    DROME_id_converter.write('%s\t%s\t%s\n' % (line_split[0], line_split[1], dict_flybase[line_split[1]]))
    
DROME_FBgn_names.close()
DROME_id_converter.close()