# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:54:09 2017

@author: Claivaz
this script allows to determine the different available states for each pair species 
"""
#DROMEvsDROAN
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROAN', 'r')

DROME_list = []
DROAN_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROAN' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROAN_list:
            DROAN_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROAN_list

#DROMEvsDROMO
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROMO', 'r')

DROME_list = []
DROMO_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROMO' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROMO_list:
            DROMO_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROMO_list

#DROMEvsDROPS
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROPS', 'r')

DROME_list = []
DROPS_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROPS' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROPS_list:
            DROPS_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROPS_list

#DROMEvsDROSI
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROSI', 'r')

DROME_list = []
DROSI_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROSI' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROSI_list:
            DROSI_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROSI_list

#DROMEvsDROVI
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROVI', 'r')

DROME_list = []
DROVI_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROVI' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROVI_list:
            DROVI_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROVI_list

#DROMEvsDROYA
dataset_exp = open('D:/UNIL/Master/Master_Project/Data/expression_analysis/data_expression_DROME_DROYA', 'r')

DROME_list = []
DROYA_list = []

for dataset_line in dataset_exp:
    if 'DROME' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROME_list:
            DROME_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])
    elif 'DROYA' in dataset_line:
        if [dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]] not in DROYA_list:
            DROYA_list.append([dataset_line.split('\t')[3], dataset_line.split('\t')[4], dataset_line.split('\t')[5]])

dataset_exp.close()

DROME_list
DROYA_list