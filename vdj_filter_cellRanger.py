# -*- coding: utf-8 -*-
"""
@author: Kai Fung (kaitious.fung@gmail.com)
"""

import pandas as pd
import numpy as np
import scipy.io as sp
import csv
import os
import argparse


#python3 vdj_filter_cellRanger.py -i test_CellRangerVDJ.csv -o TEST_OUT.csv
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfile',required=True,type=str)
parser.add_argument('-o', '--outputfile',required=True,type=str)
args = parser.parse_args()

#'''reads in cellranger output into pandas dataframe'''
cellranger_out = pd.read_csv(args.inputfile)

#print(cellranger_out.loc[:, ['barcode', 'contig_id','chain','v_gene','d_gene','j_gene','c_gene','productive']])

#'''creates a boolean array that is then used to filter out the nonproductive genes'''
is_productive = cellranger_out['productive']==("True" or "TRUE")
#print(cellranger_out['productive'])
df_productive = cellranger_out[is_productive]

#'''prints a truncated matrix with the relevant information'''
#print(df_productive.loc[:, ['barcode', 'contig_id','chain','v_gene','d_gene','j_gene','c_gene','productive']])


inserted_barcodes_IGH = set()
inserted_barcodes_IGL = set()
inserted_barcodes_IGK = set()
temp_input = []
for i in range(1,len(df_productive)):
    if df_productive.chain.iloc[i-1] == 'IGH':
        dict_IGH = {}
        dict_IGH['barcode'] = df_productive.barcode.iloc[i-1]
        if (df_productive.barcode.iloc[i-1] not in inserted_barcodes_IGH):
            inserted_barcodes_IGH.add(df_productive.barcode.iloc[i-1])
            if df_productive.contig_id.iloc[i-1] != 'None':
                dict_IGH['contig_id_IGH'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
                dict_IGH['CDR3_IGH'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_v'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_d'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_j'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_c'] = df_productive.c_gene.iloc[i-1]
        else:
            if df_productive.contig_id.iloc[i-1] != 'None':
                dict_IGH['contig_id_IGH2'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
                dict_IGH['CDR3_IGH2'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_v2'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_d2'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_j2'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGH['IGH_c2'] = df_productive.c_gene.iloc[i-1]
        temp_input.append(dict_IGH)

    elif df_productive.chain.iloc[i-1] == 'IGK':
        dict_IGK = {}
        dict_IGK['barcode'] = df_productive.barcode.iloc[i-1]
        if (df_productive.barcode.iloc[i-1] not in inserted_barcodes_IGK):
            inserted_barcodes_IGK.add(df_productive.barcode.iloc[i-1])
            if df_productive.contig_id.iloc[i-1] != 'None':
                dict_IGK['contig_id_IGK'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
                dict_IGK['CDR3_IGK'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_v'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_d'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_j'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_c'] = df_productive.c_gene.iloc[i-1]
        else:
            if df_productive.contig_id.iloc[i-1] != 'None':
                dict_IGK['contig_id_IGK2'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
                dict_IGK['CDR3_IGK2'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_v2'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_d2'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_j2'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGK['IGK_c2'] = df_productive.c_gene.iloc[i-1]
        temp_input.append(dict_IGK)
    elif df_productive.chain.iloc[i-1] == 'IGL':
        dict_IGL = {}
        dict_IGL['barcode'] = df_productive.barcode.iloc[i-1]
        if (df_productive.barcode.iloc[i-1] not in inserted_barcodes_IGL):
            inserted_barcodes_IGL.add(df_productive.barcode.iloc[i-1])
            if df_productive.contig_id.iloc[i-1] != 'None':
             dict_IGL['contig_id_IGL'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
             dict_IGL['CDR3_IGL'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_v'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_d'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_j'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_c'] = df_productive.c_gene.iloc[i-1]
        else:
            if df_productive.contig_id.iloc[i-1] != 'None':
             dict_IGL['contig_id_IGL2'] = df_productive.contig_id.iloc[i-1]
            if df_productive.cdr3.iloc[i-1] != 'None':
             dict_IGL['CDR3_IGL2'] = df_productive.cdr3.iloc[i-1]
            if df_productive.v_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_v2'] = df_productive.v_gene.iloc[i-1]
            if df_productive.d_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_d2'] = df_productive.d_gene.iloc[i-1]
            if df_productive.j_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_j2'] = df_productive.j_gene.iloc[i-1]
            if df_productive.c_gene.iloc[i-1] != 'None':
                dict_IGL['IGL_c2'] = df_productive.c_gene.iloc[i-1]
        temp_input.append(dict_IGL)
#print(temp_input)

cellrange = pd.DataFrame(temp_input)
cellrange = cellrange.replace(np.nan, '', regex=True)
#print(cellrange)

cellrange_reformated = cellrange.groupby('barcode').agg(lambda x: ' '.join(x.unique())).reset_index()

print(cellrange_reformated)
cellrange_reformated.to_csv('cellranger_reformated.csv')


#----------------------------
''' CITE-SEQ COUNT
with open('features.tsv', 'r') as f:
    reader = csv.reader(f)
    ugh = list(reader)

features = [item for sublist in ugh for item in sublist]

#print(features)

with open('barcodes.tsv', 'r') as b:
    read = csv.reader(b)
    ugh1 = list(read)

barcodes = [item for sublist in ugh1 for item in sublist]
#print(barcodes)



matrix = sp.mmread('matrix.mtx')
B = matrix.todense()
count_out = pd.DataFrame(B)
count_out.columns = barcodes
count_out.index = features
count_out_transposed = count_out.T
count_out_transposed = count_out_transposed.rename_axis('barcode').reset_index()


result = pd.merge(count_out_transposed, cellrange_reformated, on='barcode')
print(result)

result.to_csv('CiteSeqCountTestResult.csv')
'''
#----------------------------

with open('features.tsv', 'r') as f:
    reader = csv.reader(f)
    ugh = list(reader)
features = [item for sublist in ugh for item in sublist]


with open('barcodes.tsv', 'r') as b:
    read = csv.reader(b)
    ugh1 = list(read)
barcodes = [item for sublist in ugh1 for item in sublist]
#print(barcodes)


matrix = sp.mmread('matrix.mtx')
B = matrix.todense()
count_out = pd.DataFrame(B)
count_out.columns = barcodes
count_out.index = features
count_out_transposed = count_out.T
count_out_transposed = count_out_transposed.rename_axis('barcode').reset_index()

features_NoGeneExp = [c for c in count_out_transposed.columns if c[-15:] != 'Gene Expression']
#print(features_NoGeneExp)

count_out_transposed = count_out_transposed[features_NoGeneExp]
result = pd.merge(count_out_transposed, cellrange_reformated, on='barcode')
print(result)


#result.to_csv('count_CELLRANGER_Result_CDR3_CONTIG_ID_v2.csv')
#result.to_csv('TEST_OUT.csv')

#result.to_csv(args.outputfile)
