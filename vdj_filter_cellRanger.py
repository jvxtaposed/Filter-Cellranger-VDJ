# -*- coding: utf-8 -*-
"""
@author: Kai Fung (kaitious.fung@gmail.com)
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfile', required=True, type=str)
parser.add_argument('-o', '--outputfile', required=True, type=str)
args = parser.parse_args()

# Read cellranger output into pandas dataframe
cellranger_out = pd.read_csv(args.inputfile)

# Filter out nonproductive genes
df_productive = cellranger_out[cellranger_out['productive'].str.lower() == 'true']

# Initialize dictionary to store barcode info for each chain
chain_info = {'IGH': {}, 'IGL': {}, 'IGK': {}}

# Iterate over rows in the dataframe
for _, row in df_productive.iterrows():
    chain = row['chain']
    barcode = row['barcode']
    if barcode not in chain_info[chain]:
        chain_info[chain][barcode] = {}
        for col in ['contig_id', 'cdr3', 'v_gene', 'd_gene', 'j_gene', 'c_gene']:
            if row[col] != 'None':
                chain_info[chain][barcode][f'{chain}_{col}'] = row[col]
    else:
        prefix = f'{chain}_{col}'
        for col in ['contig_id', 'cdr3', 'v_gene', 'd_gene', 'j_gene', 'c_gene']:
            if row[col] != 'None':
                chain_info[chain][barcode][f'{prefix}2'] = row[col]

# Convert dictionary to dataframe
cellrange_reformated = pd.DataFrame.from_dict({(chain, barcode): values 
                                                for chain, barcodes in chain_info.items() 
                                                for barcode, values in barcodes.items()},
                                                orient='index').reset_index()

# Rename columns and write result to CSV file
cellrange_reformated.rename(columns={'level_0': 'chain', 'level_1': 'barcode'}, inplace=True)
cellrange_reformated.to_csv(args.outputfile, index=False)
