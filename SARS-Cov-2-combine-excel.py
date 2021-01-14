#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 22:20:33 2020

@author: mikkel
"""

# Example of script combining all of the excel files outputted from NetMHCpan-4.1.

import os
import pandas as pd
cwd = os.path.abspath('') 
files = os.listdir(cwd)

# Collect data from all of the excel files and store in a dataframe
df = pd.DataFrame()
for file in files:
    if file.endswith('.xls'):
        current_file = pd.read_excel(file, header = None)
        current_hla = current_file.iat[0,3]
        current_file = current_file.drop(current_file.index[0:2])
        current_file['HLA'] = current_hla
        df = df.append(current_file, ignore_index=True)
        
# Give column names and export to csv file
df.columns = ['Pos', 'Peptide', 'ID', 'Core', 'Icore', 'EL-score', 'EL_rank', 'Ave', 'NB', 'HLA']
df.to_csv('SARS-Cov-2_combined.csv')
