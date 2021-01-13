#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 23:43:58 2020

@author: mikkel
"""

import os
import pandas as pd

cwd = os.path.abspath('') 
files = os.listdir(cwd)

df = pd.DataFrame()

for file in files:
    if file.endswith('.txt'):
        
        with open(file, "r") as f:
            lines = f.readlines()
        
         
        with open(file, "w") as f:
            for line in lines:
                if not line.startswith('#'):
                    f.write(line)
            
            current_df = pd.read_csv(file, sep = " ", header=None)
            current_hla = file[27:37]
            current_df['HLA'] = current_hla
            df = df.append(current_df, ignore_index=True)


df.columns = ["Rank", "SB in SARS-Cov-2", "SB in HCov-HKU1", "BLOSUM score", "HLA"]
df.to_csv('blosum_SARS-Cov-2_HCov-HKU1.csv')
    
    