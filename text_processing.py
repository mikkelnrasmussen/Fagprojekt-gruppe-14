#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 23:43:58 2020

@author: mikkel
"""

import os
import pandas as pd

# Example of processing the txt files outputted from pep2score (BLOSUM analysis)

cwd = os.path.abspath('') 
files = os.listdir(cwd)

# Creating dataframe for content of each txt file to be stored in
df = pd.DataFrame()

for file in files:
    # Select the .txt files
    if file.endswith('.txt'):
        
        # Reading all lines in the file
        with open(file, "r") as f:
            lines = f.readlines()
        
        # For every line in the file only write the lines that does not start with #
        # Removed text in the start of every txt file
        with open(file, "w") as f:
            for line in lines:
                if not line.startswith('#'):
                    f.write(line)
            
            # Store information in dataframe and collect HLA allele info from title of the file
            current_df = pd.read_csv(file, sep = " ", header=None)
            current_hla = file[27:37]
            current_df['HLA'] = current_hla
            df = df.append(current_df, ignore_index=True)


# Name columns and export to a csv file
df.columns = ["Rank", "SB in SARS-Cov-2", "SB in HCov-HKU1", "BLOSUM score", "HLA"]
df.to_csv('blosum_SARS-Cov-2_HCov-HKU1.csv')
    
    
