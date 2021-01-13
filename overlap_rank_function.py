#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:18:25 2020

@author: mikkel
"""

import pandas as pd
import numpy as np


def overlap_rank(virus1, virus2, virus1_name, virus2_name):
    
    # List with all the 38 HLA alleles
    hla_list = ['HLA-A02:01', 'HLA-A01:01', 'HLA-A03:01', 'HLA-A24:02', 'HLA-A11:01', 
                'HLA-A26:01', 'HLA-A32:01', 'HLA-A68:01', 'HLA-A25:01', 'HLA-A31:01', 
                'HLA-A29:02', 'HLA-A23:01', 'HLA-B07:02', 'HLA-B08:01', 'HLA-B15:01', 
                'HLA-B51:01', 'HLA-B44:02', 'HLA-B18:01', 'HLA-B35:01', 'HLA-B44:03', 
                'HLA-B40:01', 'HLA-B13:02', 'HLA-B27:05', 'HLA-B57:01', 'HLA-B35:03', 
                'HLA-B38:01', 'HLA-B58:01', 'HLA-C07:01', 'HLA-C04:01', 'HLA-C07:02', 
                'HLA-C06:02', 'HLA-C12:03', 'HLA-C05:01', 'HLA-C02:02', 'HLA-C03:04', 
                'HLA-C03:03', 'HLA-C01:02', 'HLA-C15:02']
    
    # Create a dataframe with the needed columns 
    # (strong binders (SBs) in the two viruses, actual overlap and proportional overlap)
    df_virus_compare = pd.DataFrame(columns = [virus1_name + ' SB', virus2_name + 
                                               ' SB', 'Identical Icores', 
                                               'Proportional overlap'], 
                                    index = hla_list)
    

    for i in hla_list:
        
        # Restrict to only unique icores (peptide sequences) for a given HLA
        temp_virus1_hla = virus1.loc[virus1['HLA'] == i]['Icore'].unique()
        temp_virus2_hla = virus2.loc[virus2['HLA'] == i]['Icore'].unique()
        
        # Calculate the number of SBs by taking the length of two list
        size_of_SB_virus1 = np.size(temp_virus1_hla)
        size_of_SB_virus2 = np.size(temp_virus2_hla)
        
        # Calculate the actual overlap (identical icores)
        identical_icores = np.in1d(temp_virus1_hla, temp_virus2_hla)
        size_of_identical_icores = np.sum(identical_icores)
        
        # Find the peptide space by finding the smallest value of the two SBs 
        peptide_space = min(size_of_SB_virus1, size_of_SB_virus2)
        
        # Calculate the proportional overlap
        prop_overlap =  size_of_identical_icores / peptide_space * 100
        
        # Store in dataframe
        df_virus_compare.loc[i] = [size_of_SB_virus1 , size_of_SB_virus2, 
                                   size_of_identical_icores, prop_overlap]
 
    return df_virus_compare
    
    