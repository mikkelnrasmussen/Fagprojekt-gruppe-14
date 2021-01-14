#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:04:23 2020

@author: mikkel
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from venn import venn

# Import csv files containing predicted HLA-ligands from NetMHCpan for SARS-CoV-2 and the four common HCoVs
df_SARS = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')

# Restrict data to only HLA-ligands (strong binders (SB)) with EL_rank < 0.5
df_SARS_EL = df_SARS.loc[df_SARS['EL_rank'] < 0.5]
df_229E_EL = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_HKU1_EL = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_NL63_EL = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_EL = df_OC43.loc[df_OC43['EL_rank'] < 0.5]


hla_list = ['HLA-A02:01', 'HLA-A01:01', 'HLA-A03:01', 'HLA-A24:02', 
            'HLA-A11:01', 'HLA-A26:01', 'HLA-A32:01', 'HLA-A68:01', 
            'HLA-A25:01', 'HLA-A31:01', 'HLA-A29:02', 'HLA-A23:01', 
            'HLA-B07:02', 'HLA-B08:01', 'HLA-B15:01', 'HLA-B51:01', 
            'HLA-B44:02', 'HLA-B18:01', 'HLA-B35:01', 'HLA-B44:03', 
            'HLA-B40:01', 'HLA-B13:02', 'HLA-B27:05', 'HLA-B57:01', 
            'HLA-B35:03', 'HLA-B38:01', 'HLA-B58:01', 'HLA-C07:01', 
            'HLA-C04:01', 'HLA-C07:02', 'HLA-C06:02', 'HLA-C12:03', 
            'HLA-C05:01', 'HLA-C02:02', 'HLA-C03:04', 'HLA-C03:03', 
            'HLA-C01:02', 'HLA-C15:02']


for i in hla_list:
    
    # Selecting only unique HLA ligands (SB) with EL_rank < 0.5 for each HLA allele
    temp_SARS_hla = df_SARS_EL.loc[df_SARS_EL['HLA'] == i]['Icore'].unique()
    temp_229E_hla = df_229E_EL.loc[df_229E_EL['HLA'] == i]['Icore'].unique()
    temp_HKU1_hla = df_HKU1_EL.loc[df_HKU1_EL['HLA'] == i]['Icore'].unique()
    temp_NL63_hla = df_NL63_EL.loc[df_NL63_EL['HLA'] == i]['Icore'].unique()
    temp_OC43_hla = df_OC43_EL.loc[df_OC43_EL['HLA'] == i]['Icore'].unique()
    
    # If anchors points need to be excluded the # can be removed from the next lines of code
    # Removed the second amino acid in the peptide sequences
    #temp_list_SARS = [i[:1] + i[2:] for i in temp_SARS_hla]
    #temp_list_229E = [i[:1] + i[2:] for i in temp_229E_hla]
    #temp_list_HKU1 = [i[:1] + i[2:] for i in temp_HKU1_hla]
    #temp_list_NL63 = [i[:1] + i[2:] for i in temp_NL63_hla]
    #temp_list_OC43 = [i[:1] + i[2:] for i in temp_OC43_hla]
    
    # Removed the last amino acid in the peptide sequences        
    #temp_list_SARS = [i[:-1] for i in temp_list_SARS]
    #temp_list_229E = [i[:-1] for i in temp_list_229E]
    #temp_list_HKU1 = [i[:-1] for i in temp_list_HKU1]
    #temp_list_NL63 = [i[:-1] for i in temp_list_NL63]
    #temp_list_OC43 = [i[:-1] for i in temp_list_OC43]

    # Construct the comparison dictionary with sets, to be used by the venn function to create venn diagrams 
    comparison_dict = {'SARS-Cov-2': set(temp_list_SARS), 
                       'HCoV-229E': set(temp_list_229E), 
                       'HCoV-HKU1': set(temp_list_HKU1), 
                       'HCoV-NL63': set(temp_list_NL63), 
                       'HCoV-OC43': set(temp_list_OC43)}
    
    # Create venn diagrams and plot them
    venn(comparison_dict)
    plt.title(f'Comparison of SARS-Cov-2 with common HCoVs - {i} \n Excluding anchor residues')
    plt.savefig(f'venn_diagrams_without_anchors/Comparison SARS-Cov-2 and common HCoV without anchors - {i}')
    

