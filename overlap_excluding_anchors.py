#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:17:38 2020

@author: mikkel
"""
import pandas as pd
import numpy as np

# Function calculating the number an proportional overlap of identical HLA-ligands, where anchor points has been excluded,
# for two viruses being compared.

def overlap_excluding(virus1, virus2, virus1_name, virus2_name):
    
    
    hla_list = ['HLA-A02:01', 'HLA-A01:01', 'HLA-A03:01', 'HLA-A24:02', 'HLA-A11:01', 'HLA-A26:01', 
                'HLA-A32:01', 'HLA-A68:01', 'HLA-A25:01', 'HLA-A31:01', 'HLA-A29:02', 'HLA-A23:01', 
                'HLA-B07:02', 'HLA-B08:01', 'HLA-B15:01', 'HLA-B51:01', 'HLA-B44:02', 'HLA-B18:01', 
                'HLA-B35:01', 'HLA-B44:03', 'HLA-B40:01', 'HLA-B13:02', 'HLA-B27:05', 'HLA-B57:01', 
                'HLA-B35:03', 'HLA-B38:01', 'HLA-B58:01', 'HLA-C07:01', 'HLA-C04:01', 'HLA-C07:02', 
                'HLA-C06:02', 'HLA-C12:03', 'HLA-C05:01', 'HLA-C02:02', 'HLA-C03:04', 'HLA-C03:03', 
                'HLA-C01:02', 'HLA-C15:02']
    
    # Creating dataframe to store calculated values for every HLA allele
    df_virus_compare = pd.DataFrame(columns = [virus1_name + ' SB', virus2_name + ' SB', 'Identical Icores', 'Proportional overlap'], index = hla_list)
    

    for i in hla_list:
        
        # Selecting the rows containing the unique I-cores (HLA-ligands) for each HLA allele
        temp_virus1_hla = virus1.loc[virus1['HLA'] == i]['Icore'].unique()
        temp_virus2_hla = virus2.loc[virus2['HLA'] == i]['Icore'].unique()
        
        # Remove anchor points
        # Remove the second amino acid in the peptide sequences
        temp_list_virus1 = [i[:1] + i[2:] for i in temp_virus1_hla]
        temp_list_virus2 = [i[:1] + i[2:] for i in temp_virus2_hla]
        
        # Remove the last amino acid in the peptide sequences
        temp_list_virus1 = [i[:-1] for i in temp_list_virus1]
        temp_list_virus2 = [i[:-1] for i in temp_list_virus2]

        # Converting to numpy array and calculating the number og HLA-ligands (SB) 
        # in the two viruses
        temp_virus1_hla = np.array(temp_list_virus1)
        temp_virus2_hla = np.array(temp_list_virus2)
        size_of_SB_virus1 = np.size(temp_virus1_hla)
        size_of_SB_virus2 = np.size(temp_virus2_hla)
        
        # Calculating the number of identical overlap, as testing if each sequences in virus 1 
        # are identical with a seqeunce in virus 2
        identical_icores = np.in1d(temp_virus1_hla, temp_virus2_hla)
        size_of_identical_icores = np.sum(identical_icores)
        
        # Selecting the peptide space as the minimun of HLA-ligands found in either 
        # virus 1 or virus 2
        peptide_space = min(size_of_SB_virus1, size_of_SB_virus2)
        
        # Calculating the proportional overlap
        prop_overlap =  size_of_identical_icores / peptide_space * 100
        
        # Store values as a row in the dataframe
        df_virus_compare.loc[i] = [size_of_SB_virus1 , size_of_SB_virus2, size_of_identical_icores, prop_overlap]
 
    return df_virus_compare
    
    
