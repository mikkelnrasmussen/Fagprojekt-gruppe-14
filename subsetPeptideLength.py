#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 19:02:34 2020

@author: mikkel
"""

import pandas as pd
import numpy as np

# Function dividing the HLA-ligands into lengths of 8-11 amino acids and 
# outputting them in csv files

def subsetPeptideLength(virus, virus_name):
    
     # List with all the 38 HLA alleles
     hla_list = ['HLA-A02:01', 'HLA-A01:01', 'HLA-A03:01', 'HLA-A24:02', 'HLA-A11:01', 'HLA-A26:01', 'HLA-A32:01', 
                 'HLA-A68:01', 'HLA-A25:01', 'HLA-A31:01', 'HLA-A29:02', 'HLA-A23:01', 'HLA-B07:02', 'HLA-B08:01', 
                 'HLA-B15:01', 'HLA-B51:01', 'HLA-B44:02', 'HLA-B18:01', 'HLA-B35:01', 'HLA-B44:03', 'HLA-B40:01', 
                 'HLA-B13:02', 'HLA-B27:05', 'HLA-B57:01', 'HLA-B35:03', 'HLA-B38:01', 'HLA-B58:01', 'HLA-C07:01', 
                 'HLA-C04:01', 'HLA-C07:02', 'HLA-C06:02', 'HLA-C12:03', 'HLA-C05:01', 'HLA-C02:02', 'HLA-C03:04', 
                 'HLA-C03:03', 'HLA-C01:02', 'HLA-C15:02']
    
     for i in hla_list:
        # Selecting the unique HLA-ligands for each HLA allele
        temp_virus_hla = virus.loc[virus['HLA'] == i]['Icore'].unique()
        
        # Vectorization of the "len" function
        peptide_len = np.vectorize(len)
        
        # Subdivide the HLA-ligands by length (8-11)
        virus_len_8 = temp_virus_hla[peptide_len(temp_virus_hla) == 8]
        virus_len_9 = temp_virus_hla[peptide_len(temp_virus_hla) == 9]
        virus_len_10 = temp_virus_hla[peptide_len(temp_virus_hla) == 10]
        virus_len_11 = temp_virus_hla[peptide_len(temp_virus_hla) == 11]
        
        # Store divided HLA-ligands in dataframes
        df_len_8 = pd.DataFrame({'#' + virus_name + ' Icore (8-mer)': virus_len_8})
        df_len_9 = pd.DataFrame({'#' + virus_name + ' Icore (9-mer)': virus_len_9})
        df_len_10 = pd.DataFrame({'#' + virus_name + ' Icore (10-mer)': virus_len_10})
        df_len_11 = pd.DataFrame({'#' + virus_name + ' Icore (11-mer)': virus_len_11})
        
        # Changing the ":" to a "-", e.g. HLA-A02:01 to HLA-A02-01
        i = i[0:7] + '-' + i[8:]
        
        # Creating filenames
        filename_8mer = virus_name + '_' + i + '_8-mer'
        filename_9mer = virus_name + '_' + i + '_9-mer'
        filename_10mer = virus_name + '_' + i + '_10-mer'
        filename_11mer = virus_name + '_' + i + '_11-mer'
        
        # Selecting the correct directory
        if (virus_name == 'SARS-Cov-2'):
            folder_name = 'SARS-Cov-2-length-separated/'
        elif (virus_name == 'HCov-229E'):
            folder_name = 'HCov-229E-length-separated/'  
        elif (virus_name == 'HCov-HKU1'):
            folder_name = 'HCov-HKU1-length-separated/'     
        elif (virus_name == 'HCov-NL63'):
            folder_name = 'HCov-NL63-length-separated/'
        elif (virus_name == 'HCov-OC43'):
            folder_name = 'HCov-OC43-length-separated/'
        elif (virus_name == 'Zaire-ebolavirus'):
            folder_name = 'Zaire-ebolavirus-length-separated/'
        elif (virus_name == 'Influenza-virus-A-H3N2'):
            folder_name = 'Influenza-virus-A-H3N2-length-separated/'
        elif (virus_name == 'SARS-Cov-2_spike'):
            folder_name = 'pep2score/SARS-Cov-2_spike/'
        elif (virus_name == 'SARS-Cov-2_nucleocapsid'):
            folder_name = 'pep2score/SARS-Cov-2_nucleocapsid/'
        
        # Export divided HLA-ligands in csv files
        df_len_8.to_csv('/Users/mikkel/Desktop/NetMHCpan-4.1/' + folder_name + filename_8mer + '.csv', index = False)
        df_len_9.to_csv('/Users/mikkel/Desktop/NetMHCpan-4.1/' + folder_name + filename_9mer + '.csv', index = False)
        df_len_10.to_csv('/Users/mikkel/Desktop/NetMHCpan-4.1/' + folder_name + filename_10mer + '.csv', index = False)
        df_len_11.to_csv('/Users/mikkel/Desktop/NetMHCpan-4.1/' + folder_name + filename_11mer + '.csv', index = False)
        
    
# Import of csv files containing predicted HLA-ligand from NetMHCpan for the five viruses + the two controls; 
# SARS-Cov-2, HCov-229E, HCov-HKU1, HCov-NL63 and HCov-OC43, Zaire EBOV and Influenza-A H3N2.  
df_SARS = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')
df_EBOV = pd.read_csv('Zaire-ebolavirus/Zaire-ebolavirus_combined.csv')
df_H3N2 = pd.read_csv('Influenza-virus-A-H3N2/Influenza-virus-A-H3N2_combined.csv')
df_SARS_spike = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_spike.csv')
df_SARS_nucleocapsid = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_nucleocapsid.csv')

# Restrict data to only HLA-ligands (strong binders (SB)) EL_rank < 0.5
df_SARS_SB = df_SARS.loc[df_SARS['EL_rank'] < 0.5]
df_229E_SB = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_HKU1_SB = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_NL63_SB = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_SB = df_OC43.loc[df_OC43['EL_rank'] < 0.5]
df_EBOV_SB = df_EBOV.loc[df_EBOV['EL_rank'] < 0.5]
df_H3N2_SB = df_H3N2.loc[df_H3N2['EL_rank'] < 0.5]
df_SARS_spike_SB = df_SARS_spike.loc[df_SARS_spike['EL_rank'] < 0.5]
df_SARS_nucleocapsid_SB = df_SARS_nucleocapsid.loc[df_SARS_nucleocapsid['EL_rank'] < 0.5]

# Call subsetPeptideLength in order to divide HLA-ligands based on length
subsetPeptideLength(df_SARS_SB, 'SARS-Cov-2')
subsetPeptideLength(df_229E_SB, 'HCov-229E')
subsetPeptideLength(df_HKU1_SB, 'HCov-HKU1')
subsetPeptideLength(df_NL63_SB, 'HCov-NL63')
subsetPeptideLength(df_OC43_SB, 'HCov-OC43')
subsetPeptideLength(df_EBOV_SB, 'Zaire-ebolavirus')
subsetPeptideLength(df_H3N2_SB, 'Influenza-virus-A-H3N2')
subsetPeptideLength(df_SARS_spike_SB, 'SARS-Cov-2_spike')
subsetPeptideLength(df_SARS_nucleocapsid_SB, 'SARS-Cov-2_nucleocapsid')



