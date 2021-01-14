#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 22:29:57 2020

@author: mikkel
"""

import pandas as pd
from blosumOverlap import blosumOverlap

# Import of csv files containing predicted HLA-ligands for the five viruses and two controls; SARS-Cov-2, HCov-229E, HCov-HKU1, HCov-NL63 and HCov-OC43
# Zaire-EBOV and Influenza-A H3N2.
df_SARS = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')
df_EBOV = pd.read_csv('Zaire-ebolavirus/Zaire-ebolavirus_combined.csv')
df_H3N2 = pd.read_csv('Influenza-virus-A-H3N2/Influenza-virus-A-H3N2_combined.csv')

# Restrict data to only HLA-ligands (strong binders (SB)) - (EL_rank < 0.5) 
df_SARS_SB = df_SARS.loc[df_SARS['EL_rank'] < 0.5]
df_HKU1_SB = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_229E_SB = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_NL63_SB = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_SB = df_OC43.loc[df_OC43['EL_rank'] < 0.5]
df_EBOV_SB = df_EBOV.loc[df_EBOV['EL_rank'] < 0.5]
df_H3N2_SB = df_H3N2.loc[df_H3N2['EL_rank'] < 0.5] 

# Import of csv files with BLOSUM-scores obtained from the 6 comparisons with SARS-CoV-2, where anchors points were excluded
df_HKU1_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-HCov-HKU1-ex-anchors/blosum_SARS-Cov-2_HCov-HKU1_combined_ex_anchors.csv')
df_229E_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-HCov-229E-ex-anchors/blosum_SARS-Cov-2_HCov-229E_combined_ex_anchors.csv')
df_NL63_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-HCov-NL63-ex-anchors/blosum_SARS-Cov-2_HCov-NL63_combined_ex_anchors.csv')
df_OC43_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-HCov-OC43-ex-anchors/blosum_SARS-Cov-2_HCov-OC43_combined_ex_anchors.csv')
df_EBOV_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-Zaire-ebolavirus-ex-anchors/blosum_SARS-Cov-2_Zaire-ebolavirus_combined_ex_anchors.csv')
df_H3N2_blosum = pd.read_csv('pep2score/blosum-SARS-Cov-2-Influenza-virus-A-H3N2-ex-anchors/blosum_SARS-Cov-2_Influenza-virus-A-H3N2_combined_ex_anchors.csv')

# Calculate the number and proportional overlap of similar HLA-ligands 
# above a certain threshold BLOSUM-score (here 0.8) between two viruses.
df_SARS_HKU1 = blosumOverlap(df_HKU1_blosum, df_SARS_SB, df_HKU1_SB, 'SARS-Cov-2', 'HCov-HKU1')
df_SARS_229E = blosumOverlap(df_229E_blosum, df_SARS_SB, df_229E_SB, 'SARS-Cov-2', 'HCov-229E')
df_SARS_NL63 = blosumOverlap(df_NL63_blosum, df_SARS_SB, df_NL63_SB, 'SARS-Cov-2', 'HCov-NL63')
df_SARS_OC43 = blosumOverlap(df_OC43_blosum, df_SARS_SB, df_OC43_SB, 'SARS-Cov-2', 'HCov-OC43')
df_SARS_EBOV = blosumOverlap(df_EBOV_blosum, df_SARS_SB, df_EBOV_SB, 'SARS-Cov-2', 'Zaire Ebolavirus')
df_SARS_H3N2 = blosumOverlap(df_H3N2_blosum, df_SARS_SB, df_H3N2_SB, 'SARS-Cov-2', 'Influenza virus A H3N2')

# Export as csv file
df_SARS_HKU1.to_csv('SARS_HKU1_blosum_overlap_ex_anchors.csv')
df_SARS_229E.to_csv('SARS_229E_blosum_overlap_ex_anchors.csv')
df_SARS_NL63.to_csv('SARS_NL63_blosum_overlap_ex_anchors.csv')
df_SARS_OC43.to_csv('SARS_OC43_blosum_overlap_ex_anchors.csv')
df_SARS_EBOV.to_csv('SARS_EBOV_blosum_overlap_ex_anchors.csv')
df_SARS_H3N2.to_csv('SARS_H3N2_blosum_overlap_ex_anchors.csv')
