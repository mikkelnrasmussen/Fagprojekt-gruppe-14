#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 16:24:02 2020

@author: mikkel
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from overlap_excluding_anchors import overlap_excluding


# Import of csv files containing predicted HLA-ligand from NetMHCpan for the five viruses; SARS-Cov-2, HCov-229E, HCov-HKU1, HCov-NL63 and HCov-OC43.
df_SARS = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')
df_EBOV = pd.read_csv('Zaire-ebolavirus/Zaire-ebolavirus_combined.csv')
df_H3N2 = pd.read_csv('Influenza-virus-A-H3N2/Influenza-virus-A-H3N2_combined.csv')

# Restrict data to only HLA-ligands (strong binders (SB)), EL_rank < 0.5
df_SARS_SB = df_SARS.loc[df_SARS['EL_rank'] < 0.5]
df_EBOV_SB = df_EBOV.loc[df_EBOV['EL_rank'] < 0.5]
df_H3N2_SB = df_H3N2.loc[df_H3N2['EL_rank'] < 0.5] 
df_229E_EL = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_HKU1_EL = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_NL63_EL = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_EL = df_OC43.loc[df_OC43['EL_rank'] < 0.5]

# Calculate the actual overlap, proportional overlap and the number HLA-ligands 
# in the two compared viruses for a given HLA, where anchor points have been excluded
df_SARS_EBOV = overlap_excluding(df_SARS_SB, df_EBOV_SB, 'SARS-Cov-2', 'Zaire Ebolavirus')
df_SARS_H3N2 = overlap_excluding(df_SARS_SB, df_H3N2_SB, 'SARS-Cov-2', 'Influenza virus A H3N2')
df_SARS_229E = overlap_excluding(df_SARS_EL, df_229E_EL, 'SARS-Cov-2', 'HCov-229E')
df_SARS_HKU1 = overlap_excluding(df_SARS_EL, df_HKU1_EL, 'SARS-Cov-2', 'HCov-HKU1')
df_SARS_NL63 = overlap_excluding(df_SARS_EL, df_NL63_EL, 'SARS-Cov-2', 'HCov-NL63')
df_SARS_OC43 = overlap_excluding(df_SARS_EL, df_OC43_EL, 'SARS-Cov-2', 'HCov-OC43')

# Rank by proportional overlap
sorted_df_SARS_EBOV = df_SARS_EBOV.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_H3N2 = df_SARS_H3N2.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_229E = df_SARS_229E.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_HKU1 = df_SARS_HKU1.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_NL63 = df_SARS_NL63.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_OC43 = df_SARS_OC43.sort_values('Proportional overlap', ascending = False)

# Export as csv file
sorted_df_SARS_EBOV.to_csv('SARS_EBOV_rank_excluding_anchors.csv')
sorted_df_SARS_H3N2.to_csv('SARS_H3N2_rank_excluding_anchors.csv')
sorted_df_SARS_229E.to_csv('SARS_229E_rank_excluding_anchors.csv')
sorted_df_SARS_HKU1.to_csv('SARS_HKU1_rank_excluding_anchors.csv')
sorted_df_SARS_NL63.to_csv('SARS_NL63_rank_excluding_anchors.csv')
sorted_df_SARS_OC43.to_csv('SARS_OC43_rank_excluding_anchors.csv')
