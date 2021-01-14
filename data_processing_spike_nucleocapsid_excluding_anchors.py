#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:20:13 2020

@author: mikkel
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from overlap_excluding_anchors import overlap_excluding


# Import of csv files containing predicted HLA-ligands for the spike and nucleocapsid of SARS-CoV-2, the four HCoVs and two controls; 
# SARS-Cov-2, HCov-229E, HCov-HKU1, HCov-NL63 and HCov-OC43
df_SARS_spike = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_spike.csv')
df_SARS_nucleocapsid = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_nucleocapsid.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')

# Restrict data to only HLA-ligands (strong binders (SB)) with EL_rank < 0.5
df_SARS_spike_EL = df_SARS_spike.loc[df_SARS_spike['EL_rank'] < 0.5]
df_SARS_nucleocapsid_EL = df_SARS_nucleocapsid.loc[df_SARS_nucleocapsid['EL_rank'] < 0.5]
df_229E_EL = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_HKU1_EL = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_NL63_EL = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_EL = df_OC43.loc[df_OC43['EL_rank'] < 0.5]


# Spike
# =============================================================================

# Calculate the actual overlap, proportional overlap and the number HLA-ligands 
# of the SARS-CoV-2 spike protein and a compared virus for a given HLA, where anchor points has been excluded
df_SARS_spike_229E = overlap_excluding(df_SARS_spike_EL, df_229E_EL, 'SARS-Cov-2 spike', 'HCov-229E')
df_SARS_spike_HKU1 = overlap_excluding(df_SARS_spike_EL, df_HKU1_EL, 'SARS-Cov-2 spike', 'HCov-HKU1')
df_SARS_spike_NL63 = overlap_excluding(df_SARS_spike_EL, df_NL63_EL, 'SARS-Cov-2 spike', 'HCov-NL63')
df_SARS_spike_OC43 = overlap_excluding(df_SARS_spike_EL, df_OC43_EL, 'SARS-Cov-2 spike', 'HCov-OC43')

# Rank by proportional overlap
sorted_df_SARS_spike_229E = df_SARS_spike_229E.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_HKU1 = df_SARS_spike_HKU1.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_NL63 = df_SARS_spike_NL63.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_OC43 = df_SARS_spike_OC43.sort_values('Proportional overlap', ascending = False)

# Export as csv file
sorted_df_SARS_spike_229E.to_csv('SARS_spike_229E_rank_excluding_anchors.csv')
sorted_df_SARS_spike_HKU1.to_csv('SARS_spike_HKU1_rank_excluding_anchors.csv')
sorted_df_SARS_spike_NL63.to_csv('SARS_spike_NL63_rank_excluding_anchors.csv')
sorted_df_SARS_spike_OC43.to_csv('SARS_spike_OC43_rank_excluding_anchors.csv')

# Nucleocapsid
# =============================================================================

# Calculate the actual overlap, proportional overlap and the number HLA-ligands 
# of the SARS-CoV-2 nucleocapsid protein and a compared virus for a given HLA, where anchor points has been excluded
df_SARS_nucleocapsid_229E = overlap_excluding(df_SARS_nucleocapsid_EL, df_229E_EL, 'SARS-Cov-2 nucleocapsid', 'HCov-229E')
df_SARS_nucleocapsid_HKU1 = overlap_excluding(df_SARS_nucleocapsid_EL, df_HKU1_EL, 'SARS-Cov-2 nucleocapsid', 'HCov-HKU1')
df_SARS_nucleocapsid_NL63 = overlap_excluding(df_SARS_nucleocapsid_EL, df_NL63_EL, 'SARS-Cov-2 nucleocapsid', 'HCov-NL63')
df_SARS_nucleocapsid_OC43 = overlap_excluding(df_SARS_nucleocapsid_EL, df_OC43_EL, 'SARS-Cov-2 nucleocapsid', 'HCov-OC43')

# Rank by proportional overlap
sorted_df_SARS_nucleocapsid_229E = df_SARS_nucleocapsid_229E.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_HKU1 = df_SARS_nucleocapsid_HKU1.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_NL63 = df_SARS_nucleocapsid_NL63.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_OC43 = df_SARS_nucleocapsid_OC43.sort_values('Proportional overlap', ascending = False)

# Export as csv file
sorted_df_SARS_nucleocapsid_229E.to_csv('SARS_nucleocapsid_229E_rank_excluding_anchors.csv')
sorted_df_SARS_nucleocapsid_HKU1.to_csv('SARS_nucleocapsid_HKU1_rank_excluding_anchors.csv')
sorted_df_SARS_nucleocapsid_NL63.to_csv('SARS_nucleocapsid_NL63_rank_excluding_anchors.csv')
sorted_df_SARS_nucleocapsid_OC43.to_csv('SARS_nucleocapsid_OC43_rank_excluding_anchors.csv')
