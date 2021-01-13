#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 11:26:29 2021

@author: mikkel
"""

import pandas as pd
from blosumOverlap import blosumOverlap


df_SARS_spike = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_spike.csv')
df_SARS_nucleocapsid = pd.read_csv('SARS-Cov-2/SARS-Cov-2_combined_nucleocapsid.csv')
df_HKU1 = pd.read_csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv')
df_229E = pd.read_csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv')
df_NL63 = pd.read_csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv')
df_OC43 = pd.read_csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv')
df_EBOV = pd.read_csv('Zaire-ebolavirus/Zaire-ebolavirus_combined.csv')
df_H3N2 = pd.read_csv('Influenza-virus-A-H3N2/Influenza-virus-A-H3N2_combined.csv')

df_SARS_spike_SB = df_SARS_spike.loc[df_SARS_spike['EL_rank'] < 0.5]
df_SARS_nucleocapsid_SB = df_SARS_nucleocapsid.loc[df_SARS_nucleocapsid['EL_rank'] < 0.5]
df_HKU1_SB = df_HKU1.loc[df_HKU1['EL_rank'] < 0.5]
df_229E_SB = df_229E.loc[df_229E['EL_rank'] < 0.5]
df_NL63_SB = df_NL63.loc[df_NL63['EL_rank'] < 0.5]
df_OC43_SB = df_OC43.loc[df_OC43['EL_rank'] < 0.5]
df_EBOV_SB = df_EBOV.loc[df_EBOV['EL_rank'] < 0.5]
df_H3N2_SB = df_H3N2.loc[df_H3N2['EL_rank'] < 0.5]

# Spike
# =============================================================================

df_HKU1_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-HCov-HKU1/blosum_SARS-Cov-2_spike_HCov-HKU1_combined.csv')
df_229E_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-HCov-229E/blosum_SARS-Cov-2_spike_HCov-229E_combined.csv')
df_NL63_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-HCov-NL63/blosum_SARS-Cov-2_spike_HCov-NL63_combined.csv')
df_OC43_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-HCov-OC43/blosum_SARS-Cov-2_spike_HCov-OC43_combined.csv')
df_EBOV_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-Zaire-ebolavirus/blosum_SARS-Cov-2_spike_Zaire-ebolavirus_combined.csv')
df_H3N2_blosum_spike = pd.read_csv('pep2score/blosum-SARS-Cov-2_spike-Influenza-virus-A-H3N2/blosum_SARS-Cov-2_spike_Influenza-virus-A-H3N2_combined.csv')

df_SARS_HKU1_spike = blosumOverlap(df_HKU1_blosum_spike, df_SARS_spike_SB, df_HKU1_SB, 'SARS-Cov-2 spike', 'HCov-HKU1')
df_SARS_229E_spike = blosumOverlap(df_229E_blosum_spike, df_SARS_spike_SB, df_229E_SB, 'SARS-Cov-2 spike', 'HCov-229E')
df_SARS_NL63_spike = blosumOverlap(df_NL63_blosum_spike, df_SARS_spike_SB, df_NL63_SB, 'SARS-Cov-2 spike', 'HCov-NL63')
df_SARS_OC43_spike = blosumOverlap(df_OC43_blosum_spike, df_SARS_spike_SB, df_OC43_SB, 'SARS-Cov-2 spike', 'HCov-OC43')
df_SARS_EBOV_spike = blosumOverlap(df_EBOV_blosum_spike, df_SARS_spike_SB, df_EBOV_SB, 'SARS-Cov-2 spike', 'Zaire Ebolavirus')
df_SARS_H3N2_spike = blosumOverlap(df_H3N2_blosum_spike, df_SARS_spike_SB, df_H3N2_SB, 'SARS-Cov-2 spike', 'Influenza virus A H3N2')


sorted_df_SARS_spike_HKU1 = df_SARS_HKU1_spike.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_229E = df_SARS_229E_spike.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_NL63 = df_SARS_NL63_spike.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_OC43 = df_SARS_OC43_spike.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_EBOV = df_SARS_EBOV_spike.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_spike_H3N2 = df_SARS_H3N2_spike.sort_values('Proportional overlap', ascending = False)

sorted_df_SARS_spike_HKU1.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_HKU1_blosum_overlap.csv')
sorted_df_SARS_spike_229E.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_229E_blosum_overlap.csv')
sorted_df_SARS_spike_NL63.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_NL63_blosum_overlap.csv')
sorted_df_SARS_spike_OC43.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_OC43_blosum_overlap.csv')
sorted_df_SARS_spike_EBOV.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_EBOV_blosum_overlap.csv')
sorted_df_SARS_spike_H3N2.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_spike/SARS_spike_H3N2_blosum_overlap.csv')



# Nucleocapsid
# =============================================================================



df_HKU1_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-HCov-HKU1/blosum_SARS-Cov-2_nucleocapsid_HCov-HKU1_combined.csv')
df_229E_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-HCov-229E/blosum_SARS-Cov-2_nucleocapsid_HCov-229E_combined.csv')
df_NL63_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-HCov-NL63/blosum_SARS-Cov-2_nucleocapsid_HCov-NL63_combined.csv')
df_OC43_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-HCov-OC43/blosum_SARS-Cov-2_nucleocapsid_HCov-OC43_combined.csv')
df_EBOV_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-Zaire-ebolavirus/blosum_SARS-Cov-2_nucleocapsid_Zaire-ebolavirus_combined.csv')
df_H3N2_blosum_nucleocapsid = pd.read_csv('pep2score/blosum-SARS-Cov-2_nucleocapsid-Influenza-virus-A-H3N2/blosum_SARS-Cov-2_nucleocapsid_Influenza-virus-A-H3N2_combined.csv')

df_SARS_HKU1_nucleocapsid = blosumOverlap(df_HKU1_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_HKU1_SB, 'SARS-Cov-2 nucleocapsid', 'HCov-HKU1')
df_SARS_229E_nucleocapsid = blosumOverlap(df_229E_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_229E_SB, 'SARS-Cov-2 nucleocapsid', 'HCov-229E')
df_SARS_NL63_nucleocapsid = blosumOverlap(df_NL63_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_NL63_SB, 'SARS-Cov-2 nucleocapsid', 'HCov-NL63')
df_SARS_OC43_nucleocapsid = blosumOverlap(df_OC43_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_OC43_SB, 'SARS-Cov-2 nucleocapsid', 'HCov-OC43')
df_SARS_EBOV_nucleocapsid = blosumOverlap(df_EBOV_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_EBOV_SB, 'SARS-Cov-2 nucleocapsid', 'Zaire Ebolavirus')
df_SARS_H3N2_nucleocapsid = blosumOverlap(df_H3N2_blosum_nucleocapsid, df_SARS_nucleocapsid_SB, df_H3N2_SB, 'SARS-Cov-2 nucleocapsid', 'Influenza virus A H3N2')


sorted_df_SARS_nucleocapsid_HKU1 = df_SARS_HKU1_nucleocapsid.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_229E = df_SARS_229E_nucleocapsid.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_NL63 = df_SARS_NL63_nucleocapsid.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_OC43 = df_SARS_OC43_nucleocapsid.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_EBOV = df_SARS_EBOV_nucleocapsid.sort_values('Proportional overlap', ascending = False)
sorted_df_SARS_nucleocapsid_H3N2 = df_SARS_H3N2_nucleocapsid.sort_values('Proportional overlap', ascending = False)

sorted_df_SARS_nucleocapsid_HKU1.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_HKU1_blosum_overlap.csv')
sorted_df_SARS_nucleocapsid_229E.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_229E_blosum_overlap.csv')
sorted_df_SARS_nucleocapsid_NL63.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_NL63_blosum_overlap.csv')
sorted_df_SARS_nucleocapsid_OC43.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_OC43_blosum_overlap.csv')
sorted_df_SARS_nucleocapsid_EBOV.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_EBOV_blosum_overlap.csv')
sorted_df_SARS_nucleocapsid_H3N2.to_csv('Blosum-overlap-common-colds-SARS-Cov-2_nucleocapsid/SARS_nucleocapsid_H3N2_blosum_overlap.csv')