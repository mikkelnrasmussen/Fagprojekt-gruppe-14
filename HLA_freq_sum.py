#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 10:48:39 2020

@author: mikkel
"""

from collections import Counter
import pandas as pd

def HLA_freq():
    
    # Read HLA frequencies txt files from allelefrequencies.net
    with open('HLA_B_freq.txt') as f:
        hla_text = f.readlines()
    
    lines = [line.rstrip().split() for line in hla_text]
    
    # Remove excess information from file - keep HLA, freq and sample size
    for item in lines:
        del item[0:1]
        del item[1:-3]
        item.pop()
    
    # Convert freq and sample size to floats
    for item in lines:
        item[1] = float(item[1])
        item[2] = float(item[2].replace(',',''))
    
    # Create lists and store information for HLA, freq and sample size
    hla_list = []
    freq_list = []
    sample_list = []
    
    for i in range(len(lines)):
        hla_list.append(lines[i][0])
        freq_list.append(lines[i][1])
        sample_list.append(lines[i][2])

    # Create dictionary for HLA and total sample size
    hla_sum = {}
    nn = 0
    
    N = len(hla_list)
    
    # If sample size is greater than 20,000 sample size is downscaled,
    # so importance of large samples are not overestimated
    for i in range(N):
        if sample_list[i] > 20000:
            n = 20000
        else:
            n = sample_list[i]
        
        # Sum the people with a given HLA and store info in dictionary.
        # Sum the total number of people (total sample size)
        if hla_list[i][0:7] in hla_sum:
            hla_sum[hla_list[i][0:7]] += (freq_list[i] * n)
            nn += (freq_list[i] * n)
        else:
            hla_sum[hla_list[i][0:7]] = (freq_list[i] * n)
            nn += (freq_list[i] * n)
    
    # Calculate the frequency by dividing the sum of people
    # for a given HLA with the total sample size (the total number of participants)
    hla_freq_sum = {k: v / nn for k, v in hla_sum.items()}
   
    return hla_freq_sum


k = Counter(HLA_freq())
highest = k.most_common(30) 

print("Allele: Frequency")
num = 1 
for i in highest: 
    print(num, i[0],": ",round(i[1], 5)," ")
    num += 1