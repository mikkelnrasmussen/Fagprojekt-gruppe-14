#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:38:38 2020

@author: mikkel
"""

import os
import subprocess

os.chdir("/Users/mikkel/Desktop/netMHCpan-4.1/pep2score")


def blosum(virus1, virus2):
    
    hla_list = ['HLA-A02-01', 'HLA-A01-01', 'HLA-A03-01', 'HLA-A24-02', 'HLA-A11-01', 'HLA-A26-01', 
            'HLA-A32-01', 'HLA-A68-01', 'HLA-A25-01', 'HLA-A31-01', 'HLA-A29-02', 'HLA-A23-01', 
            'HLA-B07-02', 'HLA-B08-01', 'HLA-B15-01', 'HLA-B51-01', 'HLA-B44-02', 'HLA-B18-01', 
            'HLA-B35-01', 'HLA-B44-03', 'HLA-B40-01', 'HLA-B13-02', 'HLA-B27-05', 'HLA-B57-01', 
            'HLA-B35-03', 'HLA-B38-01', 'HLA-B58-01', 'HLA-C07-01', 'HLA-C04-01', 'HLA-C07-02', 
            'HLA-C06-02', 'HLA-C12-03', 'HLA-C05-01', 'HLA-C02-02', 'HLA-C03-04', 'HLA-C03-03', 
            'HLA-C01-02', 'HLA-C15-02']

    length_list = ["8-mer", "9-mer", "10-mer", "11-mer"]
    
    for i in hla_list:
        for j in length_list:
            
            command = "./pep2score_db -t 1 -blf ./BLOSUM62" + " " + virus1 + "/" + virus1 + "_" + i + "_" + j + ".csv" + " "  + virus2 + "/" + virus2 + "_"  +  i + "_" + j + ".csv" + " " + ">" + "blosum-" + virus1 + "-" + virus2 + "/" + "blosum" + virus1 + "_" + virus2 + "_" + i + "_" + j + ".txt"
        
            os.system(command)
            
            

blosum("SARS-Cov-2", "HCov-HKU1")
blosum("SARS-Cov-2", "HCov-229E")
blosum("SARS-Cov-2", "HCov-NL63")
blosum("SARS-Cov-2", "HCov-OC43")
blosum("SARS-Cov-2", 'Zaire-ebolavirus')
blosum("SARS-Cov-2", 'Influenza-virus-A-H3N2')

blosum("SARS-Cov-2_spike", "HCov-HKU1")
blosum("SARS-Cov-2_spike", "HCov-229E")
blosum("SARS-Cov-2_spike", "HCov-NL63")
blosum("SARS-Cov-2_spike", "HCov-OC43")
blosum("SARS-Cov-2_spike", 'Zaire-ebolavirus')
blosum("SARS-Cov-2_spike", 'Influenza-virus-A-H3N2')

blosum("SARS-Cov-2_nucleocapsid", "HCov-HKU1")
blosum("SARS-Cov-2_nucleocapsid", "HCov-229E")
blosum("SARS-Cov-2_nucleocapsid", "HCov-NL63")
blosum("SARS-Cov-2_nucleocapsid", "HCov-OC43")
blosum("SARS-Cov-2_nucleocapsid", 'Zaire-ebolavirus')
blosum("SARS-Cov-2_nucleocapsid", 'Influenza-virus-A-H3N2')