import pandas as pd
import os
import numpy as np
from collections import Counter
#loading of IEDB data
dfinf = pd.read_csv('influenza.csv')

dfh = pd.DataFrame(["HLA-B*35:03","HLA-B*35:03","HLA-C*07:02","HLA-A*31:01"
                   ,"HLA-B*27:05","HLA-B*27:06","HLA-B*27:05","HLA-B*27:06"
                   ,"HLA-B*07:02","HLA-C*07:02","HLA-A*03:01","HLA-A*11:01"
                   ,"HLA-A*01:01","HLA-A*01:01"], columns = ['MHC'])

dfhealthy = pd.read_csv('healthy_donors_netMHCpan_predicted_HLA.csv', sep = ";")

dfall = pd.read_csv('HLAAbias_SARSCOV2.csv',sep = ';')

dfcovid = pd.read_csv('HLACOVID19.csv',sep =';')
#frequency analysis function
def freqHLA(dataframe):
#deleting everything besides HLA alleles
    uql = pd.unique(dataframe.MHC)
    uql = np.delete(uql,np.where(uql == 'HLA class I'))
    uql = np.delete(uql,np.where(uql == 'HLA class II'))
    uql = np.delete(uql,np.where(uql == 'H2.kd'))

    
#creating empty list
    a = {}
#for-loop calculating the frequency  
    for i in uql:
        temp_hla = dataframe.loc[dataframe['MHC'] == i]
        a[i]= len(temp_hla)/len(dataframe)*100
#creating a counter        
    b = Counter(a)
    b = b.most_common(38)
    
    num = 1 
    for i in b: 
        print(num, i[0],": ",round(i[1], 5)," ")
        num += 1

#different dataframes from IEDB (alldonors, COVID patients, healthy patients etc.)
covidex = dfcovid[dfcovid.MHC != 'HLA-A*02:01']
exdf = dfinf[dfinf['MHC'].map(len) == 11]
exdf_2 = dfinf.loc[dfinf['MHC'].map(len) == 11]
exdf_2 = exdf_2.loc[dfinf['MHC'] != 'HLA class I']         
ndf = exdf_2[exdf_2.MHC != 'HLA-A*02:01']