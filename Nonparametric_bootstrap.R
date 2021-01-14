setwd('/Users/mikkel/Desktop/netMHCpan-4.1')

# Load data from identical HLA-ligand analysis, including and excluding, anchor points
# and load data from 0.8 thershold blosum analysis 

data.SARS_HKU1_rank <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2_ranked/SARS_HKU1_rank.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E_rank <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2_ranked/SARS_229E_rank.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63_rank <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2_ranked/SARS_NL63_rank.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43_rank <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2_ranked/SARS_OC43_rank.csv', 
                                header=TRUE, sep=",", as.is=TRUE)

data.SARS_HKU1_rank_ex <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2-ranked-excluding-anchors/SARS_HKU1_rank_excluding_anchors.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E_rank_ex <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2-ranked-excluding-anchors/SARS_229E_rank_excluding_anchors.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63_rank_ex <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2-ranked-excluding-anchors/SARS_NL63_rank_excluding_anchors.csv', 
                                header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43_rank_ex <- read.csv('Overlap-of-common-colds-and-SARS-Cov-2-ranked-excluding-anchors/SARS_OC43_rank_excluding_anchors.csv', 
                                header=TRUE, sep=",", as.is=TRUE)

data.SARS_HKU1_rank_blosum <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2/SARS_HKU1_blosum_overlap.csv', 
                                   header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E_rank_blosum <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2/SARS_229E_blosum_overlap.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63_rank_blosum <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2/SARS_NL63_blosum_overlap.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43_rank_blosum <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2/SARS_OC43_blosum_overlap.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)

data.SARS_HKU1_rank_blosum_ex <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2-ex-anchors/SARS_HKU1_blosum_overlap_ex_anchors.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E_rank_blosum_ex <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2-ex-anchors/SARS_229E_blosum_overlap_ex_anchors.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63_rank_blosum_ex <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2-ex-anchors/SARS_NL63_blosum_overlap_ex_anchors.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43_rank_blosum_ex <- read.csv('Blosum-overlap-common-colds-SARS-Cov-2-ex-anchors/SARS_OC43_blosum_overlap_ex_anchors.csv', 
                                       header=TRUE, sep=",", as.is=TRUE)

# Assign "HLA" as column name
names(data.SARS_HKU1_rank)[names(data.SARS_HKU1_rank) == "X"] <- "HLA"
names(data.SARS_229E_rank)[names(data.SARS_229E_rank) == "X"] <- "HLA"
names(data.SARS_NL63_rank)[names(data.SARS_NL63_rank) == "X"] <- "HLA"
names(data.SARS_OC43_rank)[names(data.SARS_OC43_rank) == "X"] <- "HLA"

names(data.SARS_HKU1_rank_ex)[names(data.SARS_HKU1_rank_ex) == "X"] <- "HLA"
names(data.SARS_229E_rank_ex)[names(data.SARS_229E_rank_ex) == "X"] <- "HLA"
names(data.SARS_NL63_rank_ex)[names(data.SARS_NL63_rank_ex) == "X"] <- "HLA"
names(data.SARS_OC43_rank_ex)[names(data.SARS_OC43_rank_ex) == "X"] <- "HLA"

names(data.SARS_HKU1_rank_blosum)[names(data.SARS_HKU1_rank_blosum) == "X"] <- "HLA"
names(data.SARS_229E_rank_blosum)[names(data.SARS_229E_rank_blosum) == "X"] <- "HLA"
names(data.SARS_NL63_rank_blosum)[names(data.SARS_NL63_rank_blosum) == "X"] <- "HLA"
names(data.SARS_OC43_rank_blosum)[names(data.SARS_OC43_rank_blosum) == "X"] <- "HLA"

names(data.SARS_HKU1_rank_blosum_ex)[names(data.SARS_HKU1_rank_blosum_ex) == "X"] <- "HLA"
names(data.SARS_229E_rank_blosum_ex)[names(data.SARS_229E_rank_blosum_ex) == "X"] <- "HLA"
names(data.SARS_NL63_rank_blosum_ex)[names(data.SARS_NL63_rank_blosum_ex) == "X"] <- "HLA"
names(data.SARS_OC43_rank_blosum_ex)[names(data.SARS_OC43_rank_blosum_ex) == "X"] <- "HLA"

# Function calculating the rank score median
rankScoreMedian <- function(hla_list){
        
        hla_a <- c()
        hla_b <- c()
        hla_c <- c()
        
        rank <- 1
        for (i in hla_list) {
                hla <- substr(i, 0, 5)
                if (hla == 'HLA-A'){
                        hla_a <- c(hla_a, rank)
                        
                } else if (hla == 'HLA-B'){
                        hla_b <- c(hla_b, rank)
                
                } else if (hla == 'HLA-C'){
                        hla_c <- c(hla_c, rank)
                }
                
                rank <- rank + 1
                        
        }
        
        return(c(median(hla_a), median(hla_b), median(hla_c)))
}

hla_list <- c("HLA-B", "HLA-B", "HLA-C", "HLA-C", "HLA-C", "HLA-B", "HLA-C", "HLA-B", "HLA-C", "HLA-C", 
              "HLA-A", "HLA-C", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", 
              "HLA-A", "HLA-B", "HLA-C", "HLA-C", "HLA-C", "HLA-A", "HLA-A", "HLA-B", "HLA-A", "HLA-A", 
              "HLA-A", "HLA-C", "HLA-A", "HLA-A", "HLA-B", "HLA-A", "HLA-A", "HLA-A")

# Perform the non-parametric boostrap with k = 100,000, replace sat to false, 
# as the we need the same HLA allele cannot occur more than once in the ranked HLA list
k <- 100000
simsamples <- replicate(k, sample(hla_list, replace = FALSE))
simmedians <- apply(simsamples, 2, rankScoreMedian)
row.names(simmedians) <- c("HLA-A", "HLA-B", "HLA-C")

# Confidence intervals for the rank score medians for the analysis with anchor points 
quantile(simmedians['HLA-A', ], c(0.025,0.975))
quantile(simmedians['HLA-B', ], c(0.025,0.975))
quantile(simmedians['HLA-C', ], c(0.025,0.975))

rankScoreMedian(data.SARS_HKU1_rank[['HLA']])
rankScoreMedian(data.SARS_229E_rank[['HLA']])
rankScoreMedian(data.SARS_NL63_rank[['HLA']])
rankScoreMedian(data.SARS_OC43_rank[['HLA']])

# P-values with anchor points
sum(simmedians['HLA-A', ] > 28.5) / k
sum(simmedians['HLA-B', ] < 16) / k
sum(simmedians['HLA-C', ] < 17) / k

# Confidence intervals for the rank score medians for the analysis without anchor points 
rankScoreMedian(data.SARS_HKU1_rank_ex[['HLA']])
rankScoreMedian(data.SARS_229E_rank_ex[['HLA']])
rankScoreMedian(data.SARS_NL63_rank_ex[['HLA']])
rankScoreMedian(data.SARS_OC43_rank_ex[['HLA']])

# P-values without anchor points
sum(simmedians['HLA-A', ] > 31) / k
sum(simmedians['HLA-B',] < 13) / k

# Confidence intervals for the rank score medians after BLOSUM for the analysis with anchor points 
rankScoreMedian(data.SARS_HKU1_rank_blosum[['HLA']])
rankScoreMedian(data.SARS_229E_rank_blosum[['HLA']])
rankScoreMedian(data.SARS_NL63_rank_blosum[['HLA']])
rankScoreMedian(data.SARS_OC43_rank_blosum[['HLA']])

# P-values, after BLOSUM with anchor points
sum(simmedians['HLA-A', ] > 29.5) / k
sum(simmedians['HLA-B',] < 12) / k

# Confidence intervals for the rank score medians after BLOSUM for the analysis without anchor points 
rankScoreMedian(data.SARS_HKU1_rank_blosum_ex[['HLA']])
rankScoreMedian(data.SARS_229E_rank_blosum_ex[['HLA']])
rankScoreMedian(data.SARS_NL63_rank_blosum_ex[['HLA']])
rankScoreMedian(data.SARS_OC43_rank_blosum_ex[['HLA']])

# P-values, after BLOSUM without anchor points
sum(simmedians['HLA-A', ] > 32) / k
sum(simmedians['HLA-B', ] < 13) / k




