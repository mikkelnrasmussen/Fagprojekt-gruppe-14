setwd('/Users/mikkel/Desktop/netMHCpan-4.1/pep2score')

library('Hmisc')

data.SARS_HKU1_ex <- read.csv('blosum-SARS-Cov-2-HCov-HKU1-ex-anchors/blosum_SARS-Cov-2_HCov-HKU1_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E_ex <- read.csv('blosum-SARS-Cov-2-HCov-229E-ex-anchors/blosum_SARS-Cov-2_HCov-229E_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63_ex <- read.csv('blosum-SARS-Cov-2-HCov-NL63-ex-anchors/blosum_SARS-Cov-2_HCov-NL63_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43_ex <- read.csv('blosum-SARS-Cov-2-HCov-OC43-ex-anchors/blosum_SARS-Cov-2_HCov-OC43_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_EBOV_ex <- read.csv('blosum-SARS-Cov-2-Zaire-ebolavirus-ex-anchors/blosum_SARS-Cov-2_Zaire-ebolavirus_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_H3N2_ex <- read.csv('blosum-SARS-Cov-2-Influenza-virus-A-H3N2-ex-anchors/blosum_SARS-Cov-2_Influenza-virus-A-H3N2_combined_ex_anchors.csv', header=TRUE, sep=",", as.is=TRUE)

# Histograms

par(mfrow = c(2,2))
hist(data.SARS_HKU1_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-HKU1 comparison \n - excluding anchors', col = 'coral3', breaks = 250)
hist(data.SARS_229E_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-229E comparison\n - excluding anchors', col = 'chartreuse3', breaks = 250)
hist(data.SARS_NL63_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-NL63 comparison\n - excluding anchors', col = 'cyan4', breaks = 250)
hist(data.SARS_OC43_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-OC43 comparison\n - excluding anchors', col = 'darkorchid3', breaks = 250)
hist(data.SARS_EBOV_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/EBOV comparison\n - excluding anchors', col = 'black', breaks = 250, add = TRUE)
hist(data.SARS_H3N2_ex$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/H3N2 comparison\n - excluding anchors', col = 'white', breaks = 250, add = TRUE)

# Boxplots

par(mfrow = c(1,1))
boxplot(data.SARS_HKU1_ex$BLOSUM.score, data.SARS_229E_ex$BLOSUM.score, data.SARS_NL63_ex$BLOSUM.score, data.SARS_OC43_ex$BLOSUM.score, data.SARS_EBOV_ex$BLOSUM.score, data.SARS_H3N2_ex$BLOSUM.score, 
        main = "Boxplot of BLOSUM score\n - excluding anchors", names=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=(c("coral3","chartreuse3", "cyan4", "darkorchid3", "black", "white")),
        xlab="BLOSUM SARS-Cov-2 comparison with", ylab="BLOSUM score", range = 0)


# Reversed cumulative distribution

xat <- seq(0, 1, by = 0.2)
xlabels <- seq(1, 0, by=-0.2)

par(mfrow = c(1,1))
Ecdf(data.SARS_HKU1_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="coral3", main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-HKU1 comparison', xaxt = 'n')
Ecdf(data.SARS_EBOV_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(data.SARS_H3N2_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(data.SARS_229E_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="chartreuse3", main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-229E comparison', xaxt = 'n')
Ecdf(data.SARS_EBOV_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(data.SARS_H3N2_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-229E", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("chartreuse3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(data.SARS_NL63_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="cyan4", main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-NL63 comparison', xaxt = 'n')
Ecdf(data.SARS_EBOV_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(data.SARS_H3N2_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-NL63", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("cyan4","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(data.SARS_OC43_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="darkorchid3", main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-OC43 comparison', xaxt = 'n')
Ecdf(data.SARS_EBOV_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(data.SARS_H3N2_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


# Combined reversed cumulative distribution 

xat <- seq(0, 1, by = 0.2)
xlabels <- seq(1, 0, by=-0.2)

par(mfrow = c(1,1))
Ecdf(1-data.SARS_HKU1_ex$BLOSUM.score, what = 'F',col="darkorchid3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2 and common HCoVs \n Excluding anchor residues', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_229E_ex$BLOSUM.score, what = 'F', col="chartreuse3", add = TRUE)
Ecdf(1-data.SARS_NL63_ex$BLOSUM.score, what = 'F', col="cyan4", add = TRUE)
Ecdf(1-data.SARS_OC43_ex$BLOSUM.score, what = 'F', col="coral3", add = TRUE)
Ecdf(1-data.SARS_EBOV_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2_ex$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3", "chartreuse3", "cyan4", "coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)
