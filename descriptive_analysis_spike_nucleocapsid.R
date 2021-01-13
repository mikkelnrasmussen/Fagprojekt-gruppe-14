setwd('/Users/mikkel/Desktop/netMHCpan-4.1/pep2score')

library('Hmisc')

# SARS-Cov-2 spike protein 
data.SARS_spike_HKU1 <- read.csv('blosum-SARS-Cov-2_spike-HCov-HKU1/blosum_SARS-Cov-2_spike_HCov-HKU1_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_spike_229E <- read.csv('blosum-SARS-Cov-2_spike-HCov-229E/blosum_SARS-Cov-2_spike_HCov-229E_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_spike_NL63 <- read.csv('blosum-SARS-Cov-2_spike-HCov-NL63/blosum_SARS-Cov-2_spike_HCov-NL63_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_spike_OC43 <- read.csv('blosum-SARS-Cov-2_spike-HCov-OC43/blosum_SARS-Cov-2_spike_HCov-OC43_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_spike_EBOV <- read.csv('blosum-SARS-Cov-2_spike-Zaire-ebolavirus/blosum_SARS-Cov-2_spike_Zaire-ebolavirus_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_spike_H3N2 <- read.csv('blosum-SARS-Cov-2_spike-Influenza-virus-A-H3N2/blosum_SARS-Cov-2_spike_Influenza-virus-A-H3N2_combined.csv', header=TRUE, sep=",", as.is=TRUE)

# SARS-Cov-2 nucleocapsid protein 
data.SARS_nucleocapsid_HKU1 <- read.csv('blosum-SARS-Cov-2_nucleocapsid-HCov-HKU1/blosum_SARS-Cov-2_nucleocapsid_HCov-HKU1_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_nucleocapsid_229E <- read.csv('blosum-SARS-Cov-2_nucleocapsid-HCov-229E/blosum_SARS-Cov-2_nucleocapsid_HCov-229E_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_nucleocapsid_NL63 <- read.csv('blosum-SARS-Cov-2_nucleocapsid-HCov-NL63/blosum_SARS-Cov-2_nucleocapsid_HCov-NL63_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_nucleocapsid_OC43 <- read.csv('blosum-SARS-Cov-2_nucleocapsid-HCov-OC43/blosum_SARS-Cov-2_nucleocapsid_HCov-OC43_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_nucleocapsid_EBOV <- read.csv('blosum-SARS-Cov-2_nucleocapsid-Zaire-ebolavirus/blosum_SARS-Cov-2_nucleocapsid_Zaire-ebolavirus_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_nucleocapsid_H3N2 <- read.csv('blosum-SARS-Cov-2_nucleocapsid-Influenza-virus-A-H3N2/blosum_SARS-Cov-2_nucleocapsid_Influenza-virus-A-H3N2_combined.csv', header=TRUE, sep=",", as.is=TRUE)


# Histograms - spike protein
par(mfrow = c(2,2))
hist(data.SARS_spike_HKU1$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/HCov-HKU1 comparison', col = 'coral3', breaks = 250)
hist(data.SARS_spike_229E$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/HCov-229E comparison', col = 'chartreuse3', breaks = 250)
hist(data.SARS_spike_NL63$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/HCov-NL63 comparison', col = 'cyan4', breaks = 250)
hist(data.SARS_spike_OC43$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/HCov-OC43 comparison', col = 'darkorchid3', breaks = 250)
hist(data.SARS_spike_EBOV$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/EBOV comparison', col = 'black', breaks = 250, add = TRUE)
hist(data.SARS_spike_H3N2$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Spike/H3N2 comparison', col = 'white', breaks = 250, add = TRUE)

# Histograms - nucleocapsid protein
par(mfrow = c(2,2))
hist(data.SARS_nucleocapsid_HKU1$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/HCov-HKU1 comparison', col = 'coral3', breaks = 250)
hist(data.SARS_nucleocapsid_229E$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/HCov-229E comparison', col = 'chartreuse3', breaks = 250)
hist(data.SARS_nucleocapsid_NL63$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/HCov-NL63 comparison', col = 'cyan4', breaks = 250)
hist(data.SARS_nucleocapsid_OC43$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/HCov-OC43 comparison', col = 'darkorchid3', breaks = 250)
hist(data.SARS_nucleocapsid_EBOV$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/EBOV comparison', col = 'black', breaks = 250, add = TRUE)
hist(data.SARS_nucleocapsid_H3N2$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2 Nucleocapsid/H3N2 comparison', col = 'white', breaks = 250, add = TRUE)




# Boxplot - spike protein
par(mfrow = c(1,1))
boxplot(data.SARS_spike_HKU1$BLOSUM.score, data.SARS_spike_229E$BLOSUM.score, 
        data.SARS_spike_NL63$BLOSUM.score, data.SARS_spike_OC43$BLOSUM.score, 
        data.SARS_spike_EBOV$BLOSUM.score, data.SARS_spike_H3N2$BLOSUM.score, 
        main = "Boxplot of BLOSUM score", names=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", 
                                                  "Zaire-EBOV", "Influenza-A-H3N2"), 
        col=(c("coral3","chartreuse3", "cyan4", "darkorchid3", "red", "green")),
        xlab="BLOSUM SARS-Cov-2 comparison with", ylab="BLOSUM score")

# Boxplot - nucleocapsid protein
par(mfrow = c(1,1))
boxplot(data.SARS_nucleocapsid_HKU1$BLOSUM.score, data.SARS_nucleocapsid_229E$BLOSUM.score, 
        data.SARS_nucleocapsid_NL63$BLOSUM.score, data.SARS_nucleocapsid_OC43$BLOSUM.score, 
        data.SARS_nucleocapsid_EBOV$BLOSUM.score, data.SARS_nucleocapsid_H3N2$BLOSUM.score, 
        main = "Boxplot of BLOSUM score", names=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", 
                                                  "Zaire-EBOV", "Influenza-A-H3N2"), 
        col=(c("coral3","chartreuse3", "cyan4", "darkorchid3", "red", "green")),
        xlab="BLOSUM SARS-Cov-2 comparison with", ylab="BLOSUM score")


xat <- seq(0, 1, by = 0.2)
xlabels <- seq(1, 0, by=-0.2)

# Reversed cumulative distribution - spike protein

par(mfrow = c(1,1))
Ecdf(1-data.SARS_spike_HKU1$BLOSUM.score, what = 'F', col="coral3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/HCov-HKU1 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_spike_229E$BLOSUM.score, what = 'F', col="chartreuse3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/HCov-229E \n comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-229E", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("chartreuse3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_spike_NL63$BLOSUM.score, what = 'F', col="cyan4", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/HCov-NL63 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-NL63", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("cyan4","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_spike_OC43$BLOSUM.score, what = 'F', col="darkorchid3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/HCov-OC43 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


# Reversed cumulative distribution - nucleocapsid protein

par(mfrow = c(1,1))
Ecdf(1-data.SARS_nucleocapsid_HKU1$BLOSUM.score, what = 'F', col="coral3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 nucleocapsid/HCov-HKU1 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_nucleocapsid_229E$BLOSUM.score, what = 'F', col="chartreuse3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 nucleocapsid/HCov-229E \n comparison', xaxt = "n")
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-229E", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("chartreuse3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_nucleocapsid_NL63$BLOSUM.score, what = 'F', col="cyan4", main = 'Reverse cumulative distibution of \n SARS-Cov-2 nucleocapsid/HCov-NL63 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-NL63", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("cyan4","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)

Ecdf(1-data.SARS_nucleocapsid_OC43$BLOSUM.score, what = 'F', col="darkorchid3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 nucleocapsid/HCov-OC43 \n comparison', xaxt = "n")
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


# Comparison of spike and nucleocapsid
# HKU1
par(mfrow = c(1,1))
Ecdf(1-data.SARS_spike_HKU1$BLOSUM.score, what = 'F', col="coral3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/nucleocapsid \n HCov-HKU1 comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_HKU1$BLOSUM.score, what = 'F', col="coral4", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="red", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="orange", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1 (spike)", "HCov-HKU1 (nucleocapsid)", "Zaire-EBOV (spike)", "Influenza-A-H3N2 (spike)", "Zaire-EBOV (nucleocapsid)", "Influenza-A-H3N2 \n (nucleocapsid)"), 
       col=c("coral3","coral4", "black", "grey", "red", "orange"), lty = 1, cex = 0.75, text.width = 0.1)

# 229E
par(mfrow = c(1,1))
Ecdf(1-data.SARS_spike_229E$BLOSUM.score, what = 'F', col="chartreuse3", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/nucleocapsid \n HCov-229E comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_229E$BLOSUM.score, what = 'F', col="chartreuse4", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="red", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="orange", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-229E spike", "HCov-229E nucleocapsid", "Zaire-EBOV (spike)", "Influenza-A-H3N2 (spike)", "Zaire-EBOV (nucleocapsid)", "Influenza-A-H3N2 \n (nucleocapsid)"), col=c("chartreuse3","chartreuse4", "black", "grey", "red", "orange"), lty = 1, cex = 0.75, text.width = 0.1)

# NL63
par(mfrow = c(1,1))
Ecdf(1-data.SARS_spike_NL63$BLOSUM.score, what = 'F', col="cyan4", main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/nucleocapsid \n HCov-NL63 comparison', xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_NL63$BLOSUM.score, what = 'F', col="cyan3", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="red", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="orange", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-NL63 spike", "HCov-NL63 nucleocapsid", "Zaire-EBOV (spike)", "Influenza-A-H3N2 (spike)", "Zaire-EBOV (nucleocapsid)", "Influenza-A-H3N2 \n (nucleocapsid)"), col=c("cyan4","cyan3", "black", "grey", "red", "orange"), lty = 1, cex = 0.75, text.width = 0.1)

# OC43
par(mfrow = c(1,1))
Ecdf(1-data.SARS_spike_OC43$BLOSUM.score, what = 'F', col="darkorchid3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2 spike/nucleocapsid \n HCov-OC43 comparison', 
     xaxt = "n")
Ecdf(1-data.SARS_spike_EBOV$BLOSUM.score, what = 'F', col="black", add = TRUE)
Ecdf(1-data.SARS_spike_H3N2$BLOSUM.score, what = 'F', col="grey", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_OC43$BLOSUM.score, what = 'F', col="darkorchid4", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_EBOV$BLOSUM.score, what = 'F', col="red", add = TRUE)
Ecdf(1-data.SARS_nucleocapsid_H3N2$BLOSUM.score, what = 'F', col="orange", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-OC43 spike", "HCov-OC43 nucleocapsid", "Zaire-EBOV (spike)", "Influenza-A-H3N2 (spike)", "Zaire-EBOV (nucleocapsid)", "Influenza-A-H3N2 \n (nucleocapsid)"), col=c("darkorchid3","darkorchid4", "black", "grey", "red", "orange"), lty = 1, cex = 0.75, text.width = 0.1)

