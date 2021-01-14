setwd('/Users/mikkel/Desktop/netMHCpan-4.1/pep2score')

library('Hmisc')

# Import csv files containing the BLOSUM-scores from the 6 comparisons with SARS-CoV-2
data.SARS_HKU1 <- read.csv('blosum-SARS-Cov-2-HCov-HKU1/blosum_SARS-Cov-2_HCov-HKU1_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_229E <- read.csv('blosum-SARS-Cov-2-HCov-229E/blosum_SARS-Cov-2_HCov-229E_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_NL63 <- read.csv('blosum-SARS-Cov-2-HCov-NL63/blosum_SARS-Cov-2_HCov-NL63_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_OC43 <- read.csv('blosum-SARS-Cov-2-HCov-OC43/blosum_SARS-Cov-2_HCov-OC43_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_EBOV <- read.csv('blosum-SARS-Cov-2-Zaire-ebolavirus/blosum_SARS-Cov-2_Zaire-ebolavirus_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.SARS_H3N2 <- read.csv('blosum-SARS-Cov-2-Influenza-virus-A-H3N2/blosum_SARS-Cov-2_Influenza-virus-A-H3N2_combined.csv', header=TRUE, sep=",", as.is=TRUE)


# Histograms of the BLOSUM-score from the 6 comparisons
par(mfrow = c(2,2))
hist(data.SARS_HKU1$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-HKU1 comparison', col = 'coral3', breaks = 250)
hist(data.SARS_229E$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-229E comparison', col = 'chartreuse3', breaks = 250)
hist(data.SARS_NL63$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-NL63 comparison', col = 'cyan4', breaks = 250)
hist(data.SARS_OC43$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/HCov-OC43 comparison', col = 'darkorchid3', breaks = 250)
hist(data.SARS_EBOV$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/EBOV comparison', col = 'black', breaks = 250, add = TRUE)
hist(data.SARS_H3N2$BLOSUM.score, xlab = 'BLOSUM score', main = 'Histogram of SARS-Cov-2/H3N2 comparison', col = 'white', breaks = 250, add = TRUE)

# Boxplot of the BLOSUM-score from the 6 comparisons
par(mfrow = c(1,1))
boxplot(data.SARS_HKU1$BLOSUM.score, data.SARS_229E$BLOSUM.score, data.SARS_NL63$BLOSUM.score, data.SARS_OC43$BLOSUM.score, data.SARS_EBOV$BLOSUM.score, data.SARS_H3N2$BLOSUM.score, 
        main = "Boxplot of BLOSUM score", names=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=(c("coral3","chartreuse3", "cyan4", "darkorchid3", "black", "grey")),
        xlab="BLOSUM SARS-Cov-2 comparison with", ylab="BLOSUM score", range=0)



# Reversed cumulative distribution - one graph for each comparison with controls
xat <- seq(0, 1, by = 0.2)
xlabels <- seq(1, 0, by=-0.2)

par(mfrow = c(1,1))
Ecdf(1-data.SARS_HKU1$BLOSUM.score, what = 'F', col="coral3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-HKU1 comparison', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_EBOV$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


Ecdf(1-data.SARS_229E$BLOSUM.score, what = 'F', col="chartreuse3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-229E comparison', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_EBOV$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-229E", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("chartreuse3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


Ecdf(1-data.SARS_NL63$BLOSUM.score, what = 'F', col="cyan4", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-NL63 comparison', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_EBOV$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-NL63", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("cyan4","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


Ecdf(1-data.SARS_OC43$BLOSUM.score, what = 'F', col="darkorchid3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2/HCov-OC43 comparison', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_EBOV$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


# Combined graph of reverse cumulative distribution curve
xat <- seq(0, 1, by = 0.2)
xlabels <- seq(1, 0, by=-0.2)

par(mfrow = c(1,1))
Ecdf(1-data.SARS_HKU1$BLOSUM.score, what = 'F', col="darkorchid3", 
     main = 'Reverse cumulative distibution of \n SARS-Cov-2 and common HCoVs \n Including anchor residues', xaxt = "n", xlab = "BLOSUM-score", ylab = "Proportion of HLA ligands")
Ecdf(1-data.SARS_229E$BLOSUM.score, what = 'F', col="chartreuse3", add = TRUE)
Ecdf(1-data.SARS_NL63$BLOSUM.score, what = 'F', col="cyan4", add = TRUE)
Ecdf(1-data.SARS_OC43$BLOSUM.score, what = 'F', col="coral3", add = TRUE)
Ecdf(1-data.SARS_EBOV$BLOSUM.score, what = 'F', xlim = c(1,0), col="black", add = TRUE)
Ecdf(1-data.SARS_H3N2$BLOSUM.score, what = 'F', xlim = c(1,0), col="grey", add = TRUE)
axis(1, at=xat, labels=xlabels)
legend("topleft", legend=c("HCov-HKU1", "HCov-229E", "HCov-NL63", "HCov-OC43", "Zaire-EBOV", "Influenza-A-H3N2"), col=c("darkorchid3", "chartreuse3", "cyan4", "coral3","black", "grey"), lty = 1, cex = 0.75, text.width = 0.1)


# Summary statictics - not used in the report
summary(data.SARS_HKU1$BLOSUM.score)
summary(data.SARS_229E$BLOSUM.score)
summary(data.SARS_NL63$BLOSUM.score)
summary(data.SARS_OC43$BLOSUM.score)

quantile(data.SARS_HKU1$BLOSUM.score, probs = seq(0, 1, 0.10), type = 2)
quantile(data.SARS_229E$BLOSUM.score, probs = seq(0, 1, 0.10), type = 2)
quantile(data.SARS_NL63$BLOSUM.score, probs = seq(0, 1, 0.10), type = 2)
quantile(data.SARS_OC43$BLOSUM.score, probs = seq(0, 1, 0.10), type = 2)


Tbl <- data.frame()

### Comparison of SARS-Cov-2 with HCov-HKU1

Tbl[1, "Antal obs."] <- sum(!is.na(data.SARS_HKU1$BLOSUM.score))
Tbl[1, "Gennemsnit"] <- mean(data.SARS_HKU1$BLOSUM.score)
Tbl[1, "Varians"] <- var(data.SARS_HKU1$BLOSUM.score)
Tbl[1, "Standardafvigelse"] <- sd(data.SARS_HKU1$BLOSUM.score)
Tbl[1, "Nedre kvartil"] <- quantile(data.SARS_HKU1$BLOSUM.score, na.rm = TRUE, type = 2)[[2]]
Tbl[1, "Median"] <- median(data.SARS_HKU1$BLOSUM.score)
Tbl[1, "Øvre kvartil"] <- quantile(data.SARS_HKU1$BLOSUM.score, na.rm = TRUE, type = 2)[[4]]

### Comparison of SARS-Cov-2 with HCov-229E

Tbl[2, "Antal obs."] <- sum(!is.na(data.SARS_229E$BLOSUM.score))
Tbl[2, "Gennemsnit"] <- mean(data.SARS_229E$BLOSUM.score)
Tbl[2, "Varians"] <- var(data.SARS_229E$BLOSUM.score)
Tbl[2, "Standardafvigelse"] <- sd(data.SARS_229E$BLOSUM.score)
Tbl[2, "Nedre kvartil"] <- quantile(data.SARS_229E$BLOSUM.score, na.rm = TRUE, type = 2)[[2]]
Tbl[2, "Median"] <- median(data.SARS_229E$BLOSUM.score)
Tbl[2, "Øvre kvartil"] <- quantile(data.SARS_229E$BLOSUM.score, na.rm = TRUE, type = 2)[[4]]

### Comparison of SARS-Cov-2 with HCov-NL63

Tbl[3, "Antal obs."] <- sum(!is.na(data.SARS_NL63$BLOSUM.score))
Tbl[3, "Gennemsnit"] <- mean(data.SARS_NL63$BLOSUM.score)
Tbl[3, "Varians"] <- var(data.SARS_NL63$BLOSUM.score)
Tbl[3, "Standardafvigelse"] <- sd(data.SARS_NL63$BLOSUM.score)
Tbl[3, "Nedre kvartil"] <- quantile(data.SARS_NL63$BLOSUM.score, na.rm = TRUE, type = 2)[[2]]
Tbl[3, "Median"] <- median(data.SARS_NL63$BLOSUM.score)
Tbl[3, "Øvre kvartil"] <- quantile(data.SARS_NL63$BLOSUM.score, na.rm = TRUE, type = 2)[[4]]

### Comparison of SARS-Cov-2 with HCov-OC43

Tbl[4, "Antal obs."] <- sum(!is.na(data.SARS_OC43$BLOSUM.score))
Tbl[4, "Gennemsnit"] <- mean(data.SARS_OC43$BLOSUM.score)
Tbl[4, "Varians"] <- var(data.SARS_OC43$BLOSUM.score)
Tbl[4, "Standardafvigelse"] <- sd(data.SARS_OC43$BLOSUM.score)
Tbl[4, "Nedre kvartil"] <- quantile(data.SARS_OC43$BLOSUM.score, na.rm = TRUE, type = 2)[[2]]
Tbl[4, "Median"] <- median(data.SARS_OC43$BLOSUM.score)
Tbl[4, "Øvre kvartil"] <- quantile(data.SARS_OC43$BLOSUM.score, na.rm = TRUE, type = 2)[[4]]


row.names(Tbl) <- c("SARS-Cov-2/HCov-HKU1", "SARS-Cov-2/HCov-229E","SARS-Cov-2/HCov-NL63", "SARS-Cov-2/HCov-OC43")

Tbl


