setwd('/Users/mikkel/Desktop/netMHCpan-4.1/')
library(ggplot2)
library(tidyr)
library(reshape)

data.SARS <- read.csv('SARS-Cov-2/SARS-Cov-2_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.HKU1 <- read.csv('Human-coronavirus-HKU1/Human-coronavirus-HKU1_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.229E <- read.csv('Human-coronavirus-229E/Human-coronavirus-229E_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.NL63 <- read.csv('Human-coronavirus-NL63/Human-coronavirus-NL63_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.OC43 <- read.csv('Human-coronavirus-OC43/Human-coronavirus-OC43_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.EBOV <- read.csv('Zaire-ebolavirus/Zaire-ebolavirus_combined.csv', header=TRUE, sep=",", as.is=TRUE)
data.H3N2 <- read.csv('Influenza-virus-A-H3N2/Influenza-virus-A-H3N2_combined.csv', header=TRUE, sep=",", as.is=TRUE)


data.SARS_SB <- subset(data.SARS, EL_rank < 0.5)
data.HKU1_SB <- subset(data.HKU1, EL_rank < 0.5)
data.229E_SB <- subset(data.229E, EL_rank < 0.5)
data.NL63_SB <- subset(data.NL63, EL_rank < 0.5)
data.OC43_SB <- subset(data.OC43, EL_rank < 0.5)
data.EBOV_SB <- subset(data.EBOV, EL_rank < 0.5)
data.H3N2_SB <- subset(data.H3N2, EL_rank < 0.5)



# Creating barplots using ggplot2 package - SARS-CoV-2
viruses_SARS <- c('SARS-CoV-2')

names_SARS <- c(names(table(data.SARS_SB$HLA)))
data.hlaCombined_SARS <- data.frame(as.vector(table(data.SARS_SB$HLA)))

data_SARS=data.frame(names,cbind(data.hlaCombined_SARS))
colnames(data_SARS) <- c("names", viruses_SARS)
data_SARS.m <- melt(data_SARS, id.vars='names')
ggplot(data_SARS.m, aes(names, value)) +   
        geom_bar(aes(fill = variable), position = position_dodge(width=0.8), stat="identity", width = 0.7) + theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
        ggtitle("Number of HLA-ligands per HLA allele") + labs(y="Number of HLA-ligands", x = "HLA alleles") + labs(colour = "Viruses") +
        scale_fill_manual("Viruses", values = c("SARS-CoV-2" = "red"))


# Creating barplots using ggplot2 package - HCoVs
viruses <- c('SARS-CoV-2', "HCoV-HKU1", "HCoV-229E", "HCoV-NL63", "HCoV-OC43")

names <- c(names(table(data.SARS_SB$HLA)))
data.hlaCombined <- data.frame(as.vector(table(data.SARS_SB$HLA)),
                               as.vector(table(data.HKU1_SB$HLA)), 
                               as.vector(table(data.229E_SB$HLA)), 
                               as.vector(table(data.NL63_SB$HLA)),
                               as.vector(table(data.OC43_SB$HLA)))
                                
data=data.frame(names, cbind(data.hlaCombined))
colnames(data) <- c("names", viruses)
data.m <- melt(data, id.vars='names')
ggplot(data.m, aes(names, value)) +   
        geom_bar(aes(fill = variable), position = position_dodge(width=0.8), stat="identity", width = 0.7) + theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
        ggtitle("Number of HLA-ligands per HLA allele") + labs(y="Number of HLA-ligands", x = "HLA alleles") + labs(colour = "Viruses") +
        scale_fill_manual("Viruses", values = c("SARS-CoV-2" = "red", "HCoV-HKU1" = "coral", "HCoV-229E" = "chartreuse3", "HCoV-NL63" = "cyan3", "HCoV-OC43" = "darkorchid3"))




# Creating barplots using ggplot2 package - Controls, EBOV and Influenza-A H3N2
control_viruses <- c("Zaire-EBOV", "Influenza-A-H3N2")

names <- c(names(table(data.EBOV_SB$HLA)))
data.hlaCombined_controls <- data.frame(as.vector(table(data.EBOV_SB$HLA)), 
                               as.vector(table(data.H3N2_SB$HLA)))

data_controls=data.frame(names, cbind(data.hlaCombined_controls))
colnames(data_controls) <- c("names", control_viruses)
data_controls.m <- melt(data_controls, id.vars='names')
ggplot(data_controls.m, aes(names, value)) +   
        geom_bar(aes(fill = variable), position = position_dodge(width=0.8), stat="identity", width = 0.7) + theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
        ggtitle("Number of HLA-ligands per HLA allele") + labs(y="Number of HLA-ligands", x = "HLA alleles") + labs(colour = "Viruses") +
        scale_fill_manual("Viruses", values = c("Zaire-EBOV" = "black", "Influenza-A-H3N2" = "grey"))





