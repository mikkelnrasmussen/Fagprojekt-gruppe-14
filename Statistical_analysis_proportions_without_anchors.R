# CI for the highest proportional overlap obtained from the analysis of identical HLA-ligands
# excluding anchor points, for the four common HCoVs.

# HLA-B44:02
p_hat_229E <- 317/27986
(p_hat_229E + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_HKU1 * (1 - p_hat_HKU1))/27986))*100
prop.test(x = c(317, 45), n = c(27986, 5436), correct = FALSE)

# HLA-B44:03
p_hat_HKU1 <- 818/30520
(p_hat_HKU1 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_HKU1 * (1 - p_hat_HKU1))/30520))*100
prop.test(x = c(818, 90), n = c(30520, 5205), correct = FALSE)

# HLA-B44:03
p_hat_NL63 <- 385/28619
(p_hat_NL63 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_NL63 * (1 - p_hat_NL63))/28619))*100
prop.test(x = c(51, 385), n = c(5759, 28619), correct = FALSE)

# HLA-B44:03
p_hat_OC43 <- 813/30520
(p_hat_OC43 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_OC43 * (1 - p_hat_OC43))/30520))*100
prop.test(x = c(813, 84), n = c(30520, 5205), correct = FALSE)


# CI and comparison of the lowest value found in the overlap analysis without anchors
p_hat_229E_low_ex <- 45/5436
(p_hat_229E_low_ex + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_229E_low_ex * (1 - p_hat_229E_low_ex))/5436))*100

# Proportion test of the lowest proportional overlap in the overlap analysis without anchors (HCov-229E HLA-A29:02)
# with the highest proportion in the control viruses
prop.test(x = c(45, 1), n = c(5436, 2583), correct = FALSE)


