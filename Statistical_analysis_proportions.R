# CI for the highest proportional overlap obtained from the analysis of identical HLA-ligands
# including anchor points, for the four common HCoVs.

p_hat_229E <- 118/22196
(p_hat_229E + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_HKU1 * (1 - p_hat_HKU1))/22196))*100
prop.test(x = c(118, 21), n = c(22196, 5436), correct = FALSE)

p_hat_HKU1 <- 423/26482
(p_hat_HKU1 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_HKU1 * (1 - p_hat_HKU1))/26482))*100
prop.test(x = c(423, 56), n = c(26482, 5205), correct = FALSE)

p_hat_NL63 <- 199/28588
(p_hat_NL63 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_NL63 * (1 - p_hat_NL63))/28588))*100
prop.test(x = c(30, 199), n = c(5837, 28588), correct = FALSE)


p_hat_OC43 <- 478/30431
(p_hat_OC43 + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_OC43 * (1 - p_hat_OC43))/30431))*100
prop.test(x = c(478, 54), n = c(30431, 5205), correct = FALSE)


# CI and comparison of the lowest value found in the overlap analysis with anchors (HCov-229E HLA-A29:02)
p_hat_229E_low <- 21/5436
(p_hat_229E_low + c(-1, 1) * qnorm(0.975) * sqrt((p_hat_229E_low * (1 - p_hat_229E_low))/5436))*100

# Proportion test of the lowest proportional overlap in the overlap analysis with anchors (HCov-229E HLA-A29:02)
prop.test(x = c(21, 0), n = c(5436, 2216), correct = FALSE)
