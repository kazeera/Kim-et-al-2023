#' Purpose:
#' To compute paired wilcoxin significance testing, which is a two tailed nonparametric test.
#' 
#' This script was not used for final analysis

# Create file to direct console output
sink(sprintf("%s/Wilcox.test_GDSC.IC50.PCC_cell.lines.x.hm.sig%sfold_Top%s.txt", output_folder, FC, TOP_WHAT))

# NORMALITY TEST: 
# It's possible to use a significance test comparing the sample distribution to a normal one in order to ascertain whether data 
# show or not a serious deviation from normality.
# There are several methods for normality test such as Kolmogorov-Smirnov (K-S) normality test and Shapiro-Wilk's test.
# The null hypothesis of these tests is that "sample distribution is normal". If the test is significant, the distribution is non-normal.
# Shapiro-Wilk's method is widely recommended for normality test and it provides better power than K-S. 
# It is based on the correlation between the data and the corresponding normal scores.
# From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
# In other words, we can assume the normality.

# Compute normality
x <- shapiro.test(unlist(pcc_matrix))
print(x)
# Print result
x <- sprintf("If p.value %s is greater 0.05, we can assume normality", x$p.value)
print(x)

# SIGNIFICANCE TESTINGG
# Perform Wilcoxin signed rank test - a non-parametric statistical hypothesis test used to compare two related samples, matched samples, 
# or repeated measurements on a single sample to assess whether their population mean ranks differ 
# (i.e. it is a paired difference test). It can be use as an alternative to the paired Student's t-test 
# (also known as "t-test for matched pairs" or "t-test for dependent samples") 
# when the sample size is small and the population cannot be assumed to be normally distributed

# Compute significance
x <- wilcox.test(pcc_matrix$BC, pcc_matrix$LP, paired = TRUE, alternative = "two.sided")
print(x)
# Print result
x <- t.test(pcc_matrix$BC, pcc_matrix$LP, paired = FALSE, alternative = "two.sided")
print(x)

# End file
sink()

