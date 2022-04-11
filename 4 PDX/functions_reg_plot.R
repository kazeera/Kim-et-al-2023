# Make linear model/ linear regression and extract linear formula (y=mx+b), R^2 and p vale as a string to plot on graph
# Note: the output is in expression format which allows for italicizing, exponent, and other formatting
# Source: https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df){
  m <- lm(formula = angle ~ value, data = df);
  x <- summary(m)
  
  # Make an expression for the plot
  # Refer https://www.r-bloggers.com/math-notation-for-r-plot-titles-expression-and-bquote/
  eq <- substitute(italic(y) == b%.%italic(x) + a~","~~italic(R)~"="~pearson_corr*","~~italic(r)^2~"="~r2*~","~~italic(p)~"="~pval, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pval = format(summary(m)$coefficients[2,4], digits = 3),
                        pearson_corr = format(cor(df$angle, df$value), digits = 3)) #https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R5_Correlation-Regression/R5_Correlation-Regression_print.html
  )
  as.character(as.expression(eq));
}

# How do I change the number of decimal places on axis labels in ggplot2?
# https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.1f", x)


# If R^2 goes up, better fit
# lm( angle ~ ssGSEA)
# lm(angle ~ ssGSEA + cancer_type)
# lm(angle ~ ssGSEA * cancer_type)
