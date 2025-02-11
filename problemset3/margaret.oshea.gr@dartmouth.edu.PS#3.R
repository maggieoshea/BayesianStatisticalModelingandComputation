############################################
##  file: margaret.oshea.gr@dartmouth.edu_PS3.R
## Written on R Version 2024.12.0+467
######################################################
##  Maggie O'Shea
##  copyright by the author
##  distributed under the GNU general public license
##  https://www.gnu.org/licenses/gpl.html
##  no warranty (see license details at the link above)
######################################################
##   Course: Bayesian Statistical Modeling & Computation
##   Professor Klaus Keller
##   February 14, 2025
##   Problem Set #3
##################################################
# contact: margaret.oshea.gr@dartmouth.edu
##################################################
# sources:
# Density plots in R: https://www.geeksforgeeks.org/histograms-and-density-plots-in-r/
# Integration in R: https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
# R 'norm' funtions: https://seankross.com/notes/dpqr/
# Qian, S. S., Stow, C. A., & Borsuk, M. E. (2003). On Monte Carlo methods for Bayesian inference. Ecological Modelling, 159(2-3), 269–277.
# D’Agostini, G. (2003). Bayesian reasoning in data analysis: A critical introduction. Singapore: World Scientific Publishing. (Chapter 6 only).
# Ruckert, K. L., Guan, Y., Bakker, A. M. R., Forest, C. E., & Keller, K. (2017). The effects of time-varying observation errors on semi-empirical sea-level projections. Climatic Change, 140(3-4), 349–360. https://doi.org/10.1007/s10584-016-1858-z
# Kim, Y., Bang, H., Kim, Y. & Bang, H. Introduction to Kalman Filter and Its Applications. in Introduction and Implementations of the Kalman Filter (IntechOpen, 2018). doi:10.5772/intechopen.80600.
# Conversation with Anna Valentine

####################################################
# Clear any existing variables and plots. 
rm(list = ls())
graphics.off()

## Packages ##
# If packages ggplot2 and dplyr not already downloaded, un-tag (remove the #) to download packages before running the rest of the script.
#install.packages("ggplot2")
#install.packages("dplyr")
library(ggplot2)
library(dplyr)

## Q2 ##
# Draw the probability density function for the usable fuel in the tank without 
#any other information besides the fuel gauge reading. Determine the expected 
#value of available fuel, the most likely value of available fuel, and probability 
#of negative fuel in the tank. 

# define a seed for reproducibility
set.seed(930)

# Number of trials 
n_trials <- 10^5

# sample from normal distribution with observed fuel as mean, and sd as error in fuel sensor readings
samples <- rnorm(n_trials, 34, 20)


expected_value <- mean(samples)
density_estimate <- density(samples)
df <- approxfun(density(samples))

# Plot PDF
plot(density(samples), main="", 
     xlab="Fuel in Tank (liters)", ylab="PDF")
abline(v=expected_value,col="blue",lty=1,lwd=4)
abline(v=0,col="black",lty=2,lwd=2)
legend("topright", c("Expected Value","Zero Fuel"),
       lwd=c(1,2,2), lty=c(1,2,2), col=c("blue","black"), cex=0.75)
title(main="PDF of Usable Fuel based on\nFuel Gauge Reading of 34 Liters")

# Find the probability of negative fuel

# density function (interpolating across points in density() values to get function to integrate)
density_function <- approxfun(density_estimate$x, density_estimate$y, rule=2)
# integrate from negative infinity to 0 
prob_less_than_zero <- integrate(density_function, lower=-Inf, upper=0)$value
print(prob_less_than_zero)

# Find the most probable value 
max_probability = max(density_estimate$y)
index <- which(density_estimate$y == max_probability) 
mostprobablevalue <- density_estimate$x[index] 
print(mostprobablevalue)

## Q4 ##
#Use a grid-based method to determine your Bayesian update from your prior 
#and the likelihood function. Add this posterior to the plot produced above. 
# Determine now the probability of negative fuel. Has this fixed the issue? If so, how?

# New Prior is based on uniform distribution representing the capacity of the tank
# max is based on total fuel capacity of tank
max=182
min=0
observed = 34
tanksd = 20
# grid-based method
gridvalues = seq(min, max, length.out=182)

# prior is uniform distribution from 0 to 182
uniform_prior <- dunif(gridvalues, min=min, max=max)

# likelihood 
likelihood_values <- dnorm(gridvalues, mean=observed, sd=tanksd)

posterior <- likelihood_values*uniform_prior
posterior_normed <- posterior / sum(posterior)

# Plot CDF 
dx = (gridvalues[2] - gridvalues[1])
cdf_posterior <- cumsum(posterior_normed) * dx
plot(gridvalues, cdf_posterior, type = "l", lwd = 5, col = "blue", 
     main = "CDF of Posterior Distribution", xlab = "Fuel in Tank (liters)", ylab = "CDF")
abline(v=0, col='black', lty=2, lwd=2)
legend("bottomright", legend = c("CDF", '0 Fuel'),
       col = c("blue", 'black'), lwd = 2,
       lty = c(1,2),
       cex = 1)



# Find probability of negative fuel 
zero_grid_values_index <- gridvalues<0
probnegfuel <- sum(posterior_normed[zero_grid_values_index])
print(probnegfuel)

# Plotting all together
plot(density(samples), col="blue",
     lwd=2,
     ylab= "PDF", 
     xlab='Fuel (liters)',
     main="Comparing Updated Posterior with Original PDF")
lines(gridvalues, posterior_normed, col = "red", lwd = 2)

# Add uniform distribution 
lines(gridvalues, uniform_prior, col = "brown", lwd = 3) 
segments(0, 0, 0, max(uniform_prior), col = "brown", lwd = 3)
segments(182, 0, 182, max(uniform_prior), col = "brown", lwd = 3)

# Add vertical lines for expected values
abline(v = 0, col = "black", lty = 2, lwd = 2)
abline(v = observed, col = "green", lty = 2, lwd = 2)

legend("topright", legend = c("Q2 PDF", "Q4 Posterior", "Uniform Prior", "Observed 34 Liters", "0 Liters"),
       col = c("blue", "red", "brown", "green", "black"), 
       lwd = 2, 
       lty = c(1, 1, 1, 2, 2),  # Add line types here, 2 for dashed
       cex = 0.75)

## Q5 ##
# Repeat the step above using a Bayes Monte Carlo method. 

# define a seed for reproducibility
set.seed(930)

# randomly drawing n samples from parameter distribution 
# min/max from tank capacity
max=182
min=0
n_trials = 30*10^5
prior_samples <- runif(n_trials, min=min, max=max)

# Likelihood -- 
# likelihood is dnorm, so combining the samples with dnorm is posterior
posterior_values <- dnorm(prior_samples, mean = observed, sd = tanksd)

# Normalize the posterior
BMC_posterior <- posterior_values / sum(posterior_values)

# Posterior Sample
BMC_posterior_samples <- sample(prior_samples, size = n_trials, replace = TRUE, prob = BMC_posterior)

plot(density(BMC_posterior_samples), col = "purple", lwd = 2, 
     xlab = "Fuel in Tank (liters)", ylab = "Density", main = "Posterior from Bayes Monte Carlo (BMC)")


# Test convergence of BMC #
convergence_test <- 2 * qt(0.975, length(BMC_posterior_samples) - 1) * sd(BMC_posterior_samples) / sqrt(length(BMC_posterior_samples))

# Now check if the CI is sufficiently small (<= 0.05)
if (convergence_test <= 0.05) {
  print("Confidence Interval is sufficiently small")
} else {
  print("Confidence Interval is too large")
}

# Plotting all together
plot(density(samples), col="blue",
     lwd=2,
     ylab= "PDF", 
     xlab='Fuel (liters)',
     main="Comparing Posteriors: PDF, Updated Posterior, & BMC",
     ylim = c(0, 0.021))
lines(gridvalues, posterior_normed, col = "red", lwd = 2)
lines(density(BMC_posterior_samples), col = "purple", lwd = 4)

# Add vertical lines for expected values
abline(v = 0, col = "black", lty = 2, lwd = 2)
abline(v = observed, col = "green", lty = 2, lwd = 2)
legend("topright", legend = c("Q2 PDF","Q4 Posterior", "Q5 BMC Posterior", "Observed 34 Liters", "0 Liters"),
       col = c("blue", "red", "purple", "green",  "black"), 
       lty = c(1, 1, 1, 2, 2),
       lwd = 2, cex = 0.75)


# Find probability of negative fuel 
zero_grid_values_index_BMC <- BMC_posterior_samples<0
probnegfuel_BMC <- sum(BMC_posterior[zero_grid_values_index_BMC])


## Q7 ##
# Produce a plot of the estimated available flight time. 
#What is the probability that you make an airport that is 100 minutes flight time
# away with at least 30 min reserve fuel required by regulations?
#What is the probability that you run out of fuel trying to make it to this airport?
fuelperminute = 18/60
sd_fuelperminute = 2/60

# 100 minutes plus 30 minutes reserve
fuel_needed_100min = fuelperminute*130
fuel_needed_100min_upper = (fuelperminute*130)+(sd_fuelperminute*130)
fuel_needed_100min_lower = (fuelperminute*130)-(sd_fuelperminute*130)


# Probability of having enough fuel for 130 minutes
min130_grid_values_index <- gridvalues<=fuel_needed_100min
probenoughfuel <- sum(posterior_normed[min130_grid_values_index])
print(probenoughfuel)

min130_grid_values_index_upper <- gridvalues<=fuel_needed_100min_upper
probenoughfuel_upper <- sum(posterior_normed[min130_grid_values_index_upper])
print(probenoughfuel_upper)

min130_grid_values_index_lower <- gridvalues<=fuel_needed_100min_lower
probenoughfuel_lower <- sum(posterior_normed[min130_grid_values_index_lower])
print(probenoughfuel_lower)



# Probability of running out of fuel 
minimum_fuel_needed_100min = fuelperminute*100
minimum_needed_100min_upper = (fuelperminute*100)+(sd_fuelperminute*100)
minimum_needed_100min_lower = (fuelperminute*100)-(sd_fuelperminute*100)

fuelmin_grid_values_index <- gridvalues<=minimum_fuel_needed_100min
probrunningout <- sum(posterior_normed[fuelmin_grid_values_index])
print(probrunningout)

fuelmin_grid_values_index_upper <- gridvalues<=minimum_needed_100min_upper
probrunningout_upper <- sum(posterior_normed[fuelmin_grid_values_index_upper])
print(probrunningout_upper)

fuelmin_grid_values_index_lower <- gridvalues<=minimum_needed_100min_lower
probrunningout_lower <- sum(posterior_normed[fuelmin_grid_values_index_lower])
print(probrunningout_lower)


# Plot results on PDFs
plot(density(samples), col="blue",
     lwd=2,
     ylab= "PDF", 
     xlab='Fuel (liters)',
     main="Comparing Posteriors: PDF, Updated Posterior, & BMC\nWith Needed Fuel Estimates",
     ylim=c(0, 0.021))
lines(gridvalues, posterior_normed, col = "red", lwd = 2)
lines(density(BMC_posterior_samples), col = "purple", lwd = 2)

# Add the vertical lines
abline(v = 0, col = "black", lty = 2, lwd = 2)
abline(v = observed, col = "green", lty = 2, lwd = 3)
abline(v = fuel_needed_100min , col = "orange", lty = 2, lwd = 2)
abline(v = minimum_fuel_needed_100min, col = "cyan1", lty = 2, lwd = 2)

# Fill the area between the upper and lower values (uncertainty) around estimates
polygon(c(fuel_needed_100min_lower, fuel_needed_100min_lower, fuel_needed_100min_upper, fuel_needed_100min_upper), 
        c(0, 0.055, 0.055, 0), col = rgb(1, 0.647, 0, alpha = 0.3), border = NA)
polygon(c(minimum_needed_100min_lower, minimum_needed_100min_lower, minimum_needed_100min_upper, minimum_needed_100min_upper), 
        c(0, 0.055, 0.055, 0), col = rgb(0, 1, 1, alpha = 0.3), border = NA)

# Legend
legend("topright", legend = c("Q2 PDF", "Q4 Posterior", "Q5 BMC Posterior", "Observed 34 Liters", "0 Liters",
                              "Fuel Needed for 100 min + Reserve", "Minimum Fuel for 100 minutes"),
       col = c("blue", "red", "purple", "green", "black", "orange", "cyan1"), 
       lwd = 2, 
       lty = c(1, 1, 1, 2, 2, 2, 2),
       cex = 0.75)


minuteperfuel = 1/fuelperminute
sd_minuteperfuel = 1/sd_fuelperminute

available_flighttime = gridvalues*minuteperfuel
available_flighttime_upper = (gridvalues*minuteperfuel)+(gridvalues*sd_fuelperminute)
available_flighttime_lower = (gridvalues*minuteperfuel)-(gridvalues*sd_fuelperminute)

# Plot Estimated available flight time relative to fuel 
plot(x = gridvalues, y=(available_flighttime/60),
     type='l', 
     lwd=1,
     ylab= "Flight Time (hours)", 
     xlab='Fuel (liters)',
     main="Flight Time vs Fuel Needed")
polygon(c(gridvalues, rev(gridvalues)),
        c((available_flighttime_upper / 60), rev(available_flighttime_lower / 60)),
        col = rgb(1, 0, 0, 0.3, alpha=0.7), border = NA)
lines(x = gridvalues, y=(available_flighttime/60), col = "black", lwd = 2)
legend("bottomright", legend = c("Estimate", "Upper and Lower Bound"),
       col = c("black", "red"), lwd = 2, cex = 0.75)






















