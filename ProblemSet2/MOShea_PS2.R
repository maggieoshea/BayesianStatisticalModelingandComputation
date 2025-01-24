############################################
##  file: MOShea_PS2.R
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
##   January 24, 2025
##   Problem Set #2
##################################################
# contact: margaret.oshea.gr@dartmouth.edu
##################################################
# sources:
# - Applegate, P. J., & Keller, K. (Eds.). (2016). Risk analysis in the Earth Sciences: A Lab manual. 2nd edition. Leanpub. Retrieved from https://leanpub.com/raes
# - Callahan, C. W. & Mankin, J. S. Globally unequal effect of extreme heat on economic growth. Science Advances 8, eadd3726 (2022).
# - Gottlieb, A. R. & Mankin, J. S. Evidence of human influence on Northern Hemisphere snow loss. Nature 625, 293–300 (2024).
# - coin.R file from 'Bayesian Stats' course with Dr. Klaus Keller
# - R help files accessed through R-studio for syntax
# - Class discussions in Bayesian Statistical Modeling & Computation, Dartmouth College, Winter 2025
# - Online R coding help, namely: 
#   - For loops: https://www.w3schools.com/r/r_for_loop.asp
#   - Binary Variable: https://stackoverflow.com/questions/19970009/create-binary-column-0-1-based-on-condition-in-another-column
#   - Slice_min: https://stackoverflow.com/questions/24070714/extract-row-corresponding-to-minimum-value-of-a-variable-by-group
# - Conversation with Alexis Hudes (specifically for B) and Anna Valentine (about convergence)
####################################################

## Packages ##
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)

### Instructions: ###
## Question 4A. Use Monte Carlo method to determine the mean and 95 Percentile from known univariate normal distribution: Mean of zero and SD of one with one seed and one sample size:

# define a seed for reproducibility
set.seed(930)

# Number of trials 
n_trials <- 10^5

samples <- rnorm(n_trials, 0, 1)

# Find the mean #
mean = mean(samples)
print(mean)

# 95% confidence interval: #
standard_deviation = sd(samples)
sqrt_n = sqrt(n_trials)
CI_95_upper = mean+(1.96*(standard_deviation/sqrt_n))
CI_95_lower = mean-(1.96*(standard_deviation/sqrt_n))
print(paste("95% Confidence Interval Upper Bound: ", CI_95_upper))
print(paste("95% Confidence Interval Lower Bound: ", CI_95_lower))

# Plot Results from one Sample #
hist(samples, main = "", xlab = "Samples from Normal Distribution") 
abline(v=mean,col="blue",lty=1,lwd=4)
abline(v=0,col="yellow",lty=2,lwd=4)
abline(v=CI_95_lower, col="red", lty=2, lwd=2)
abline(v=CI_95_upper, col="red", lty=2, lwd=2)
legend("topleft", c("Mean","Prior Expectation",
                    "95% confidence interval"),
       lwd=c(1,2,2), lty=c(1,2,2), col=c("blue","yellow","red"), cex=0.75)
title(main="Sample Means from Random Sampling\n of Normal Distribution N(0,1)")

## Characterize Uncertainty from Sample Size and Seed ##

# Vary seed and sample size
seeds <- seq.int(930, 1030, by = 1) # 100 seeds
sample_sizes <- seq.int(1000, 50000, by = 1000) # 500 sample sizes

# Prepare dataframe to save results from each seed and sample size combination
mean_df <- data.frame(
  seed = integer(),
  sample_size = integer(),
  mu = numeric(),
  CI_upper = numeric(),
  CI_lower = numeric()
)

# for each seed
for (seed in seeds) {
  set.seed(seed)
  
  # for each sample size
  for (size in sample_sizes){
    samples <- rnorm(size, 0, 1)
    
    # calculate the mean
    mean = mean(samples)
    
    # calculate CI
    standard_deviation = sd(samples)
    sqrt_n = sqrt(size)
    CI_95_upper = mean+(1.96*(standard_deviation/sqrt_n))
    CI_95_lower = mean-(1.96*(standard_deviation/sqrt_n))
    
    # append this data to the dataframe as a new row
    mean_df <- rbind(mean_df, data.frame(seed = seed, sample_size = size, mu = mean, CI_upper = CI_95_upper,  CI_lower = CI_95_lower))
    
    
  }
}

# Plot with ggplot2 the results from 500 sample sizes and 100 seeds
ggplot(data = mean_df, aes(x = sample_size, y =mu, color = factor(seed), group = seed)) + 
  geom_line() +  
  labs(title = "Characterizing Seed & Sample Size Uncertainty\n in Univariate Normal Distribution Sample Means",
       subtitle = "Colored by Unique Seed",
       x = "Number of Samples",
       y = "Estimate",
       color = "Seed") + 
  theme_minimal() +  
  theme(legend.position = "none")  

## Creating a "mean" dataframe which takes the average CI and mean across all seeds
CI_mean <- mean_df %>%
  group_by(sample_size) %>%
  summarise(avg_upperCI = mean(CI_upper),
            avg_lowerCI = mean(CI_lower),
            avg_mu = mean(mu))

CI_mean$differencefromtruth = (0-CI_mean$avg_mu)

colors <- c("Mean" = "red", "Confidence Interval" = "black")

ggplot(data = CI_mean, aes(x = sample_size)) + 
  geom_line(data = CI_mean, aes(y = avg_upperCI))+
  geom_line(data = CI_mean, aes(y = avg_lowerCI, color = 'Confidence Interval'))+
  geom_line(data = CI_mean, aes(y = avg_mu, color = 'Mean'))+

  labs(title = "Characterizing Sample Size Uncertainty",
       subtitle = "Averaged Across 100 Seeds",
       x = "Number of Samples",
       y = "Estimate",
       color = "Legend") + 
  theme_minimal() +  
  theme(legend.position = "right")   +
  scale_color_manual(values = colors)

### Instructions: ###
## Quesion 4B. Determine the value of pi with your estimated uncertainties ##
# Sources of uncertainty to characterize: sample size and seeds

# Vary number of samples and seeds
seeds <- seq.int(930, 1130, by = 1) # 200 seeds
sample_sizes <- seq.int(100, 10000, by = 100) # number of trials in each seed sample

# Prep dataframe to hold results with different seeds and sample sizes 
results_df <- data.frame(
  seed = integer(),
  sample_size = integer(),
  pi_estimate = numeric()
)

# for each seed
for (seed in seeds) {
  set.seed(seed)
  # for each sample size
  for (sample in sample_sizes){
    # sample size is the number of trials in each seed run
    # The number of trials is the number of samples in the circle/square space
    sample_size=sample
    n_trials <- seq.int(0, sample_size, by = 1) 
    
    # prep circle and square variables to be able to add up the number of circle points vs square points
    circle <- 0
    square <- 0
    
    
    for (trial in n_trials) {
      # include print statement if interested in tracking which trial run
      #print(paste("running trial #", trial))
      n = 1 # 1 xy combo per trial
      
      # random number in uniform distribution from 1 to -1
      x_point <- runif(n, -1, 1)
      y_point <- runif(n, -1,1)
      
      # calculate the distance between random xy point and origin
      distance_from_origin = sqrt(x_point^2 + y_point^2)
      
      # if distance is less than radius 1, then it's in the circle, otherwise in square 
      if (distance_from_origin <= 1){
        circle = circle +1}
      else
        square = square +1 }
    
    
    # proportion of circle in this space. 
    # The denominator is square + circle because anywhere in circle is also in square
    circle_prop = circle/(square+circle)
    pi_value = circle_prop*4
    
    # append results from this run to the dataframe for results
    results_df <- rbind(results_df, data.frame(seed = seed, sample_size = sample_size, pi_estimate = pi_value))
    
    
  }
  }


# Estimating Convergence...
# Find the smallest sample size at which point the different between the truth is within 0.001
earliest_meeting_threshold <- results_df %>%
  group_by(seed) %>%  # Group by the seed column
  filter(pi_estimate >= pi - 0.001 & pi_estimate <= pi + 0.001) %>% 
  slice_min(order_by = sample_size) # filter for only columns that meet threshold


# plot results from each seed and sample size using ggplot2
ggplot(data = results_df, aes(x = sample_size, y = pi_estimate, color = factor(seed), group = seed)) + 
  geom_line() +  
  labs(title = "Characterizing Seed & Sample Size Uncertainty in Estimates of Pi",
       subtitle = "Colored by Unique Seed",
       x = "Number of Samples",
       y = "Estimated Pi Value",
       color = "Seed") +  # Label for the color legend
  theme_minimal() +  
  theme(legend.position = "none")  

# calculate mean pi estimates and CI across all seeds, grouped by sample sizes
mean_pi_df <- results_df %>%
  group_by(sample_size) %>%
  summarise(sd_pi = sd(pi_estimate), 
            mean_pi = mean(pi_estimate),
            # sqrt(100) is the sqrt of n, the number of seeds per sample size
            avg_upperCI = mean_pi+(1.96*sd_pi/(sqrt(100))),
            avg_lowerCI = mean_pi-(1.96*sd_pi/(sqrt(100))))

# Find where the pi estimate is within a threshold of true pi (0.001)
mean_pi_df$difference_from_truepi = abs(pi-mean_pi_df$mean_pi)
mean_pi_df$meetsthreshold <- ifelse(mean_pi_df$difference_from_truepi<=0.001, 1, 0)


# Plot results from average mean and CI by sample size
colors <- c("Mean" = "red", "Confidence Interval" = "black")

ggplot(data = mean_pi_df, aes(x = sample_size)) + 
  geom_line(data = mean_pi_df, aes(y = avg_upperCI))+
  geom_line(data = mean_pi_df, aes(y = avg_lowerCI, color = 'Confidence Interval'))+
  geom_line(data = mean_pi_df, aes(y = mean_pi, color = 'Mean'))+
  # Add line for true pi
  #geom_hline(yintercept= pi, color = "blue", size=0.5) + 
  
  labs(title = "Characterizing Sample Size Uncertainty",
       subtitle = "Averaged Across 200 Seeds",
       x = "Number of Samples",
       y = "Pi Estimate",
       color = "Legend") + 
  theme_minimal() +  
  theme(legend.position = "right")   +
  scale_color_manual(values = colors)


