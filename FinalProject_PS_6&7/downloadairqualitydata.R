############################################
##  file: margaret.oshea.gr@dartmouth.edu#4.R
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
##   Final Project
##################################################
# contact: margaret.oshea.gr@dartmouth.edu
##################################################
# sources:
# Lobell D., & Tsiang, M. "EESS 260 Class Notes" Advanced Statistical Methods for Earth System Analysis. (2012). 
# https://rpubs.com/Rich521/airquality
# https://stackoverflow.com/questions/53089219/specify-path-in-write-csv-function
####################################################

# If packages tidyverse not already downloaded, un-tag (remove the #) to download packages before running the rest of the script.
#install.packages("tidyverse")
library(tidyverse)

data(airquality)
head(airquality)

path <- "/Users/f006z55/Desktop/Classes/BayesianStats/"

filename <- "airquality_data.csv" 
write.csv(airquality, file = paste(path, filename, sep = ""), row.names = FALSE) 
