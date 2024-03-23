#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE 1. DOMINANT EIGENVALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script calculates the dominant eigenvalues in plant and animal databases
#
# Author: Haoran Wu (haoran.wu@wolfson.ox.ac.uk)
# Institute: Environmental Change Institute, University of Oxford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(dplyr)
require(ggplot2)
require(sads)
rm(list = ls())
source("0-function.R")

#Plant Matrix
load("D:/0_Files/2_ZeroSum/COMPADRE_v.6.23.5.0.RData")
dat <- compadre %>% compadre_filter() 
eigens <- dat %>% compadre_calc_domin_eigen()
print(eigens %>% eigen_plot())     #FIGURE 1A
print(SAC_plot(eigens$domin))      #FIGURE 1C
geometric_mean_test(eigens$domin, delta = 0.01)

#Animal Matrix
load("D:/0_Files/2_ZeroSum/COMADRE_v.4.23.3.1.RData")
ani <- comadre %>% compadre_filter()
ani_eigens <- ani %>% compadre_calc_domin_eigen()
print(ani_eigens %>% eigen_plot()) #FIGURE 1B
print(SAC_plot(ani_eigens$domin))
geometric_mean_test(ani_eigens$domin, delta = 0.01)