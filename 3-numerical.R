#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE 4. NUMERICAL BEHAVIOUR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script performs numerical simulations to see whether a zero sum assumption
#   is met.
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
dat %>% compadre_calc_simulation() %>% simulation_plot() +ylim(-60,0) #FIGURE4A


#Animal Matrix
load("D:/0_Files/2_ZeroSum/COMADRE_v.4.23.3.1.RData")
ani <- comadre %>% compadre_filter()
ani %>% compadre_calc_simulation %>% simulation_plot() + ylim(-60,0) #FIGURE4B



