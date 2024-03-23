#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE 2-3. VARIATIONS IN ASYMPOTIC BEHAVIOUR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script calculates the dominant eigenvalues in plant and animal databases
#   with respect to different life forms 
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
print(eigen_plot_Type(eigens$domin, dat$metadata$OrganismType))     #FIGURE 2A
geometric_mean_test_Type(eigens$domin, dat$metadata$OrganismType)   #FIGURE 2A-test
print(SAC_plotType(eigens$domin, dat$metadata$OrganismType))        #FIGURE 3A   

#Animal Matrix
load("D:/0_Files/2_ZeroSum/COMADRE_v.4.23.3.1.RData")
ani <- comadre %>% compadre_filter()
ani_eigens <- ani %>% compadre_calc_domin_eigen()
print(eigen_plot_Type(ani_eigens$domin, ani$metadata$OrganismType))   #FIGURE 2B
geometric_mean_test_Type(ani_eigens$domin, ani$metadata$OrganismType) #FIGURE 2B-test
print(SAC_plotType(ani_eigens$domin, ani$metadata$OrganismType))          #FIGURE 3B   



