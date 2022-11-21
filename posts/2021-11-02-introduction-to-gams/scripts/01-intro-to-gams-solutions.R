# try out some things in R

library(mgcv)
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)

#loading the data ####

trawls <- read.csv("data/trawl_nl.csv")

head(trawls)

# ---- Exercise 1 ---- ####

# Try modelling the natural log of total trawl biomass as a function
# log10(depth). plot the function, look at the model summary, and determine,
# using edf, if you think that the default k = 10 is sufficient degrees of
# freedom to model this relationship

# Exercise 1 Solution ####
model_biomass_k10 <- gam(log(total) ~ s(log10(depth), k = 10,bs="tp"),
                      data=trawls,
                      method="REML")

draw(model_biomass_k10)
summary(model_biomass_k10)

#Try increasing k here
model_biomass_k30 <- gam(log(total) ~ s(log10(depth), k = 30,bs="tp"),
                         data=trawls,
                         method="REML")

draw(model_biomass_k30)
summary(model_biomass_k30)


# ---- Exercise 2 ---- ####

# Try extending your model of the log of total trawl biomass by adding a smooth
# for temperature and year. plot the functions, look at the model summary, and
# determine, using edf, if you think that the default k = 10 is sufficient
# degrees of freedom to model this relationship


model_biomass_3d <- gam(log(total) ~ s(log10(depth), bs="tp", k = 30)+
                          s(temp_bottom)+
                           s(year),
                         data=trawls,
                         method="REML")

draw(model_biomass_3d)
summary(model_biomass_3d)


