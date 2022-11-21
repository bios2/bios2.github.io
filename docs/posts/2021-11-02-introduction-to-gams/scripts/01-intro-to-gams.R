# try out some things in R

library(mgcv)
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)

#loading the data ####

trawls <- read.csv("data/trawl_nl.csv")

head(trawls)

# This data is from the Newfoundland and Labrador fall trawl survey (2J, 3K,
# and 3L), from 2005-2014. This spans the period where cod biomass started
# increasing again and Northern Shrimp started declining. The data includes the
# following variables:
#
# year: survey year, from 2005-2014
# trip_id: The unique id for the trawl. Includes year, vessel, trip, and set
# stratum: survey stratum
# area_trawled: The estimated area swept by the gear (km^2)
# long: longitude, in decimal-degrees
# lat: latitude.
# temp_bottom: average bottom temperature of trawl, in degrees C
# depth: Average trawl depth in meters.
# total: Total biomass caught by the trawl, scaled by area swept (kg/km^2)
# shrimp: Northern Shrimp (Pandalus borealis) biomass per unit area caught in the trawl (kg/km^2)
# cod: Atlantic Cod (Gadus morhua) biomass per unit area caught in the trawl (kg/km^2)
# richness: Total number of distinct taxonomic units caught in the trawl
# x: longitude in UTM coordinates (meters)
# y: latitude in UTM coordinates (meters)


# 1. 1d models####

#  - simple models with `s()`

# inspect data
ggplot(trawls, aes(x = depth, y= temp_bottom))+
  geom_point()+
  scale_x_log10()

ggplot(trawls, aes(x = depth, y= temp_bottom))+
  geom_point()+
  geom_smooth(method = lm, formula = y~x)+
scale_x_log10()

# try a simple model ####
model_temp_k10 <- gam(temp_bottom ~  s(log10(depth), k = 10, bs="tp"),
                      data=trawls,
                      method="REML")

# k = 10 sets the number of basis functions for the smoother to 10. This sets
# the maximum variability of this smoother. For technical reasons which we'll
# talk about in the afternoon, when k=10, the final smoother only has 9 basis
# functions when fit.

# bs = "tp" says we want to use a thin-plate smoother, which is a different type
# of spline smoother than the cubic splines we talked about in class. k = 10 and
# bs = "tp" are the default settings, so you could just write s(log10(depth))
# here.

# Note that we're using method = "REML". We really recommend specifically
# setting method to REML here!


#  - `summary` and `plot`
# what does that look like?
plot(model_temp_k10)

# We can also use the gratia package
draw(model_temp_k10)

# is it flexible enough?
summary(model_temp_k10)

# Things to look for when evaluating the model: edf: the effective degrees of
# freedom: basically, how much effect is the penalty actually having on the shape
# of the curve? If edf is very close to the number of basis functions (here 9),
# it may be a sign that you are not using enough basis functions to model the
# actual wiggliness of the function.


# increase k ####

model_temp_k30 <- gam(temp_bottom~s(log10(depth), k=30),
                      data=trawls,
                      method="REML")
draw(model_temp_k30)


# is it flexible enough?
summary(model_temp_k30)
# better



# Let's also see what role the smoothing parameter is playing here: ####

model_temp_k30$sp
#should be  about 0.03

#what happens if we increase that by 100?

model_temp_k30a <- gam(temp_bottom~s(log10(depth), k=30,sp = 3),
                       data=trawls,
                       method="REML")
draw(model_temp_k30a)

#What if we reduce the penalty to zero?
model_temp_k30b <- gam(temp_bottom~s(log10(depth), k=30,sp = 0),
                       data=trawls,
                       method="REML")
draw(model_temp_k30b)

#What if we really increase the penalty?
model_temp_k30c <- gam(temp_bottom~s(log10(depth), k=30,sp = 3000),
                       data=trawls,
                       method="REML")
draw(model_temp_k30c)




# ---- Exercise 1 ----

# Try modelling the natural log of total trawl biomass as a function
# log10(depth). plot the function, look at the model summary, and determine,
# using edf, if you think that the default k = 10 is sufficient degrees of
# freedom to model this relationship

# Here's a starting point for your model:


model_biomass <- gam(log(total) ~ [...],
                      data=trawls,
                      method="REML")





#  Adding more than one smooth to your model ####

# add a year effect:
model_tempyear <- gam(temp_bottom~s(log10(depth), k=20) + s(year, k=10),
                      data=trawls,
                      method="REML")


summary(model_tempyear)
# what does that look like?
plot(model_tempyear, pages=1)

draw(model_tempyear)

#Also useful to look at the penalty terms here:
model_tempyear$sp


# We can also use a different smoother here; for instance, we could treat year
# as a random effect. This means year-effects are not smoothed toward nearby
# years, but all the average temperature across all years will be pulled toward
# the overall average value:

# First we need to create a factor variable for years:

trawls$yearf = factor(trawls$year)

model_tempyearf <- gam(temp_bottom~s(log10(depth), k=20) + s(yearf, bs = "re"),
                      data=trawls,
                      method="REML")

plot(model_tempyearf, pages = 1)

draw(model_tempyearf)

summary(model_tempyearf)


# ---- Exercise 2 ----

# Try extending your model of the log of total trawl biomass by adding a smooth
# for temperature and year. plot the functions, look at the model summary, and
# determine, using edf, if you think that the default k = 10 is sufficient
# degrees of freedom to model this relationship



model_biomass_3d <- gam(log(total) ~ s(log10(depth), bs="tp", k = 30)+
                           [...],
                         data=trawls,
                         method="REML")

draw(model_biomass_3d)
summary(model_biomass_3d)




