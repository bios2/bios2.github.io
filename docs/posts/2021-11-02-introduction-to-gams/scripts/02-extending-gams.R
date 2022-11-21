# example and practical scripts for part 02!

# load requisite packages ####
library('here')
library('mgcv')
library('tidyr')
library('dplyr')
library('ggplot2')
library('gratia')
library('patchwork')

# Loading data ####

#load the trawls data set from last session
trawls <- read.csv("data/trawl_nl.csv")

# Let's look at a new data set as well: counts of positive influenza tests from
# California in 4 years 

flu_data <- read.csv("data/california-flu-data.csv") %>%
  mutate(season = factor(season))

#This data consists of the following columns:

# season: the flu season. Starts in early October of the first year, goes to
#         September of the next
# year: the year a given observation was taken in
# week: the week of the year (from 1-52)
# week_centered: week centered around January of the given season. Negative
#                values correspond to weeks before January, pos to weeks after
# tests_pos: number of tests positive for at least 1 strain of influenza that week
# tests_total: total number of tests ran in a given week



# A. Modelling non-normal data ####

ggplot(trawls, aes(x =depth, y=richness)) +
  geom_point()+
  scale_y_log10()+
  labs(x="Depth (m)", y="Number of Species")

# Does richness vary with depth?


rich_depthl10 <- gam(richness~s(log10(depth), k=30), 
                     family=poisson,
                     method="REML", 
                     data=trawls)

plot(rich_depthl10, shade=TRUE)

summary(rich_depthl10)

# Are we using the right distribution?

par(mfrow=c(2,2))
gam.check(rich_depthl10)

appraise(rich_depthl10)



# ---- Exercise 3 ----

# Try modelling the species richness of the trawl data using a negative binomial family 
# as a function of depth. Use the appraise function to evaluate the model fit

# Here's a starting point for your model:


rich_depthl10_negbin <- gam(richness ~ s(log10(depth), k=30),
                            data=trawls,
                            method="REML",
                            family = nb(theta = NULL,link = "log"))

# Non-normal flu-data example ####

# plot the data

flu_plot_raw <- ggplot(flu_data, 
                       aes(x =  week_centered, 
                           y= tests_pos,color=season))+
  facet_wrap(~season)+
  geom_line()+
  geom_point()+
  scale_color_brewer(palette = "Set1")

#filter to just one year for now
flu_data_2010 <- flu_data %>%
  filter(season=="2010-2011")

# Let's try smoothing this using a Poisson model:

flu_2010_tests_mod <- gam(tests_pos ~ s(week_centered, k = 15),
                          data= flu_data_2010, 
                          family = poisson,
                          method = "REML")

summary(flu_2010_tests_mod)
draw(flu_2010_tests_mod)

# However, we should account for the fact that the number of tests being done
# We can do this with an offset:

flu_2010_posrate_mod <- gam(tests_pos ~ s(week_centered, k = 15) + offset(log(tests_total)),
                            data= flu_data_2010, 
                            family = poisson,
                            method = "REML")

summary(flu_2010_posrate_mod)
draw(flu_2010_posrate_mod)

# Finally, we might want to use a cyclic smoother here, to account for the fact
# that the start and end of the flu season may have similar positivity rates:


flu_2010_posrate_cc_mod <- gam(tests_pos ~ s(week_centered, k = 15,bs="cc") + offset(log(tests_total)),
                               data= flu_data_2010, 
                               family = poisson,
                               method = "REML")

summary(flu_2010_posrate_cc_mod)
draw(flu_2010_posrate_cc_mod)





# B. Two-dimensional models  ####

# subset to just 2010's trawls
trawls_2010 <- filter(trawls, year==2010)
range(trawls_2010$shrimp)
head(trawls_2010$shrimp, 3)

hist(trawls_2010$shrimp, xlab="Biomass", main="")

spatial_shrimp <- gam(shrimp ~ s(x, y),
                      data=trawls_2010,
                      method="REML")

plot(spatial_shrimp, asp=1)


spatial_shrimp_tw <- gam(shrimp ~ s(x, y, k=100),
                      data=trawls_2010,
                      family="tw",
                      method="REML")


draw(spatial_shrimp_tw, asp=1)
appraise(spatial_shrimp_tw)

# Models with interactions ####

shrimp_te <- gam(shrimp ~ te(log10(depth), temp_bottom),
                 data=trawls_2010,
                 family=tw,
                 method="REML")


summary(shrimp_te)

draw(shrimp_te)

shrimp_ti <- gam(shrimp ~ ti(log10(depth), temp_bottom)  +
                   s(temp_bottom) + s(log10(depth)),
                 data=trawls_2010,
                 family=tw,
                 method="REML")


draw(shrimp_ti)


## Exercise 4 ####

# we tried fitting a Gaussian to the shrimp in the lecture
# what about other distributions?

# 1. copy this model and try out using the tw() and the quasi()  distribution in place
#    of __OTHER_DISTRIBUTIONS__. Remember to rename your new models!
# 2. look at the summary() output of these models and compare to those in
#    the slides
# 3. if you have time, plot the models too

# spatial shrimp biomass
b_shrimp_template <- gam(shrimp ~ s(x, y, k =100),
                         data=trawls_2010,
                         family=__OTHER_DISTRIBUTIONS__,
                         method="REML")




# ---

## Exercise 5: creating a simple species distribution model ####

# we looked at te() and ti() in the slides, let's try a model with
# space and other effects

# 1. modify the model below to:
#   a. have a smooth of x,y and a tensor of temperature and depth, you
#       can choose whether the tensor is using te() or ti()
#   b. use tw(), nb() or quasipoisson() as the response distribution
# 2. look at the summary of this model and compare to those in the slides
# 3. plot the results!

mega_shrimp <- gam(shrimp ~ __YOUR_FORMULA_HERE__,
                   data=trawls_2010,
                   family=__OTHER_DISTRIBUTIONS__,
                   method="REML")


# Spatiotemporal models ####


shrimp_xyt <- gam(shrimp ~ ti(y,x, year, d=c(2,1), 
                              bs=c("tp", "cr"), k=c(20, 5)) +
                    s(x, y, bs="tp", k=20) +
                    s(year,  bs="cr", k=5),
                  data=trawls,
                  family=tw,
                  method="REML")


summary(shrimp_xyt)
plot(shrimp_xyt, select=1, scheme=2, asp=1)

par(mfrow=c(1, 2))
plot(shrimp_xyt, select=2, scheme=2, asp=1)
plot(shrimp_xyt, select=3, scheme=2, scale=0, seWithMean=TRUE)

plot(shrimp_xyt, select=1, scheme=2, asp=1)


