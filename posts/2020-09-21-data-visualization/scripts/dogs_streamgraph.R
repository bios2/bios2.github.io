# Script to make a streamgraph of the top 10 most popular dog breeds in 
# New York City from 1999 to 2015

# load libraries
library(magrittr) # piping
library(lubridate) # dealing with dates
library(dplyr) # data manipulation
library(streamgraph) #devtools::install_github("hrbrmstr/streamgraph")
library(htmlwidgets) # to save the widget!

# load the dataset
# more information about this dataset can be found here:
# https://www.kaggle.com/smithaachar/nyc-dog-licensing-clean
nyc_dogs <- read.csv("data/nyc_dogs.csv")

# convert birth year to date format (and keep only the year)
nyc_dogs$AnimalBirthYear <- mdy_hms(nyc_dogs$AnimalBirthMonth) %>% year()

# identify 10 most common dogs
topdogs <- nyc_dogs %>% count(BreedName) 
topdogs <- topdogs[order(topdogs$n, decreasing = TRUE),]
# keep 10 most common breeds (and remove last year of data which is incomplete)
df <- filter(nyc_dogs, BreedName %in% topdogs$BreedName[2:11] & AnimalBirthYear < 2016) %>% 
  group_by(AnimalBirthYear) %>% 
  count(BreedName) %>% ungroup()

# get some nice colours from viridis (magma)
cols <- viridis::viridis_pal(option = "magma")(length(unique(df$BreedName)))

# make streamgraph!
pp <- streamgraph(df, 
                  key = BreedName, value = n, date = AnimalBirthYear, 
                  height="600px", width="1000px") %>%
  sg_legend(show=TRUE, label="names: ") %>%
  sg_fill_manual(values = cols) 
saveWidget(pp, file=paste0(getwd(), "_posts/2020-09-21-data-visualization/figures/dogs_streamgraph.html"))
