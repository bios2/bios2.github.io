# Script to make animated plot of volcano eruptions over time

# Load libraries:
library(dplyr) # data manipulation
library(ggplot2) # plotting
library(gganimate) # animation
library(gifski) # creating gifs

# set ggplot theme
theme_set(theme_classic() +
            theme(axis.title = element_text(size = 11, face = "bold"),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 13, face = "bold"),
                  legend.title = element_text(size = 11, face = "bold"),
                  legend.text = element_text(size = 10)))

# function to floor a year to the decade
floor_decade = function(value){return(value - value %% 10)}

# read data 
eruptions <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-05-12/eruptions.csv')

# select top 5 most frequently exploding volcanoes
temp <- group_by(eruptions, volcano_name) %>% tally() 
temp <- temp[order(temp$n, decreasing = TRUE),]

# make a time series dataset (number of eruptions per year)
eruptions$start_decade = floor_decade(eruptions$start_year)

# filter dataset to subset we want to visualize
df <- eruptions %>% 
  filter(between(start_decade, 1900, 2019)) %>%
  filter(volcano_name %in% temp$volcano_name[1:5]) %>%
  group_by(start_decade) %>%
  count(volcano_name) %>% ungroup()

# plot!
p <- ggplot(df, aes(x = start_decade, y = n, fill = volcano_name)) +
  geom_area() +
  geom_vline(aes(xintercept = start_decade)) + # line that follows the current decade
  scale_fill_viridis_d(option = "magma", end = .8) +
  labs(x = "", y = "Number of eruptions", fill = "Volcano",
       title = 'Eruptions of the top 5 most frequently erupting volcanos worldwide') +
  # gganimate part: reveals each decade
  transition_reveal(start_decade) 
animate(p, duration = 5, fps = 20, width = 800, height = 300, renderer = gifski_renderer())
anim_save("figures/volcano_eruptions.gif")
