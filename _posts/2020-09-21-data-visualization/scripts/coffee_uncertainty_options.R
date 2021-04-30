# Script to generate plots with various ways of representing uncertainty, based 
# Coffee & Code dataset from https://www.kaggle.com/devready/coffee-and-code/data

# set-up & data manipulation ---------------------------------------------------

# load packages
library(ggplot2) # for plots, built layer by layer
library(dplyr) # for data manipulation
library(magrittr) # for piping
library(ggridges) # for density ridge plots
library(patchwork) # great package for "patching" plots together!

# set ggplot theme
theme_set(theme_classic() +
            theme(axis.title = element_text(size = 11, face = "bold"),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 13, face = "bold"),
                  legend.title = element_text(size = 11, face = "bold"),
                  legend.text = element_text(size = 10)))

# import data
df <- read.csv("data/coffee_code.csv")

# set labels to be used in all plots
coffee_labels <- labs(title = "Does coffee help programmers code?",
                      x = "Coffee while coding", 
                      y = "Time spent coding \n(hours/day)") 

# the variable CodingWithoutCoffee is negative, which is harder to understand
# (i.e. "No" means they drink coffee...). So, let's transform it into a more 
# intuitive variable!
df$CodingWithCoffee <- gsub("No", "Usually", df$CodingWithoutCoffee)
df$CodingWithCoffee <- gsub("Yes", "Rarely\n or never", df$CodingWithCoffee)
# convert to factor and set levels so they show up in a logical order
df$CodingWithCoffee <- factor(df$CodingWithCoffee,
                              levels = c("Rarely\n or never", 
                                         "Sometimes", 
                                         "Usually"))

# calculate summary statistics for the variable of choice
df_summary <- group_by(df, CodingWithCoffee) %>%
  summarise(
    # mean
    mean_codinghours = mean(CodingHours), 
    # standard deviation
    sd_codinghours = sd(CodingHours), 
    # standard error
    se_codinghours = sd(CodingHours)/sqrt(length(CodingHours)))


# 1. Error bars (standard error) -----------------------------------------------

ggplot(df_summary) +
  geom_errorbar(aes(x = CodingWithCoffee, 
                    ymin = mean_codinghours - se_codinghours,
                    ymax = mean_codinghours + se_codinghours), 
                width = .2) +
  geom_point(aes(x = CodingWithCoffee, y = mean_codinghours), 
             size = 3) +
  coffee_labels + ylim(0,10)
ggsave("figures/coffee_errorbars.png", width = 5, height = 3, units = "in")

# 2. Boxplot -------------------------------------------------------------------

ggplot(df) +
  geom_boxplot(aes(x = CodingWithCoffee, y = CodingHours)) +
  coffee_labels
ggsave("figures/coffee_boxplot.png", width = 5, height = 3, units = "in")


# 3. Error bar demonstration ---------------------------------------------------

# get some nice colours from viridis (magma)
error_cols <- viridis::viridis_pal(option = "magma")(5)[2:4]
# set labels to be used in the palette
error_labels = c("standard deviation","95% confidence interval","standard error")

ggplot(df_summary) +
  # show the raw data
  geom_jitter(data = df, aes(x = CodingWithCoffee, 
                             y = CodingHours),
              alpha = .5, width = .05, col = "grey") +
  # standard deviation
  geom_errorbar(aes(x = CodingWithCoffee, 
                    ymin = mean_codinghours - sd_codinghours,
                    ymax = mean_codinghours + sd_codinghours,
                    col = "SD"), width = .2, lwd = 1) +
  # 95% confidence interval
  geom_errorbar(aes(x = CodingWithCoffee, 
                    ymin = mean_codinghours - 1.96*se_codinghours,
                    ymax = mean_codinghours + 1.96*se_codinghours, 
                    col = "CI"), width = .2, lwd = 1) +
  # standard error
  geom_errorbar(aes(x = CodingWithCoffee, 
                    ymin = mean_codinghours - se_codinghours,
                    ymax = mean_codinghours + se_codinghours, 
                    col = "SE"), width = .2, lwd = 1) +
  geom_point(aes(x = CodingWithCoffee, y = mean_codinghours), 
             size = 3) +
  coffee_labels + ylim(c(0,11)) +
  # manual palette/legend set-up!
  scale_colour_manual(name = "Uncertainty metric", 
                      values = c(SD = error_cols[1], 
                                 CI = error_cols[2], 
                                 SE = error_cols[3]),
                      labels = error_labels) +
  theme(legend.position = "top")
ggsave("figures/coffee_bars_demo.png", width = 7, height = 5, units = "in")


# 4. Jitter plot with violin ---------------------------------------------------

ggplot() +
  geom_jitter(data = df, aes(x = CodingWithCoffee, 
                             y = CodingHours),
              alpha = .5, width = .05, col = "grey") +
  geom_violin(data = df, aes(x = CodingWithCoffee, 
                             y = CodingHours), alpha = 0) +
  geom_linerange(data = df_summary,
                 aes(x = CodingWithCoffee, 
                     ymin = mean_codinghours - se_codinghours,
                     ymax = mean_codinghours + se_codinghours)) +
  geom_point(data = df_summary, 
             aes(x = CodingWithCoffee, 
                 y = mean_codinghours), size = 3) +
  coffee_labels
ggsave("figures/coffee_violin_jitter.png", width = 5, height = 3, units = "in")


# 5. Density ridge plot --------------------------------------------------------

ggplot(df) + 
  aes(y = CodingWithCoffee, x = CodingHours, fill = stat(x)) +
  geom_density_ridges_gradient(scale = 1.9, size = .2, rel_min_height = 0.005) +
  # colour palette (gradient according to CodingHours)
  scale_fill_viridis_c(option = "magma", direction = -1) +
  # remove legend - it's not necessary here!
  theme(legend.position = "none") +
  labs(title = coffee_labels$title, 
       x = coffee_labels$y, 
       y = "Coffee \nwhile coding") + 
  theme(axis.title.y = element_text(angle=0, hjust = 1, vjust = .9, 
                                    margin = margin(t = 0, r = -50, b = 0, l = 0)))
ggsave("figures/coffee_density_ridges.png", width = 5, height = 3, units = "in")

# 6. Jitter vs. Rug plot ------------------------------------------------------------------

jitterplot <- ggplot(df, aes(x = CoffeeCupsPerDay, y = CodingHours)) +
  geom_jitter(alpha = .2) +
  geom_smooth(fill = error_cols[1], col = "black", method = lm, lwd = .7) +
  coffee_labels + ylim(c(0,11)) + labs(x = "Cups of coffee (per day)")

rugplot <- ggplot(df, aes(x = CoffeeCupsPerDay, y = CodingHours)) +
  geom_smooth(fill = error_cols[1], col = "black", method = lm, lwd = .7) +
  geom_rug(position="jitter", alpha = .7) + ylim(c(0,11)) +
  coffee_labels + labs(x = "Cups of coffee (per day)")

# patch the two plots together
jitterplot + rugplot
ggsave("figures/coffee_jitter_vs_rugplot.png", width = 10, height = 4, units = "in")



