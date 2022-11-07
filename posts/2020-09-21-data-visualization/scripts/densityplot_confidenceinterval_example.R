# Script to generate 95% confidence intervals of a generated random normal distribution
# as an example in Day 2: Visualizing uncertainty.

# load library
library(ggplot2)
library(magrittr)
library(dplyr)

# set ggplot theme
theme_set(theme_classic() +
            theme(axis.title = element_text(size = 11, face = "bold"),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 13, face = "bold"),
                  legend.title = element_text(size = 11, face = "bold"),
                  legend.text = element_text(size = 10)))

# set random seed
set.seed(22)

# generate population (random normal distribution)
df <- data.frame("value" = rnorm(50, mean = 0, sd = 1))

# descriptive stats for each distribution
desc_stats = df %>% 
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val = sqrt(var(value)/length(value)))

# build density plot!
p <- ggplot(data = df, aes(x = value, y = ..count..)) +
  geom_density(alpha = .2, lwd = .3) +
  xlim(c(min(df$value-1), max(df$value+1))) 
# extract plotted values
base_p <- ggplot_build(p)$data[[1]]
# shade the 95% confidence interval
p + 
  geom_area(data = subset(base_p, 
                          between(x, 
                                  left = (desc_stats$mean_val - 1.96*desc_stats$se_val),
                                  right = (desc_stats$mean_val + 1.96*desc_stats$se_val))),
            aes(x = x, y = y), fill = "cadetblue3", alpha = .6) +
  # add vertical line to show population mean
  geom_vline(aes(xintercept = 0), lty = 2) +
  annotate("text", x = 0.9, y = 19, label = "True mean", fontface = "italic") +
  # label axis!
  labs(x = "Variable of interest", y = "") 
ggsave("figures/confidenceinterval_example.png", width = 5, height = 3.5, units = "in")

