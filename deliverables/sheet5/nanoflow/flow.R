library(ggplot2)
library(dplyr)

data <- read.csv("flow.csv")

processed_data <- data |>
  mutate(iteration = as.factor(iteration))

processed_data |>
  ggplot(aes(bin, avg_velocity_y, color = iteration)) +
  geom_point() +
  geom_line() +
  labs(
    x = "Bin Position (X-Axis)",
    y = "Average Y-Velocity",
    title = "Flow through a nanotube (Gravity: -3, End time: 20, Time delta: 0.0005)"
  )
