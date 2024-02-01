library(ggplot2)
library(dplyr)

g3 <- read.csv("flowg3.csv") |> mutate(g = 3)
g6 <- read.csv("flowg6.csv") |> mutate(g = 6)
g9 <- read.csv("flowg9.csv") |> mutate(g = 9)

data <- bind_rows(g3, g6, g9)


processed_data <- data |>
  mutate(iteration = as.factor(iteration), g = as.factor(g)) |>
  filter(iteration == 40000)

processed_data |>
  ggplot(aes(bin, avg_velocity_y, color = g)) +
  geom_point() +
  geom_line() +
  labs(x="Bin Position (X-Axis)", y="Average Y-Velocity", title="Average Y-Velocity After 20 Time Units")
 