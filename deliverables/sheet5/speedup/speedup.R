library(ggplot2)
library(dplyr)

data <- read.csv("thread.csv")

base <- data |> filter(threads == 1)
base_time <- base$time
base_mups <- base$mups

processed_data <- data |>
  mutate(speedup = base_time / time)

processed_data |>
  ggplot(aes(threads, speedup)) +
  geom_point() +
  geom_line() +
  labs(
    x = "#Threads",
    y = "Speedup T(1)/T(#Threads)",
    title = "Speedup with Threads"
  )