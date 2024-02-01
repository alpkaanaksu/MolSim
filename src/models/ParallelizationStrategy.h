//
// Created by Alp Kaan Aksu on 01.01.24.
//
#pragma once

#include <string>

/**
 * @brief Enum for parallelization strategy
 */
enum class ParallelizationStrategy {
    Cells,
    Neighbors,
    None
};

/**
 * @brief Converts parallelization strategy to string
 *
 * @param parallelizationStrategy
 * @return String representation of parallelization strategy
 */
std::string parallelizationStrategyToString(ParallelizationStrategy parallelizationStrategy);

/**
 * @brief Converts string to parallelization strategy
 *
 * @param parallelizationStrategy
 * @return ParallelizationStrategy represemted by a string
 */
ParallelizationStrategy stringToParallelizationStrategy(const std::string &parallelizationStrategy);