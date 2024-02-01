//
// Created by Alp Kaan Aksu on 02.12.23.
//

#include "ParallelizationStrategy.h"

std::string parallelizationStrategyToString(ParallelizationStrategy parallelizationStrategy) {
    switch (parallelizationStrategy) {
        case ParallelizationStrategy::Cells:
            return "Cells";
        case ParallelizationStrategy::Neighbors:
            return "Neighbors";
        case ParallelizationStrategy::None:
            return "None";
    }
}

ParallelizationStrategy stringToParallelizationStrategy(const std::string &parallelizationStrategy) {
    if (parallelizationStrategy == "cells") {
        return ParallelizationStrategy::Cells;
    } else if (parallelizationStrategy == "neighbors") {
        return ParallelizationStrategy::Neighbors;
    } else if (parallelizationStrategy == "none") {
        return ParallelizationStrategy::None;
    } else {
        return ParallelizationStrategy::Cells;
    }
}

