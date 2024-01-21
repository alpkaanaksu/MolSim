//
// Created by Feryal Ezgi on 16.11.2023.
//

#include "Generator.h"
#include "ArrayUtils.h"
#include <array>
#include <iostream>
#include <spdlog/spdlog.h>

void Generator::cuboid(ParticleContainer &container, std::array<double, 3> position, std::array<int, 3> size,
                       double meshWidth, std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma) {
    for (int x = 0; x < size[0]; x++) {
        for (int y = 0; y < size[1]; y++) {
            for (int z = 0; z < size[2]; z++) {
                container.add(Particle{{
                                               position[0] + x * meshWidth,
                                               position[1] + y * meshWidth,
                                               position[2] + z * meshWidth
                                       },
                                       velocity,
                                       mass,
                                       epsilon,
                                       sigma,
                                       typeId}
                );
            }
        }
    }
}

void Generator::membrane(ParticleContainer &container, std::array<double, 3> position, std::array<int, 3> size,
                         double meshWidth, std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma, double avgBondLength, int stiffnessFactor) {

    std::vector<Particle*> particleIndexVector;

    int nextId = 0;

    for (int x = 0; x < size[0]; x++) {
        for (int y = 0; y < size[1]; y++) {
            for (int z = 0; z < size[2]; z++) {
                Particle* newParticle = new Particle {
                        {
                                position[0] + x * meshWidth,
                                position[1] + y * meshWidth,
                                position[2] + z * meshWidth
                        },
                        velocity,
                        mass,
                        epsilon,
                        sigma,
                        typeId,
                        avgBondLength,
                        stiffnessFactor,
                        //(x == 2 && y == 2) || (x == 2 && y == 3)
                        (x == 17 && y == 24) || (x == 17 && y == 25) || (x == 18 && y == 24) || (x == 18 && y == 25)
                };

                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        // Skip the current particle
                        if (dx == 0 && dy == 0) {
                            continue;
                        }

                        int nx = nextId % size[1] + dx;
                        int ny = nextId / size[1] + dy;

                        // Check if the neighbor is within the valid range
                        if (nx >= 0 && nx < size[0] && ny >= 0 && ny < size[1]) {
                            int neighborParticleIndex = ny * size[0] + nx;

                            if (neighborParticleIndex >= 0 && neighborParticleIndex < (size[0] * size[1])) {
                                if (dx == 0 || dy == 0) {
                                    newParticle->addDirectNeighbor(neighborParticleIndex);
                                } else {
                                    newParticle->addDiagonalNeighbor(neighborParticleIndex);
                                }
                            }
                        }
                    }
                }

                newParticle->setId(nextId++);
                particleIndexVector.push_back(newParticle);
                container.add(*newParticle);
            }
        }
    }

    // Print out all particles' positions and their neighbors' positions (not ids)
    for (int i = 0; i < particleIndexVector.size(); i++) {
        std::cout << std::endl;
        std::cout << "Particle " << i << " position: " << particleIndexVector[i]->getX()[0] << " " << particleIndexVector[i]->getX()[1] << " " << particleIndexVector[i]->getX()[2] << std::endl;
        std::cout << "Particle " << i << " direct neighbors: " << std::endl;
        for (int j = 0; j < particleIndexVector[i]->getDirectNeighbors().size(); j++) {
            int neighborIndex = particleIndexVector[i]->getDirectNeighbors()[j];
            std::cout << "Particle " << neighborIndex << " position: "
                      << particleIndexVector[neighborIndex]->getX()[0] << " "
                      << particleIndexVector[neighborIndex]->getX()[1] << " "
                      << particleIndexVector[neighborIndex]->getX()[2] << std::endl;
        }
        std::cout << "Particle " << i << " diagonal neighbors: " << std::endl;
        for (int j = 0; j < particleIndexVector[i]->getDiagonalNeighbors().size(); j++) {
            int neighborIndex = particleIndexVector[i]->getDiagonalNeighbors()[j];
            std::cout << "Particle " << neighborIndex << " position: "
                      << particleIndexVector[neighborIndex]->getX()[0] << " "
                      << particleIndexVector[neighborIndex]->getX()[1] << " "
                      << particleIndexVector[neighborIndex]->getX()[2] << std::endl;
        }
    }
}


// Iterate over a cubic area around the sphere with the given parameters and add a particle to container if it is inside the sphere boundaries
void Generator::sphere(ParticleContainer &container, std::array<double, 3> center, int radius, double meshWidth,
                       std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma) {

    // Distance from the center to the edge of the sphere
    double dis = radius * meshWidth;
    // Calculate the bounds for iteration
    double minBoundX = center[0] - dis;
    double maxBoundX = center[0] + dis;
    double minBoundY = center[1] - dis;
    double maxBoundY = center[1] + dis;
    double minBoundZ = center[2] - dis;
    double maxBoundZ = center[2] + dis;

    // Iterate over the cubic area around the sphere and add particles
    for (double x = minBoundX; x <= maxBoundX; x += meshWidth) {
        for (double y = minBoundY; y <= maxBoundY; y += meshWidth) {
            for (double z = minBoundZ; z <= maxBoundZ; z += meshWidth) {
                std::array<double, 3> position = {x, y, z};
                double normalizedDistance = ArrayUtils::L2Norm(position - center) / meshWidth;

                // Check if the particle is inside the sphere boundaries
                if (normalizedDistance <= radius) {
                    container.add(Particle{position, velocity, mass, epsilon, sigma, typeId});
                }
            }
        }
    }
}

// Iterate over a square area around the disk with the given parameters and add a particle to container if it is inside the disk boundaries
void Generator::disk(ParticleContainer &container, std::array<double, 3> center, int radius, double meshWidth,
                       std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma) {

    double dis = radius * meshWidth;
    double minBoundX = center[0] - dis;
    double maxBoundX = center[0] + dis;
    double minBoundY = center[1] - dis;
    double maxBoundY = center[1] + dis;

    for (double x = minBoundX; x <= maxBoundX; x += meshWidth) {
        for (double y = minBoundY; y <= maxBoundY; y += meshWidth) {
            std::array<double, 3> position = {x, y, center[2]};
            double normalizedDistance = ArrayUtils::L2Norm(position - center) / meshWidth;
            if (normalizedDistance <= radius) {
                container.add(Particle{position, velocity, mass, epsilon, sigma,typeId});
            }
        }
    }
}