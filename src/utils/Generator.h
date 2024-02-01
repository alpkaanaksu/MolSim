//
// Created by Alp Kaan Aksu on 04.11.23.
// Updated by Feryal Ezgi Aşkın on 19.11.23

#pragma once

#include "models/Particle.h"
#include "models/ParticleContainer.h"

namespace Generator {
    /**
     * @brief Generates a cuboid of particles
     *
     * @param container
     * @param position
     * @param size
     * @param meshWidth
     * @param velocity
     * @param mass
     * @param typeId
     */
    void cuboid(ParticleContainer &container, std::array<double, 3> position, std::array<int, 3> size, double meshWidth,
           std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma, bool fixed);

    /**
     * @brief Generates a sphere of particles
     * @param container
     * @param center
     * @param radius
     * @param meshWidth
     * @param velocity
     * @param mass
     * @param typeId
     */
    void sphere(ParticleContainer &container, std::array<double, 3> center, int radius, double meshWidth,
                std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma, bool fixed);

    /**
     * @brief Generates a disk of particles
     * @note Equivalent to a '2d sphere'
     * @param container
     * @param center
     * @param radius
     * @param meshWidth
     * @param velocity
     * @param mass
     * @param typeId
     */
    void disk(ParticleContainer &container, std::array<double, 3> center, int radius, double meshWidth,  std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma, bool fixed);



    /**
     * @brief Generates a cuboid structure of particles with specified attributes and establishes neighbor relationships.
     *
     * The `membrane` method initializes particles within a cuboid structure based on the provided parameters, including
     * position, velocity, mass, and interaction parameters. It establishes neighbor relationships between particles, distinguishing
     * between direct and diagonal neighbors within the cuboid structure.
     *
     * @param container A reference to the ParticleContainer where generated particles will be added.
     * @param position An array representing the starting position of the cuboid structure.
     * @param size An array specifying the dimensions (x, y, z) of the cuboid structure in terms of particles.
     * @param meshWidth The distance between neighboring particles in the cuboid structure.
     * @param velocity An array representing the initial velocity of the particles.
     * @param mass The mass of each particle.
     * @param typeId The type identifier for particles.
     * @param epsilon The interaction parameter epsilon for Lennard-Jones potential.
     * @param sigma The interaction parameter sigma for Lennard-Jones potential.
     * @param avgBondLength The average bond length between particles.
     * @param stiffnessFactor The stiffness factor for particle interactions.
     *
     * The cuboid structure is defined by iterating over the specified size dimensions and creating particles at
     * calculated positions. Neighbors are determined based on the cuboid structure, with direct and diagonal relationships
     * established. Particle information, including positions and neighbor positions, is printed for debugging purposes.
     */

    void membrane(ParticleContainer &container, std::array<double, 3> position, std::array<int, 3> size,
                       double meshWidth, std::array<double, 3> velocity, double mass, int typeId, double epsilon, double sigma, double avgBondLength, int stiffnessFactor);

}
