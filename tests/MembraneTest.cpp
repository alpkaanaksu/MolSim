//
// Created by Feryal Ezgi on 1.02.2024.
//
#include <gtest/gtest.h>
#include <random>
#include <cmath>
#include <spdlog/spdlog.h>
#include "../src/models/ParticleContainer.h"
#include "../src/utils/Generator.h"
#include "../src/utils/ArrayUtils.h"

class MembraneTest : public ::testing::Test {
protected:
    void TeardownParticleContainer() {
        for (auto &particle: particles) {
            particleContainer.remove(particle);
        }
        particles.clear();
    }

    std::vector<Particle> particles;
    ParticleContainer particleContainer;

};


TEST_F(MembraneTest, MembraneNeighboursTest) {
    spdlog::info("Starting MembraneTest");

    std::array<double, 3> position = {0.0, 0.0, 0.0};
    std::array<int, 3> size = {3, 3, 1};
    double meshWidth = 1.0;
    std::array<double, 3> velocity = {0., 0., 0.};
    double mass = 1.0;
    int typeId = 1;
    double epsilon = 5.0;
    double sigma = 1.0;
    double avgBondLength = 1.0;
    int stiffnessFactor = 2;
    Generator::membrane(particleContainer, position, size, meshWidth, velocity, mass, typeId, epsilon, sigma, avgBondLength, stiffnessFactor);

    int expectedParticleCount = size[0] * size[1] * size[2];
    EXPECT_EQ(particleContainer.size(), expectedParticleCount);
    particles = particleContainer.getParticles();

    // Iterate over particles and check their neighbors
    for (const Particle &particle: particles) {
    // Check if particles have correct neighbors set by the membrane logic.
    Particle nonConstParticle = particle;
    int particlePosition = nonConstParticle.getId();

    if(particlePosition == 4){
        EXPECT_EQ(nonConstParticle.isDiagonalNeighbor(0), true);
        EXPECT_EQ(nonConstParticle.isDirectNeighbor(1), true);
        EXPECT_EQ(nonConstParticle.isDiagonalNeighbor(2), true);
    }
    }

    TeardownParticleContainer();
    ASSERT_TRUE(particleContainer.size() == 0);
    spdlog::info("MembraneTest completed");
}

