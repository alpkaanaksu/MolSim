//
// Created by Alp Kaan Aksu on 29.10.23.
//

#include "ParticleContainer.h"

ParticleContainer::ParticleContainer() {
    particles = std::vector<Particle>();

}

void ParticleContainer::applyToAll(const std::function<void(Particle&)>& function) {
    for (auto &p: particles) {
        function(p);
    }
}

void ParticleContainer::applyToAllPairs(const std::function<void(Particle&, Particle&)>& function) {
    for (auto &p1: particles) {
        for (auto &p2 : particles) {
            if (p1 == p2) {
                // same particle. skip.
                continue;
            }

            function(p1, p2);
        }
    }
}

void ParticleContainer::plot(outputWriter::VTKWriter &writer) {
    for (auto &p: particles) {
        writer.plotParticle(p);
    }
}

void ParticleContainer::add(const Particle &particle) {
    particles.push_back(particle);
}

size_t ParticleContainer::count() {
    return particles.size();
}