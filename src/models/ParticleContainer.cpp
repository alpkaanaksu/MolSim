//
// Created by Alp Kaan Aksu, Berke Saylan on 29.10.23.
//

#include <sstream>
#include "ParticleContainer.h"

#include "utils/Generator.h"
#include "io/reader/JSONReader.h"

#include <spdlog/spdlog.h>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

ParticleContainer::ParticleContainer() {
    particles = std::vector<Particle>();

}

void ParticleContainer::applyToAll(const std::function<void(Particle &)> &function) {
    for (auto &p: particles) {
        function(p);
    }
}

void ParticleContainer::applyToAllPairs(const std::function<void(Particle &, Particle &)> &function) {
    for (auto &p1: particles) {
        for (auto &p2: particles) {
            if (p1 == p2) {
                // same particle. skip.
                continue;
            }

            function(p1, p2);
        }
    }
}

void ParticleContainer::applyToAllPairsOnce(const std::function<void(Particle &, Particle &)> &function) {
    for (int i = 0; i < particles.size(); i++) {
        for (int j = i + 1; j < particles.size(); j++) {
            function(particles[i], particles[j]);
        }
    }
}


void ParticleContainer::add(const Particle &particle) {
    particles.push_back(particle);
}

void ParticleContainer::add(const nlohmann::json &objects) {
    for (auto &object: objects) {
        if (object["type"] == "particle") {
            std::array<double, 3> f = {0., 0., 0.};
            std::array<double, 3> old_f = {0., 0., 0.};

            if (object.contains("force")) {
                f = object["force"];
            }

            if (object.contains("old_force")) {
                old_f = object["old_force"];
            }

            add(Particle{object["position"], object["velocity"], f, old_f, object["mass"], object["epsilon"], object["sigma"], object["type_id"], 0.0, 0});
        } else if (object["type"] == "cuboid") {
            Generator::cuboid(*this, object["position"], object["size"], object["mesh_width"], object["velocity"],
                              object["mass"], object["type_id"], object["epsilon"], object["sigma"]);
        } else if (object["type"] == "sphere") {
            Generator::sphere(*this, object["center"], object["radius"], object["mesh_width"], object["velocity"],
                              object["mass"], object["type_id"], object["epsilon"], object["sigma"]);
        } else if (object["type"] == "disk") {
            Generator::disk(*this, object["center"], object["radius"], object["mesh_width"], object["velocity"],
                              object["mass"], object["type_id"], object["epsilon"], object["sigma"]);
        } else if (object["type"] == "checkpoint") {
            resolveCheckpoint(object["path"]);
        } else if (object["type"] == "membrane") {
            Generator::membrane(*this, object["position"], object["size"], object["mesh_width"], object["velocity"],
                              object["mass"], object["type_id"], object["epsilon"], object["sigma"], object["avg_bond_length"], object["stiffness_factor"]);
        }
    }
}

void ParticleContainer::resolveCheckpoint(const std::string &path) {
    // Get the absolute path for source
    fs::path absolutePath = fs::absolute(source);

    // Append the referencedFile to the directory of path
    fs::path resolvedPath = absolutePath.parent_path() / path;

    // Check if the resolved path exists
    if (fs::exists(resolvedPath)) {
        spdlog::info("Resolved path for " + path + ": " + resolvedPath.string());
    } else {
        spdlog::error("File " + path + " does not exist relative to " + source);
        exit(-1);
    }

    auto j = JSONReader::readFile(resolvedPath.string());
    add(j);
}

int ParticleContainer::size() {
    return static_cast<int>(particles.size());
}

void ParticleContainer::remove(Particle &particle) {
    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i] == particle) {
            particles.erase(particles.begin() + i);
            break;
        }
    }
}

std::vector<Particle> &ParticleContainer::getParticles() {
    return particles;
}

std::string ParticleContainer::toString() {
    std::stringstream stream;
    stream << "\n--- Container ---------"
           << "\nType: naive"
           << "\nNumber of particles: " << size()
           << "\n------------------------\n";

    return stream.str();
}

nlohmann::ordered_json ParticleContainer::json() {
    nlohmann::ordered_json j;
   for (auto &p : particles) {
         j.push_back(p.json());
    }

   return j;
}

void ParticleContainer::setSource(const std::string src) {
    this->source = src;
}

std::string ParticleContainer::getSource() {
    return this->source;
}

std::ostream &operator<<(std::ostream &stream, ParticleContainer &simulation) {
    stream << simulation.toString();
    return stream;
}

void ParticleContainer::applyToAllHalo(const std::function<void(Particle &)> &function) {
    applyToAll(function);
}
