//
// Created by Alp Kaan Aksu on 01.11.23.
//

#include <iostream>

#include "nlohmann/json.hpp"

#include "Simulation.h"
#include "io/reader/FileReader.h"
#include "io/reader/JSONReader.h"
#include "io/outputWriter/JSONWriter.h"
#include "io/outputWriter/Writer.h"
#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include <spdlog/spdlog.h>
#include <chrono>

#include "models/LinkedCellParticleContainer.h"

using json = nlohmann::json;

Simulation::Simulation(const std::string &filepath) {
    json definition = JSONReader::readFile(filepath);

    if (definition["simulation"]["particle_container"]["type"] == "basic") {
        particles = std::make_shared<ParticleContainer>();
    } else if (definition["simulation"]["particle_container"]["type"] == "linked_cell") {
        if (definition["simulation"]["particle_container"].contains("boundary")) {
            auto top = BoundaryBehavior::Reflective;
            auto bottom = BoundaryBehavior::Reflective;
            auto left = BoundaryBehavior::Reflective;
            auto right = BoundaryBehavior::Reflective;
            auto front = BoundaryBehavior::Reflective;
            auto back = BoundaryBehavior::Reflective;

            if (definition["simulation"]["particle_container"]["boundary"].contains("all")) {
                auto behavior = stringToBoundaryBehavior(
                        definition["simulation"]["particle_container"]["boundary"]["all"]);
                top = behavior;
                bottom = behavior;
                left = behavior;
                right = behavior;
                front = behavior;
                back = behavior;
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("top")) {
                top = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["top"]);
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("bottom")) {
                bottom = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["bottom"]);
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("left")) {
                left = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["left"]);
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("right")) {
                right = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["right"]);
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("front")) {
                front = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["front"]);
            }

            if (definition["simulation"]["particle_container"]["boundary"].contains("back")) {
                back = stringToBoundaryBehavior(definition["simulation"]["particle_container"]["boundary"]["back"]);
            }

            particles = std::make_shared<LinkedCellParticleContainer>(
                    definition["simulation"]["particle_container"]["dimensions"][0],
                    definition["simulation"]["particle_container"]["dimensions"][1],
                    definition["simulation"]["particle_container"]["dimensions"][2],
                    definition["simulation"]["particle_container"]["cutoff_radius"],
                    definition["simulation"]["time_delta"],
                    top,
                    bottom,
                    right,
                    left,
                    front,
                    back
            );

        } else {
            particles = std::make_shared<LinkedCellParticleContainer>(
                    definition["simulation"]["particle_container"]["dimensions"][0],
                    definition["simulation"]["particle_container"]["dimensions"][1],
                    definition["simulation"]["particle_container"]["dimensions"][2],
                    definition["simulation"]["particle_container"]["cell_size"],
                    definition["simulation"]["time_delta"]
            );
        }
    }

    particles->setSource(filepath);

    endTime = definition["simulation"]["end_time"];
    deltaT = definition["simulation"]["time_delta"];
    videoDuration = definition["simulation"]["video_duration"];
    fps = definition["simulation"]["frame_rate"];
    out = definition["simulation"]["output_path"];
    in = filepath;
    outputType = outputWriter::stringToOutputType(definition["simulation"]["output_type"]);
    gravity = definition["simulation"].contains("gravity")
              ? (double) definition["simulation"]["gravity"]
              : 0.0;

    particles->add(definition["objects"]);

    if (definition["simulation"]["model"] == "basic") {
        model = Model::gravityModel(deltaT);
    } else if (definition["simulation"]["model"] == "lennard_jones") {
        model = Model::lennardJonesModel(deltaT);
    } else if (definition["simulation"]["model"] == "membrane") {
        model = Model::membraneModel(deltaT);
        membrane = true;
    }

    if (definition["simulation"].contains("thermostat")) {
        thermostat = Thermostat(definition["simulation"]["thermostat"]["initial_temperature"],
                                definition["simulation"]["thermostat"]["interval"],
                                definition["simulation"]["thermostat"]["dimension"],
                                definition["simulation"]["thermostat"]["brownian"]);

        if (definition["simulation"]["thermostat"].contains("target_temperature") &&
            definition["simulation"]["thermostat"].contains("delta_temperature")) {
            // Target temperature and ∆T are specified
            thermostat.setMaxTemperatureChange(definition["simulation"]["thermostat"]["delta_temperature"]);
            thermostat.setTargetTemperature(definition["simulation"]["thermostat"]["target_temperature"]);

        } else if (definition["simulation"]["thermostat"].contains("target_temperature") &&
                   !definition["simulation"]["thermostat"].contains("delta_temperature")) {
            // Target temperature specified and ∆T not
            thermostat.setTargetTemperature(definition["simulation"]["thermostat"]["target_temperature"]);

        } else if (!definition["simulation"]["thermostat"].contains("target_temperature") &&
                   definition["simulation"]["thermostat"].contains("delta_temperature")) {
            // ∆T specified and target temperature not
            thermostat.setMaxTemperatureChange(definition["simulation"]["thermostat"]["delta_temperature"]);
        } else {
            // ∆T = ∞ and targetTemperature = initialTemperature per default
            thermostat = Thermostat(definition["simulation"]["thermostat"]["initial_temperature"],
                                    definition["simulation"]["thermostat"]["interval"],
                                    definition["simulation"]["thermostat"]["dimension"],
                                    definition["simulation"]["thermostat"]["brownian"]);
        }

    }

    if(definition["simulation"].contains("checkpoints")){
        for(auto &checkpoint : definition["simulation"]["checkpoints"]){
            checkpoints.push(checkpoint);
        }
    }

    checkpoints.push(-1);
}

Simulation::Simulation(Model model, double endTime, double deltaT, int videoDuration, int fps, const std::string &in,
                       std::string out, outputWriter::OutputType outputType)
        : endTime(endTime), deltaT(deltaT), videoDuration(videoDuration), fps(fps), in(in), out(std::move(out)),
          model(std::move(model)), outputType(outputType) {

    FileReader::readFile(*particles, in);
}

void Simulation::run() {
    auto linkedCellParticleContainer = dynamic_cast<LinkedCellParticleContainer *>(particles.get());

    outputWriter::prepareOutputFolder(out);

    double current_time = 0;

    int iteration = 0;

    int plotInterval = static_cast<int>((endTime / deltaT) / (videoDuration * fps));

    // Avoid division by zero
    if (plotInterval == 0) {
        plotInterval = 30;
    }

    std::array<double, 3> pullingForce = {0.0, 0.8, 0.0};
    auto resetForce = Model::resetForceFunction();
    auto force = model.forceFunction();
    auto position = model.positionFunction();
    auto velocity = model.velocityFunction();

    // Calculate initial force to avoid starting with 0 force
    if (membrane && linkedCellParticleContainer != nullptr) {
        linkedCellParticleContainer->applyToAllPairsOnceMembrane(force);
    } else {
        particles->applyToAllPairsOnce(force);
    }

    // Brownian Motion with scaling factor
    if (thermostat.getNumDimensions() != -1 && thermostat.isInitializeWithBrownianMotion()) {
        thermostat.initializeTemperature(*particles);
    }

    double nextCheckpoint = checkpoints.front();
    checkpoints.pop();

    auto before = std::chrono::high_resolution_clock::now();

    long numberOfUpdates { 0 };

    // for this loop, we assume: current x, current f and current v are known
    while (current_time <= endTime) {
        // calculate new x
        // Try to cast to LinkedCellParticleContainer

        if (linkedCellParticleContainer != nullptr) {
            // particles points to a LinkedCellParticleContainer
            linkedCellParticleContainer->applyToAll(position, true);

        } else {
            particles->applyToAll(position);
        }

        // calculate new f
        particles->applyToAll([&resetForce, this, &numberOfUpdates, &pullingForce, &current_time, &linkedCellParticleContainer](Particle &p) {
            resetForce(p);
            p.setF(p.getF() + Model::verticalGravityForce(p.getM(), gravity));

            // Pull
            if(current_time < 150 && p.isPulled()) {
                p.setF(p.getF() + pullingForce);
            }

            numberOfUpdates++;
        });

        if (membrane && linkedCellParticleContainer != nullptr) {
            linkedCellParticleContainer->applyToAllPairsOnceMembrane(force);
        } else {
            linkedCellParticleContainer->applyToAllPairsOnce(force);
        }
        
    

        // calculate new v
        particles->applyToAll([&velocity, &numberOfUpdates](Particle &p){
            velocity(p);
        });

        iteration++;

        if (thermostat.getNumDimensions() != -1 && iteration % thermostat.getThermostatInterval() == 0) {
            thermostat.scaleVelocities(*particles);
        }


        if (nextCheckpoint != -1 && current_time >= nextCheckpoint) {
            std::ostringstream oss;
            oss << std::setprecision(2) <<  nextCheckpoint;
            std::string checkpointStr = oss.str();

            spdlog::info("Checkpoint " + checkpointStr + " reached. Saving simulation to file.");
            JSONWriter::writeFile(particles->json(), out + "/checkpoint_" + checkpointStr + ".cp.json");
            spdlog::info("Simulation saved.");
            nextCheckpoint = checkpoints.front();
            checkpoints.pop();
        }

        if (iteration % plotInterval == 0) {
            plotParticles(iteration);
        }

        double percentage = current_time / endTime * 100;

        std::cout << std::fixed << std::setprecision(2) << "Running simulation: [ " << current_time << " (" << current_time / endTime * 100
                  << "%) ] " << "\r" << std::flush;

        current_time += deltaT;
    }

    auto after = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(after - before);

    std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
    std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration - sec);

    spdlog::info("Running simulation: [ 100% ] Done.\n");

    spdlog::info("Saving end state of the simulation to file.\n");
    JSONWriter::writeFile(particles->json(), out + "/checkpoint_end.cp.json");

    // MUP/s
    double seconds = duration.count() / 1e6;
    double updatesPerSecond = numberOfUpdates / seconds;

    spdlog::info("Time: " + std::to_string(sec.count()) + "." + std::to_string(ms.count()));
    spdlog::info("MUP/s: " + std::to_string(updatesPerSecond));
}

void Simulation::plotParticles(int iteration) {
    outputWriter::VTKWriter vtkWriter{};
    outputWriter::XYZWriter xyzWriter{};

    outputWriter::Writer *writer = &vtkWriter;

    switch (outputType) {
        case outputWriter::VTK: {
            writer = &vtkWriter;
            break;
        }
        case outputWriter::XYZ: {
            writer = &xyzWriter;
            break;
        }
        case outputWriter::Disabled: {
            writer = nullptr;
            break;
        }
    }

    std::string out_name(out + "/MD");

    if (writer != nullptr) {
        writer->plotParticles(*particles, out_name, iteration);
    }
}

std::string Simulation::toString() const {
    std::stringstream stream;
    stream << "\n====== Simulation ======"
           << "\nEnd time: " << endTime
           << "\nTime delta: " << deltaT
           << "\nGravity: " << gravity
           << "\nVideo duration (s): " << videoDuration
           << "\nFrames per second: " << fps
           << "\n"
           << "\nReading from: " << in
           << "\nOutput to: " << out << '/'
           << "\nOutput type: " << outputWriter::outputTypeToString(outputType)
           << "\n" << particles->toString()
           << "\n========================\n";

    return stream.str();
}

double Simulation::getEndTime() const {
    return endTime;
}

double Simulation::getDeltaT() const {
    return deltaT;
}

int Simulation::getVideoDuration() const {
    return videoDuration;
}

int Simulation::getFPS() const {
    return fps;
}

std::string Simulation::getInputFilePath() const {
    return in;
}

std::string Simulation::getOutputPath() const {
    return out;
}

std::shared_ptr<ParticleContainer> Simulation::getParticles() const {
    return particles;
}

Model Simulation::getModel() const {
    return model;
}

outputWriter::OutputType Simulation::getOutputType() const {
    return outputType;
}

const Thermostat &Simulation::getThermostat() const {
    return thermostat;
}

std::ostream &operator<<(std::ostream &stream, Simulation &simulation) {
    stream << simulation.toString();
    return stream;
}