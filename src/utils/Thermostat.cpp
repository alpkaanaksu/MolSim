//
// Created by Feryal Ezgi on 9.12.2023.
//

#include "Thermostat.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include "ArrayUtils.h"

Thermostat::Thermostat() : targetTemperature(0.0),
                           maxTemperatureChange(0.0),
                           thermostatInterval(0),
                           numDimensions(-1),
                           initializeWithBrownianMotion(false) {
}

Thermostat::Thermostat(double initialTemperature, size_t thermostatInterval, int numDimensions,
                       bool initializeWithBrownianMotion)
        : initialTemperature(initialTemperature),
          targetTemperature(initialTemperature),
          maxTemperatureChange(std::numeric_limits<double>::infinity()),
          thermostatInterval(thermostatInterval),
          numDimensions(numDimensions),
          initializeWithBrownianMotion(initializeWithBrownianMotion) {
}


void Thermostat::scaleVelocities(ParticleContainer &particleContainer) {
    double currentTemperature = getCurrentTemperature(particleContainer);
    double temperatureChange = std::copysign(
            std::min(std::abs(targetTemperature - currentTemperature), maxTemperatureChange),
            targetTemperature - currentTemperature);

    double newTemperature = currentTemperature + temperatureChange;
    double scalingFactor = std::sqrt(newTemperature / currentTemperature);

    // Check for zero temperature to avoid division by zero
    if (currentTemperature == 0.0 && (std::isnan(scalingFactor) || std::isinf(scalingFactor))) {
        spdlog::warn("Current temperature is zero. Velocity scaling is skipped.");
        return;
    }

    particleContainer.applyToAll([&scalingFactor](Particle &p) {
        if(!p.isFixed()){
        p.setV(scalingFactor * p.getV());
        }
    });
}


void Thermostat::scaleVelocitiesWithAvg(ParticleContainer &particleContainer) {
    // Calculate average velocity
    std::array<double, 3> totalVelocity = {0.0, 0.0, 0.0};
    // Amount of non-fixed particles
    int N = 0;
    particleContainer.applyToAll([&totalVelocity, &N](Particle &particle) {
        if (!particle.isFixed()) {
            totalVelocity = totalVelocity + particle.getV();
            N++;
        }
    });

    std::array<double, 3> averageVelocity = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        averageVelocity[i] = totalVelocity[i] / N;
    }

    // Scale non-fixed particle velocities
    double currentTemperature = getCurrentTemperature(particleContainer);
    double temperatureChange = std::copysign(
            std::min(std::abs(targetTemperature - currentTemperature), maxTemperatureChange),
            targetTemperature - currentTemperature);

    double newTemperature = currentTemperature + temperatureChange;
    double scalingFactor = std::sqrt(newTemperature / currentTemperature);

    // Check for zero temperature to avoid division by zero
    if (currentTemperature == 0.0 && (std::isnan(scalingFactor) || std::isinf(scalingFactor))) {
        spdlog::warn("Current temperature is zero. Velocity scaling is skipped.");
        return;
    }

    particleContainer.applyToAll([&scalingFactor, &averageVelocity](Particle &p) {
        if (!p.isFixed()) {
            p.setV(averageVelocity + scalingFactor * (p.getV() - averageVelocity));
        }
    });

}


double Thermostat::getCurrentTemperature(ParticleContainer &particleContainer) const {
    double kineticEnergy = 0.0;
    particleContainer.applyToAll([&kineticEnergy](Particle &particle) {
        if (!particle.isFixed()){
            double vSquared = particle.getV()[0] * particle.getV()[0] + particle.getV()[1] * particle.getV()[1] +
                              particle.getV()[2] * particle.getV()[2];

            kineticEnergy = kineticEnergy + (0.5 * particle.getM() * vSquared);
        }

    });

    // We assume everything to be dimensionless, therefore kB = 1.
    return 2 * kineticEnergy / (particleContainer.size() * numDimensions);
}


void Thermostat::initializeTemperature(ParticleContainer &particleContainer) {
    particleContainer.applyToAll([this](Particle &p) {
        auto v = p.getV();
        p.setV(p.getV() +
               maxwellBoltzmannDistributedVelocity(std::sqrt(initialTemperature / p.getM()), numDimensions));

    });

}

double Thermostat::getInitialTemperature() const {
    return initialTemperature;
}

double Thermostat::getTargetTemperature() const {
    return targetTemperature;
}

double Thermostat::getMaxTemperatureChange() const {
    return maxTemperatureChange;
}

size_t Thermostat::getThermostatInterval() const {
    return thermostatInterval;
}

int Thermostat::getNumDimensions() const {
    return numDimensions;
}

bool Thermostat::isInitializeWithBrownianMotion() const {
    return initializeWithBrownianMotion;
}

void Thermostat::setInitialTemperature(double initialTemperature) {
    Thermostat::initialTemperature = initialTemperature;
}

void Thermostat::setTargetTemperature(double targetTemperature) {
    Thermostat::targetTemperature = targetTemperature;
}

void Thermostat::setMaxTemperatureChange(double maxTemperatureChange) {
    Thermostat::maxTemperatureChange = maxTemperatureChange;
}

void Thermostat::setThermostatInterval(size_t thermostatInterval) {
    Thermostat::thermostatInterval = thermostatInterval;
}

void Thermostat::setNumDimensions(int numDimensions) {
    Thermostat::numDimensions = numDimensions;
}

void Thermostat::setInitializeWithBrownianMotion(bool initializeWithBrownianMotion) {
    Thermostat::initializeWithBrownianMotion = initializeWithBrownianMotion;
}
