//
// Created by Feryal Ezgi on 9.12.2023.
//
#pragma once

#include <limits>
#include "../models/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"

class Thermostat {
private:
    double initialTemperature;
    double targetTemperature;
    double maxTemperatureChange;
    size_t thermostatInterval;
    int numDimensions;
    bool initializeWithBrownianMotion;

public:
    Thermostat();

    // Minimum requirement for the thermostat, ∆T = ∞ and targetTemperature = initialTemperature
    Thermostat(double initialTemperature, size_t thermostatInterval, int numDimensions,
               bool initializeWithBrownianMotion);

    void initializeTemperature(ParticleContainer &particleContainer);

    // Function to scale velocities directly or gradually
    void scaleVelocities(ParticleContainer &particleContainer);

    double getCurrentTemperature(ParticleContainer &particleContainer) const;

    double getInitialTemperature() const;

    double getTargetTemperature() const;

    double getMaxTemperatureChange() const;

    size_t getThermostatInterval() const;

    int getNumDimensions() const;

    bool isInitializeWithBrownianMotion() const;

    void setInitialTemperature(double initialTemperature);

    void setTargetTemperature(double targetTemperature);

    void setMaxTemperatureChange(double maxTemperatureChange);

    void setThermostatInterval(size_t thermostatInterval);

    void setNumDimensions(int numDimensions);

    void setInitializeWithBrownianMotion(bool initializeWithBrownianMotion);
};


