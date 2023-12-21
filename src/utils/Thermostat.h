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
    /**
     * @brief Default constructor (thermostat off)
     *
     * @note Thermostat is not applied, if it wasn't initialized correctly, this is done by checking the value of numDimensions.
     */
    Thermostat();

    /**
     * @brief Initializes the thermostat with the given initial temperature and thermostat interval
     *
     * @param initialTemperature The initial temperature of the particle container
     * @param thermostatInterval The interval at which the thermostat is applied
     * @param numDimensions The number of dimensions of the particle container
     * @param initializeWithBrownianMotion If true, the velocities of the particles are initialized with a Maxwell-Boltzmann distribution
     *
     * @note Minimum requirement for the thermostat, ∆T = ∞ and targetTemperature = initialTemperature
     */
    Thermostat(double initialTemperature, size_t thermostatInterval, int numDimensions,
               bool initializeWithBrownianMotion);

    /**
     * @brief Initializes the temperature of the particle container according to the initial temperature
     *
     * @param particleContainer The particle container to initialize the temperature of
     *
     * @note if initializeWithBrownianMotion (see constructor) is true, the velocities of the particles are initialized with a Maxwell-Boltzmann distribution
     */
    void initializeTemperature(ParticleContainer &particleContainer);

    /**
     * @brief Scales the velocities of all particles in the particle container to match the target temperature
     *
     * @param particleContainer The particle container to scale the velocities of
     */
    void scaleVelocities(ParticleContainer &particleContainer);

    /**
     * @brief Calculates the current temperature of the particle container according to the kinetic energy of the particles
     *
     * @param particleContainer The particle container to calculate the temperature of
     * @return The current temperature of the particle container
     */
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


