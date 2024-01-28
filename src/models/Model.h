//
// Created by Alp Kaan Aksu on 01.11.23.
//

#pragma once

#include <utility>

#include "functional"
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include <array>

/**
 * Model class
 * Represents a model for a simulation
 */
class Model {
private:
    /**
     * Formula that is used to calculate the forces on particles
     */
    std::function<void(Particle&, Particle&)> force;

    /**
     * Formula that is used to calculate the positions of particles
     */
    std::function<void(Particle&)> position;

    /**
     * Formula that is used to calculate the velocity of particles
     */
    std::function<void(Particle&)> velocity;

public:
    Model(std::function<void(Particle&, Particle&)> force, std::function<void(Particle&)> position, std::function<void(Particle&)> velocity) : force(std::move(force)), position(std::move(position)), velocity(std::move(velocity)) {};

    Model() : force([](Particle &p1, Particle &p2){}), position([](Particle &p){}),  velocity([](Particle &p){}) {};


    /**
     * @brief Returns a function that takes a particle as a parameter and sets the force values of that particle along all axis to 0
     *
     * @return The lambda function for resetting the force
     */

    static std::function<void(Particle&)> resetForceFunction() {
        return [](Particle &p) {
            p.updateF(std::array<double, 3> {0., 0., 0.});
        };
    };


    /**
     * @brief Returns a force vector, which if applied to a particle, would result in a force that is
     *        equal to the gravitational force working on the particle
     *
     * @param m The mass of the particle, which the force vector is calculated for
     * @param g The gravitational acceleration constant
     * @return The force vector representing the gravitational force working on the particle
     */
    static std::array<double, 3> verticalGravityForce(double m, double g) {
        return {0, m * g, 0};
    }

    /**
     * @brief Returns the function for force
     *
     * @return The function for force
     */
    std::function<void(Particle&, Particle&)> forceFunction() {
        return force;
    }

    /**
     * @brief Returns the function for position
     *
     * @return The lambda function for position
     */
    std::function<void(Particle&)> positionFunction() {
        return position;
    }

    /**
     * @brief Returns the function for velocity
     *
     * @return The lambda function for velocity
     */
    std::function<void(Particle&)> velocityFunction() {
        return velocity;
    }

    /**
     * @brief Return a basic model
     *
     * @param deltaT Time delta for the model
     * @return A Model object configured with the most basic formulas
     */
    static Model gravityModel(double deltaT) {
        auto force = [](Particle &p1, Particle &p2) {
            auto nextForce =
                    p1.getM() * p2.getM() / std::pow(p2.distanceTo(p1), 3) * p1.diffTo(p2);

            p1.setF(p1.getF() + nextForce);
        };

        auto position = [deltaT](Particle &p) {
            auto x =
                    p.getX() +
                    deltaT * p.getV() +
                    (deltaT * deltaT / (2 * p.getM())) * p.getF();

            p.setX(x);
        };

        auto velocity = [deltaT](Particle &p) {
            auto v =
                    p.getV() +
                    (deltaT / (2 * p.getM())) * (p.getOldF() + p.getF());

            p.setV(v);
        };

        return Model{force, position, velocity};
    };

    /**
     * @brief Return a model that utilizes Lennard Jones potential for force calculation between particles
     *
     * @param deltaT Time delta for the model
     * @param epsilon represents the depth of the potential well, which reflects the strength of the interaction between the particles
     * @param sigma distance at which the inter-particle potential is zero.
     * @return A Model object configured with the basic formulas for velocity, position calculation and Lennard Jones potential formula for force calculation
     */
    static Model lennardJonesModel(double deltaT) {
        auto ljForce = [](Particle &p1, Particle &p2) {
            // Lorentz-Berthelot mixing rules
            auto epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
            auto sigma = (p1.getSigma() + p2.getSigma()) / 2;

            auto distance = p2.distanceTo(p1);
            //todo remove later
            //distance = std::max(distance, 0.1);

            auto distance6 = std::pow(distance, 6);
            auto sigma6 = std::pow(sigma, 6);

            auto nextForce = (-24 * epsilon / std::pow(distance,2))
                             * ((sigma6/distance6) - 2*(std::pow(sigma6,2)/std::pow(distance6,2)))
                             * p2.diffTo(p1);

            if(!p1.isFixed()){
                 p1.setF(p1.getF() + nextForce);
            }
            if(!p2.isFixed()){
                 p2.setF(p2.getF() - nextForce);
            }


        };

        auto position = [deltaT](Particle &p) {
            auto x =
                    p.getX() +
                    deltaT * p.getV() +
                    (deltaT * deltaT / (2 * p.getM())) * p.getF();

            if(!p.isFixed()){
                p.setX(x);
            }
        };

        auto velocity = [deltaT](Particle &p) {
            auto v =
                    p.getV() +
                    (deltaT / (2 * p.getM())) * (p.getOldF() + p.getF());

            if(!p.isFixed()){
                p.setV(v);

            }
        };

        return Model{ljForce, position, velocity};
    }
};