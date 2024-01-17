/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"

#include <spdlog/spdlog.h>

Particle::Particle(int type_arg) {
    type = type_arg;
    spdlog::debug("Particle generated!");
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other)
        : x(other.x), v(other.v), f(other.f), old_f(other.old_f),
          m(other.m), epsilon(other.epsilon), sigma(other.sigma),
          type(other.type), directNeighbors(), diagonalNeighbors() {
    spdlog::debug("Particle generated by copy!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double eps, double sig, int type_arg, double avgBondLength_arg, int stiffnessFactor_arg, bool pulled_arg)
        : x(x_arg), v(v_arg), m(m_arg), type(type_arg),
          f({0., 0., 0.}), old_f({0., 0., 0.}), epsilon(eps), sigma(sig),
          directNeighbors(), diagonalNeighbors(), avgBondLength(avgBondLength_arg), stiffnessFactor(stiffnessFactor_arg),
          pulled(pulled_arg) {
    spdlog::debug("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   std::array<double, 3> f_arg, std::array<double, 3> old_f_arg,
                   double m_arg, double eps, double sig, int type_arg)
        : x(x_arg), v(v_arg), m(m_arg), f(f_arg), old_f(old_f_arg),
          type(type_arg), epsilon(eps), sigma(sig),
          directNeighbors(), diagonalNeighbors() {
    spdlog::debug("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   std::array<double, 3> f_arg, std::array<double, 3> old_f_arg,
                   double m_arg, double eps, double sig, int type_arg, double avgBondLength_arg, int stiffnessFactor_arg, bool pulled_arg)
        : x(x_arg), v(v_arg), m(m_arg), f(f_arg), old_f(old_f_arg),
          type(type_arg), epsilon(eps), sigma(sig),
          directNeighbors(), diagonalNeighbors(), avgBondLength(avgBondLength_arg), stiffnessFactor(stiffnessFactor_arg), pulled(pulled_arg) {
    spdlog::debug("Particle generated!");
}


Particle::~Particle() {
    spdlog::debug("Particle destructed!");
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

void Particle::updateF(const std::array<double, 3> &f_arg) {
    old_f = f;
    f = f_arg;
}

void Particle::setX(const std::array<double, 3> &x_arg) { x = x_arg; }

void Particle::setV(const std::array<double, 3> &v_arg) { v = v_arg; }

void Particle::setF(const std::array<double, 3> &f_arg) { f = f_arg; }

double Particle::getEpsilon() const { return epsilon; }

double Particle::getSigma() const { return sigma; }

void Particle::setEpsilon(double eps) {
    this->epsilon = eps;
}

void Particle::setSigma(double sig) {
    this->sigma = sig;
}

void Particle::setOldF(const std::array<double, 3> &old_f_arg) { old_f = old_f_arg; }

std::array<double, 3> Particle::diffTo(Particle &particle) {
    return particle.getX() - x;
}

double Particle::distanceTo(Particle &particle) {
    return ArrayUtils::L2Norm(diffTo(particle));
}

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

void Particle::setType(int type_arg) { type = type_arg; }

std::string Particle::toString() const {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f
           << " old_f: " << old_f << " type: " << type;
    return stream.str();
}

bool Particle::operator==(Particle &other) {
    return (x == other.x) and (v == other.v) and (f == other.f) and
           (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

nlohmann::ordered_json Particle::json() {
    nlohmann::ordered_json j;
    j["type"] = "particle";
    j["type_id"] = type;
    j["mass"] = m;
    j["position"] = x;
    j["velocity"] = v;
    j["force"] = f;
    j["old_force"] = old_f;
    j["epsilon"] = epsilon;
    j["sigma"] = sigma;

    return j;
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
    stream << p.toString();
    return stream;
}

void Particle::addDirectNeighbor(Particle* neighbor){
    directNeighbors.push_back(neighbor);
}

void Particle::addDiagonalNeighbor(Particle* neighbor){
    diagonalNeighbors.push_back(neighbor);
}

std::vector<Particle *> &Particle::getDirectNeighbors() {
    return directNeighbors;
}

std::vector<Particle *> &Particle::getDiagonalNeighbors() {
    return diagonalNeighbors;
}

double Particle::getAvgBondLength() {
    return avgBondLength;
}

int Particle::getStiffnessFactor() {
    return stiffnessFactor;
}

bool Particle::isPulled() {
    return pulled;
}