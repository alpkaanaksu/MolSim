//
// Created by berks on 17.11.2023.
//

#include "LinkedCellParticleContainer.h"

#include <iostream>
#include <spdlog/spdlog.h>
#include <cmath>

#include "Particle.h"
#include "../utils/ArrayUtils.h"


LinkedCellParticleContainer::LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT,
                                                         BoundaryBehavior boundaryBehaviorTop,
                                                         BoundaryBehavior boundaryBehaviorBottom,
                                                         BoundaryBehavior boundaryBehaviorRight,
                                                         BoundaryBehavior boundaryBehaviorLeft,
                                                         BoundaryBehavior boundaryBehaviorFront,
                                                         BoundaryBehavior boundaryBehaviorBack)
        :
        xSize(xSize), ySize(ySize), zSize(zSize),
        xCells(static_cast<int>(std::floor(xSize / cutoffRadius))),
        yCells(static_cast<int>(std::floor(ySize / cutoffRadius))),
        zCells(static_cast<int>(std::floor(zSize / cutoffRadius))),
        cutoffRadius(cutoffRadius),
        deltaT(deltaT), boundaryBehaviorTop(boundaryBehaviorTop),
        boundaryBehaviorBottom(boundaryBehaviorBottom), boundaryBehaviorRight(boundaryBehaviorRight),
        boundaryBehaviorLeft(boundaryBehaviorLeft), boundaryBehaviorFront(boundaryBehaviorFront),
        boundaryBehaviorBack(boundaryBehaviorBack) {

    cellXSize = xSize / xCells;
    cellYSize = ySize / yCells;
    cellZSize = zSize / zCells;

    xCells += 2;
    yCells += 2;
    zCells += 2;

    int numberOfCells = xCells * yCells * zCells;

    cells = std::vector<std::vector<Particle>>(numberOfCells);

    isHaloCellVector = std::vector<bool>(numberOfCells, true);

    boundaryCellIndices = std::set<int>();
    haloCellIndices = std::set<int>();

    // Precompute cells-vector indices of the halo cells in and store them in haloCellIndices
    for (int x = 0; x < xCells; x++) {
        for (int y = 0; y < yCells; y++) {
            haloCellIndices.insert(index3dTo1d(x, y, 0));                 // Bottom halo
            haloCellIndices.insert(index3dTo1d(x, y, zCells - 1));        // Top halo
        }
    }

    for (int x = 0; x < xCells; x++) {
        for (int z = 0; z < zCells; z++) {
            haloCellIndices.insert(index3dTo1d(x, 0, z));                 // Front halo
            haloCellIndices.insert(index3dTo1d(x, yCells - 1, z));        // Back halo
        }
    }

    for (int y = 0; y < yCells; y++) {
        for (int z = 0; z < zCells; z++) {
            haloCellIndices.insert(index3dTo1d(0, y, z));                 // Left halo
            haloCellIndices.insert(index3dTo1d(xCells - 1, y, z));        // Right halo
        }
    }


    // Precompute cells-vector indices of the boundary cells in and store them in boundaryCellIndices
    for (int x = 1; x < xCells - 1; x++) {
        for (int y = 1; y < yCells - 1; y++) {
            boundaryCellIndices.insert(index3dTo1d(x, y, 1));                 // Front boundary
            boundaryCellIndices.insert(index3dTo1d(x, y, zCells - 2));        // Back boundary
        }
    }

    for (int x = 1; x < xCells - 1; x++) {
        for (int z = 1; z < zCells - 1; z++) {
            boundaryCellIndices.insert(index3dTo1d(x, 1, z));                 // Bottom boundary
            boundaryCellIndices.insert(index3dTo1d(x, yCells - 2, z));        // Top boundary
        }
    }

    for (int y = 1; y < yCells - 1; y++) {
        for (int z = 1; z < zCells - 1; z++) {
            boundaryCellIndices.insert(index3dTo1d(1, y, z));                 // Left boundary
            boundaryCellIndices.insert(index3dTo1d(xCells - 2, y, z));        // Right boundary
        }
    }

    // Iterate through haloCellIndices and set corresponding pairs to false
    for (int haloIndex : haloCellIndices) {
        isHaloCellVector[haloIndex] = false;
    }

}

LinkedCellParticleContainer::LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT)
        : LinkedCellParticleContainer(xSize, ySize, zSize, cutoffRadius, deltaT, BoundaryBehavior::Reflective) {}

LinkedCellParticleContainer::LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT, BoundaryBehavior boundaryBehavior)
        : LinkedCellParticleContainer(xSize, ySize, zSize, cutoffRadius, deltaT, boundaryBehavior, boundaryBehavior, boundaryBehavior, boundaryBehavior, boundaryBehavior, boundaryBehavior) {}

LinkedCellParticleContainer::~LinkedCellParticleContainer() = default;

int LinkedCellParticleContainer::index3dTo1d(int x, int y, int z) {
    return (x + y * xCells + z * xCells * yCells);
}

std::array<int, 3> LinkedCellParticleContainer::index1dTo3d(int index) {
    int x = index % xCells;
    int y = (index / xCells) % yCells;
    int z = index / (xCells * yCells);

    return {x, y, z};
}

int LinkedCellParticleContainer::cellIndexForParticle(const Particle &particle) {
    int xIndex = static_cast<int>((std::floor((particle.getX()[0]) / cellXSize)));
    int yIndex = static_cast<int>((std::floor((particle.getX()[1]) / cellYSize)));
    int zIndex = static_cast<int>((std::floor((particle.getX()[2]) / cellZSize)));

    if (xIndex < 0 || xIndex >= xCells || yIndex < 0 || yIndex >= yCells || zIndex < 0 || zIndex >= zCells) {
        spdlog::info("Particle out of bounds: {}, {}, {}", particle.getX()[0], particle.getX()[1], particle.getX()[2]);
        // -1 means halo cell
        return -1;
    }

    return (xIndex+1) + (yIndex+1) * xCells + (zIndex+1) * xCells * yCells;
}


void LinkedCellParticleContainer::applyToAllPairsOnce(const std::function<void(Particle&, Particle&)>& function) {
    // Iterate through all cells in the container
    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        // Skip halo cells
        if (!isHaloCellVector[cellIndex]) continue;

        auto coords = index1dTo3d(cellIndex);
        auto &firstCell = cells[cellIndex];  // Extract the vector of particles from the pair

        // Iterate through all pairs of particles in the same cell
        for (int i = 0; i < firstCell.size(); i++) {
            for (int j = i + 1; j < firstCell.size(); j++) {
                // Check if the pair has been processed before by comparing memory addresses
                if (&firstCell[i] < &firstCell[j]) {
                    function(firstCell[i], firstCell[j]);
                }
            }
        }

        // Iterate through neighboring cells
        for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <= 1; y++) {
                for (int z = -1; z <= 1; z++) {
                    int neighborX = coords[0] + x;
                    int neighborY = coords[1] + y;
                    int neighborZ = coords[2] + z;

                    /*if (neighborX <= 0 || neighborX >= xCells - 1
                        || neighborY <= 0 || neighborY >= yCells - 1
                        || neighborZ <= 0 || neighborZ >= zCells - 1) continue;*/

                    if (x == 0 && y == 0 && z == 0) continue;

                    auto &currentCell = cells[index3dTo1d(neighborX, neighborY,
                                                          neighborZ)];  // Extract vector from the pair

                    for (auto &p1: firstCell) {
                        for (auto &p2: currentCell) {
                            // Check if the pair has been processed before by comparing memory addresses
                            if (&p1 < &p2) {
                                function(p1, p2);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*void LinkedCellParticleContainer::counterParticleReflection(Particle &particle, int cellIndex, const std::function<void(Particle&, Particle&)>& function) {
    std::array<int, 3> cellCoords = index1dTo3d(cellIndex);
    std::array<double, 3> particleCoords = particle.getX();

    if (cellCoords[0] == 1) {
        particleCoords[0] = -particleCoords[0];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    } else if (cellCoords[0] == xCells - 2) {
        particleCoords[0] = 2 * xSize - particleCoords[0];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    }

    if(cellCoords[1] == 1) {
        particleCoords[1] = -particleCoords[1];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    } else if (cellCoords[1] == yCells - 2) {
        particleCoords[1] = 2 * ySize - particleCoords[1];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    }

    if(cellCoords[2] == 1) {
        particleCoords[2] = -particleCoords[2];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    } else if (cellCoords[2] == zCells - 2) {
        particleCoords[2] = 2 * zSize - particleCoords[2];
        Particle reflectedParticle(particleCoords, particle.getV(), particle.getM(), particle.getEpsilon(),particle.getSigma(), particle.getType());
        function(particle, reflectedParticle);
    }
}*/

void LinkedCellParticleContainer::applyToAll(const std::function<void(Particle&)>& function) {
    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        if (!isHaloCellVector[cellIndex]) continue;  // Skip processing for halo cells

        for (auto& particle : cells[cellIndex]) {
            function(particle);
        }
    }
}

void LinkedCellParticleContainer::applyToAllHalo(const std::function<void(Particle&)>& function) {
    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        for (auto& particle : cells[cellIndex]) {
            particle.setType(isHaloCellVector[cellIndex] ? 2 : 1);

            function(particle);
        }
    }
}

void LinkedCellParticleContainer::applyToAll(const std::function<void(Particle&)>& function, bool updateCells) {
    deleteParticlesInHaloCells();

    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        if (!isHaloCellVector[cellIndex]) continue;  // Skip processing for halo cells

        for (auto& particle : cells[cellIndex]) {
            function(particle);
        }

        if (updateCells) {
            updateParticleCell(cellIndex);
        }
    }

    updateHaloCells();
}

void LinkedCellParticleContainer::add(const Particle &particle) {
    int cellIndex = cellIndexForParticle(particle);

    if (cellIndex != -1) {
        addParticleToCell(cellIndex, particle);
    } else {
        spdlog::info("Particle out of bounds: {}, {}, {}", particle.getX()[0], particle.getX()[1], particle.getX()[2]);
    }

}

void LinkedCellParticleContainer::addParticleToCell(int cellIndex, const Particle &particle) {
    cells[cellIndex].push_back(particle);
}

void LinkedCellParticleContainer::updateParticleCell(int cellIndex) {
    auto &cell = cells[cellIndex];

    for (auto it = cell.begin(); it != cell.end();) {

        vectorReverseReflection(*it);
        handlePeriodicBoundary(*it);

        int newCellIndex = cellIndexForParticle(*it);

        if (newCellIndex != cellIndex) {
            add(*it);
            it = cell.erase(it);  // Remove the particle from the old cell
        } else {
            ++it;
        }
    }
}

int LinkedCellParticleContainer::size() {
    int size = 0;

    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        if (!isHaloCellVector[cellIndex]) continue;  // Skip processing for halo cells

        size += static_cast<int>(cells[cellIndex].size());
    }

    return size;
}

double LinkedCellParticleContainer::updatePositionOnUpperPeriodic(const double axisPosition, int axisIndex) {
    double maxSize = axisIndex == 0 ? xSize : axisIndex == 1 ? ySize : zSize;

    return axisPosition - maxSize;
}

double LinkedCellParticleContainer::updatePositionOnLowerPeriodic(const double axisPosition, int axisIndex) {
    double maxSize = axisIndex == 0 ? xSize : axisIndex == 1 ? ySize : zSize;

    return axisPosition + maxSize;
}


void LinkedCellParticleContainer::handlePeriodicBoundary(Particle& particle) {
    std::array<double, 3> updatedPosition = particle.getX();

    if(particle.getX()[0] > xSize && boundaryBehaviorRight == BoundaryBehavior::Periodic) {
        updatedPosition[0] = updatePositionOnUpperPeriodic(updatedPosition[0], 0);
    }else if(particle.getX()[0] < 0 && boundaryBehaviorLeft == BoundaryBehavior::Periodic) {
        updatedPosition[0] = updatePositionOnLowerPeriodic(updatedPosition[0], 0);
    }

    if(particle.getX()[1] > ySize && boundaryBehaviorTop == BoundaryBehavior::Periodic) {
        updatedPosition[1] = updatePositionOnUpperPeriodic(updatedPosition[1], 1);
    } else if(particle.getX()[1] < 0 && boundaryBehaviorBottom == BoundaryBehavior::Periodic) {
        updatedPosition[1] = updatePositionOnLowerPeriodic(updatedPosition[1], 1);
    }

    if(particle.getX()[2] > zSize && boundaryBehaviorFront == BoundaryBehavior::Periodic) {
        updatedPosition[2] = updatePositionOnUpperPeriodic(updatedPosition[2], 2);
    } else if(particle.getX()[2] < 0 && boundaryBehaviorBack == BoundaryBehavior::Periodic) {
        updatedPosition[2] = updatePositionOnLowerPeriodic(updatedPosition[2], 2);
    }

    particle.setX(updatedPosition);
}

void LinkedCellParticleContainer::deleteParticlesInHaloCells() {
    // Iterate through halo cells
    for (int haloIndex : haloCellIndices) {
        // Delete particles in the halo cell
        cells[haloIndex].clear();
    }
}

void LinkedCellParticleContainer::updateHaloCells() {
    for (int boundaryCellIndex : boundaryCellIndices) {
        std::array<int, 3> boundary3d = index1dTo3d(boundaryCellIndex);

        if (boundaryBehaviorLeft == BoundaryBehavior::Periodic && boundary3d[0] == xCells - 2) {
            upperBoundaryToLowerHaloOneAxis(boundaryCellIndex, 0);
        } else if (boundaryBehaviorRight == BoundaryBehavior::Periodic && boundary3d[0] == 1) {
            lowerBoundaryToUpperHaloOneAxis(boundaryCellIndex, 0);
        }

        if (boundaryBehaviorBottom == BoundaryBehavior::Periodic && boundary3d[1] == yCells - 2) {
            upperBoundaryToLowerHaloOneAxis(boundaryCellIndex, 1);
        } else if (boundaryBehaviorTop == BoundaryBehavior::Periodic && boundary3d[1] == 1) {
            lowerBoundaryToUpperHaloOneAxis(boundaryCellIndex, 1);
        }

        if (boundaryBehaviorBack == BoundaryBehavior::Periodic && boundary3d[2] == zCells - 2) {
            upperBoundaryToLowerHaloOneAxis(boundaryCellIndex, 2);
        } else if (boundaryBehaviorFront == BoundaryBehavior::Periodic && boundary3d[2] == 1) {
            lowerBoundaryToUpperHaloOneAxis(boundaryCellIndex, 2);
        }
    }
}

void LinkedCellParticleContainer::upperBoundaryToLowerHaloOneAxis(int boundaryCellIndex, int axisIndex) {
    double maxSize = axisIndex == 0 ? xSize : axisIndex == 1 ? ySize : zSize;
    int minCells = 0;

    std::array<int, 3> halo3d = index1dTo3d(boundaryCellIndex);
    halo3d[axisIndex] = minCells;
    int haloCellIndex = index3dTo1d(halo3d[0], halo3d[1], halo3d[2]);

    for(auto &particle : cells[boundaryCellIndex]) {
        std::array<double, 3> updatedPosition = particle.getX();
        updatedPosition[axisIndex] = updatedPosition[axisIndex] - maxSize;
        Particle newParticle = Particle(updatedPosition, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());

        //spdlog::info("Adding particle to halo cell: {}, {}, {}", newParticle.getX()[0], newParticle.getX()[1], newParticle.getX()[2]);
        addParticleToCell(haloCellIndex, newParticle);
    }
}


void LinkedCellParticleContainer::lowerBoundaryToUpperHaloOneAxis(int boundaryCellIndex, int axisIndex) {
    double maxSize = axisIndex == 0 ? xSize : axisIndex == 1 ? ySize : zSize;
    int maxCells = axisIndex == 0 ? xCells : axisIndex == 1 ? yCells : zCells;

    std::array<int, 3> halo3d = index1dTo3d(boundaryCellIndex);
    halo3d[axisIndex] = maxCells - 1;
    int haloCellIndex = index3dTo1d(halo3d[0], halo3d[1], halo3d[2]);

    for(auto &particle : cells[boundaryCellIndex]) {
        std::array<double, 3> updatedPosition = particle.getX();
        updatedPosition[axisIndex] = updatedPosition[axisIndex] + maxSize;
        Particle newParticle = Particle(updatedPosition, particle.getV(), particle.getM(), particle.getEpsilon(), particle.getSigma(), particle.getType());

        addParticleToCell(haloCellIndex, newParticle);
    }
}


//Final alternative is this:
void LinkedCellParticleContainer::reflectIfNecessaryOnAxis(Particle& particle, double axisMin, double axisMax, int axisIndex) {
    std::array<double, 3> position = particle.getX();
    std::array<double, 3> velocity = particle.getV();

    while (true) {
        bool reflected = false;

        if (axisIndex == 0) {
            if (boundaryBehaviorRight == BoundaryBehavior::Reflective && position[0] > axisMax) {
                position[0] = 2 * axisMax - position[0];
                velocity[0] = -velocity[0];
                reflected = true;
            } else if (boundaryBehaviorLeft == BoundaryBehavior::Reflective && position[0] < axisMin) {
                position[0] = 2 * axisMin - position[0];
                velocity[0] = -velocity[0];
                reflected = true;
            }
        } else if (axisIndex == 1) {
            if (boundaryBehaviorTop == BoundaryBehavior::Reflective && position[1] > axisMax) {
                position[1] = 2 * axisMax - position[1];
                velocity[1] = -velocity[1];
                reflected = true;
            } else if (boundaryBehaviorBottom == BoundaryBehavior::Reflective && position[1] < axisMin) {
                position[1] = 2 * axisMin - position[1];
                velocity[1] = -velocity[1];
                reflected = true;
            }
        } else if (axisIndex == 2) {
            if (boundaryBehaviorFront == BoundaryBehavior::Reflective && position[2] > axisMax) {
                position[2] = 2 * axisMax - position[2];
                velocity[2] = -velocity[2];
                reflected = true;
            } else if (boundaryBehaviorBack == BoundaryBehavior::Reflective && position[2] < axisMin) {
                position[2] = 2 * axisMin - position[2];
                velocity[2] = -velocity[2];
                reflected = true;
            }
        }

        if (!reflected) {
            particle.setX(position);
            particle.setV(velocity);
            break;  // Break the loop if no further reflection is needed
        }
    }
}


void LinkedCellParticleContainer::vectorReverseReflection(Particle& particle) {
    reflectIfNecessaryOnAxis(particle,  0, xSize, 0);
    reflectIfNecessaryOnAxis(particle, 0, ySize, 1);
    reflectIfNecessaryOnAxis(particle,  0, zSize, 2);
}

nlohmann::ordered_json LinkedCellParticleContainer::json() {
    nlohmann::ordered_json j;

    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
        if (!isHaloCellVector[cellIndex]) continue;  // Skip processing for halo cells

        for (auto& particle : cells[cellIndex]) {
            j.push_back(particle.json());
        }
    }

    return j;
}

std::string LinkedCellParticleContainer::toString() {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2)
           << "\n--- Container ---------"
           << "\nType: linked-cell"
           << "\nDomain size: " << xSize << " • " << ySize << " • " << zSize
           << "\nCell size: " << cellXSize << " • " << cellYSize << " • " << cellZSize
           << "\nNumber of particles: " << size()
           << "\nBoundary behavior (top): " << boundaryBehaviorToString(boundaryBehaviorTop)
           << "\nBoundary behavior (bottom): " << boundaryBehaviorToString(boundaryBehaviorBottom)
           << "\nBoundary behavior (left): " << boundaryBehaviorToString(boundaryBehaviorLeft)
           << "\nBoundary behavior (right): " << boundaryBehaviorToString(boundaryBehaviorRight)
           << "\nBoundary behavior (front): " << boundaryBehaviorToString(boundaryBehaviorFront)
           << "\nBoundary behavior (back): " << boundaryBehaviorToString(boundaryBehaviorBack)
           << "\n------------------------\n";

    return stream.str();
}

std::ostream &operator<<(std::ostream &stream, LinkedCellParticleContainer &simulation) {
    stream << simulation.toString();
    return stream;
}

double LinkedCellParticleContainer::getXSize() const {
    return xSize;
}

double LinkedCellParticleContainer::getYSize() const {
    return ySize;
}

double LinkedCellParticleContainer::getZSize() const {
    return zSize;
}

int LinkedCellParticleContainer::getXCells() const {
    return xCells;
}

int LinkedCellParticleContainer::getYCells() const {
    return yCells;
}

int LinkedCellParticleContainer::getZCells() const {
    return zCells;
}

double LinkedCellParticleContainer::getCutoffRadius() const {
    return cutoffRadius;
}

double LinkedCellParticleContainer::getCellXSize() const {
    return cellXSize;
}

double LinkedCellParticleContainer::getCellYSize() const {
    return cellYSize;
}

double LinkedCellParticleContainer::getCellZSize() const {
    return cellZSize;
}

double LinkedCellParticleContainer::getDeltaT() const {
    return deltaT;
}

const std::vector<std::vector<Particle>> &LinkedCellParticleContainer::getCells() const {
    return cells;
}

const std::set<int> &LinkedCellParticleContainer::getBoundaryCellIndices() const {
    return boundaryCellIndices;
}

const std::set<int> &LinkedCellParticleContainer::getHaloCellIndices() const {
    return haloCellIndices;
}

const std::vector<bool> &LinkedCellParticleContainer::getIsHaloCellVector() const {
    return isHaloCellVector;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorTop() const {
    return boundaryBehaviorTop;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorBottom() const {
    return boundaryBehaviorBottom;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorRight() const {
    return boundaryBehaviorRight;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorLeft() const {
    return boundaryBehaviorLeft;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorFront() const {
    return boundaryBehaviorFront;
}

BoundaryBehavior LinkedCellParticleContainer::getBoundaryBehaviorBack() const {
    return boundaryBehaviorBack;
}
