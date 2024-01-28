//
// Created by berks on 17.11.2023.
//

#pragma once

#include <vector>
#include <deque>
#include <array>
#include <list>
#include <set>
#include "Particle.h"
#include "ParticleContainer.h"
#include "BoundaryBehavior.h"

/**
 * A container that stores particles in a linked cell data structure. The cells are stored in a 1D vector.
 *
 * This implementation is significantly faster than the naive implementation
 * @image html plot.png width=500px
 */
class LinkedCellParticleContainer : public ParticleContainer{
private:
    double xSize;
    double ySize;
    double zSize;

    int xCells;
    int yCells;
    int zCells;

    double  cutoffRadius;

    double cellXSize;
    double cellYSize;
    double cellZSize;

    double deltaT;
    /**
     * The vector that contains all the particles in the container
     */

    std::vector<std::vector<Particle>> cells;

    std::set<int> boundaryCellIndices;

    std::set<int> haloCellIndices;

    std::vector<bool> isHaloCellVector;

    BoundaryBehavior boundaryBehaviorTop;
    BoundaryBehavior boundaryBehaviorBottom;
    BoundaryBehavior boundaryBehaviorRight;
    BoundaryBehavior boundaryBehaviorLeft;
    BoundaryBehavior boundaryBehaviorFront;
    BoundaryBehavior boundaryBehaviorBack;

public:

    /**
     * @brief Construct a new Linked Cell Particle Container object
     * @note Delegates to the other constructor with all boundary behaviors set to 'Reflective'
     * @param xSize
     * @param ySize
     * @param zSize
     * @param cutoffRadius
     * @param deltaT
     */
    LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT);

    /**
     * @brief Construct a new Linked Cell Particle Container object
     * @param xSize
     * @param ySize
     * @param zSize
     * @param cutoffRadius
     * @param deltaT
     * @param boundaryBehaviorTop
     * @param boundaryBehaviorBottom
     * @param boundaryBehaviorRight
     * @param boundaryBehaviorLeft
     * @param boundaryBehaviorFront
     * @param boundaryBehaviorBack
     */
    LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT, BoundaryBehavior boundaryBehaviorTop, BoundaryBehavior boundaryBehaviorBottom, BoundaryBehavior boundaryBehaviorRight, BoundaryBehavior boundaryBehaviorLeft, BoundaryBehavior boundaryBehaviorFront, BoundaryBehavior boundaryBehaviorBack);

    /**
     * @brief Construct a new Linked Cell Particle Container object
     * @note Delegates to the other constructor with all boundary behaviors set to the given value
     * @param xSize
     * @param ySize
     * @param zSize
     * @param cutoffRadius
     * @param deltaT
     * @param boundaryBehaviorAll
     */
    LinkedCellParticleContainer(double xSize, double ySize, double zSize, double cutoffRadius, double deltaT, BoundaryBehavior boundaryBehaviorAll);

    ~LinkedCellParticleContainer();

    /**
     * @brief Get the index of a cell with given coordinates (3d index)
     *
     * @param x
     * @param y
     * @param z
     * @return 1d index
     */
    int index3dTo1d(int x, int y, int z);

    /**
     * @brief Get the coordinates (3d index) of a cell with given index (1d)
     *
     * @param index
     * @return Coordinates
     */
    std::array<int, 3> index1dTo3d(int index);

    /**
     * @brief Apply a lambda function to all particle pairs
     * Checks the neighbor cells
     * Each combination is iterated only once (e.g. useful for Newton's 3rd law)
     * @param function
     */
    virtual void applyToAllPairsOnce(const std::function<void(Particle &, Particle &)> &function);

    /**
     * @brief Apply a lambda function to all particles
     * @param function
     */
    virtual void applyToAll(const std::function<void(Particle &)> &function);

    /**
     * @brief Find the correct cell for a given particle
     * @param particle
     *
     * @return 1d index of the cell the given particle belongs to
     */
    int cellIndexForParticle(const Particle &particle);

    /**
     * @brief Add a particle to the container
     * @note The particle is added to the correct cell (using cellIndexForParticle)
     * @param particle
     */
    void add(const Particle &particle);

    /**
     * @brief Update particles in the cell adressed by the given index
     * Remove the particles that are not in the cell anymore, add them to new cells if necessary
     * @param cellIndex
     *
     * @param cellIndex
     */
    void updateParticleCell(int cellIndex);

    /**
     * @brief Apply a lambda function to all particles
     * @note updateParticleCell is called if updateCells is true
     *
     * @param function
     * @param updateCells
     */
    void applyToAll(const std::function<void(Particle &)> &function, bool updateCells);

    /**
     * @brief Add a particle to the cell at the given index
     *
     * @param cellIndex
     * @param particle
     */
    void addParticleToCell(int cellIndex, const Particle &particle);

    /**
     * @brief Return number of particles in the container
     *
     * @return Number of particles
     */
    int size();

    nlohmann::ordered_json json();

    std::string toString();

    double getXSize() const;

    double getYSize() const;

    double getZSize() const;

    int getXCells() const;

    int getYCells() const;

    int getZCells() const;

    double getCutoffRadius() const;

    double getCellXSize() const;

    double getCellYSize() const;

    double getCellZSize() const;

    double getDeltaT() const;

    const std::vector<std::vector<Particle>> &getCells() const;

    const std::set<int> &getBoundaryCellIndices() const;

    const std::set<int> &getHaloCellIndices() const;

    const std::vector<bool> &getIsHaloCellVector() const;

    BoundaryBehavior getBoundaryBehaviorTop() const;

    BoundaryBehavior getBoundaryBehaviorBottom() const;

    BoundaryBehavior getBoundaryBehaviorRight() const;

    BoundaryBehavior getBoundaryBehaviorLeft() const;

    BoundaryBehavior getBoundaryBehaviorFront() const;

    BoundaryBehavior getBoundaryBehaviorBack() const;

    /**
    * @brief Makes the particle that crosses upper boundary to reappear on the lower boundary
    *
    * Given the current position of the particle and the corresponding axis index, the function calculates
    * the maxSize of the axis for which the particle has crossed the boundary. Then, it returns the updated
    * position of the particle by subtracting the maxSize from the current position, making the particle reappear
    * on the lower boundary (minSize, usually 0.0).
    *
    * @param axisPosition The current position along the specified axis.
    * @param axisIndex The index of the axis for which the position is updated.
    * @return The updated position on the upper periodic boundary.
    */
    double updatePositionOnUpperPeriodic(const double axisPosition, int axisIndex);

    /**
     * @brief Makes the particle that crosses lower boundary to reappear on the upper boundary
     *
     * Given the current position of the particle and the corresponding axis index, the function calculates
     * the maxSize of the axis for which the particle has crossed the boundary. Then, it returns the updated
     * position of the particle by adding the maxSize to the current position, making the particle reappear
     * on the upper boundary (maxSize).
     *
     * @param axisPosition The current position along the specified axis.
     * @param axisIndex The index of the axis for which the position is updated.
     * @return The updated position on the lower periodic boundary.
     */
    double updatePositionOnLowerPeriodic(const double axisPosition, int axisIndex);



    /**
     * @brief Makes the particle that crosses a boundary reappear on the opposite boundary
     *
     * This method checks if a particle has crossed the simulation boundary for all axis and if in the simulation
     * periodic boundary conditions are applied for that axis. If so, the particle is moved to the opposite
     * boundary alongside the axis with the help of functions updatePositionOnUpperPeriodic and
     * updatePositionOnLowerPeriodic.
     *
     * @param particle The particle for which periodic boundary conditions are applied.
     * @return
     */
    void handlePeriodicBoundary(Particle &particle);

    /**
     * @brief Handles periodic boundary conditions for a single axis.
     *
     * Copies particles from a boundary cell to the corresponding halo cell along the specified axis.
     *
     * @param boundaryCellIndex Index of the boundary cell to be handled.
     * @param axisIndex Axis along which the boundary condition is applied (0 for x, 1 for y, 2 for z).
     * @param isLowerHalo Indicates whether it is the lower halo cell (true) or upper halo cell (false).
     */
    void handleBoundariesOneAxis(int boundaryCellIndex, int axisIndex, bool isLowerHalo);

    /**
     * @brief Handles periodic boundary conditions for two axes simultaneously.
     *
     * Copies particles from a boundary cell to the corresponding halo cell along the specified axes. Useful when 2 of the
     * 3 axes have periodic boundary conditions.
     *
     * @param boundaryCellIndex Index of the boundary cell to be handled.
     * @param axisIndex1 First axis along which the boundary condition is applied.
     * @param axisIndex2 Second axis along which the boundary condition is applied.
     * @param isLowerHalo1 Indicates whether it is the lower halo cell for the first axis (true) or upper halo cell (false).
     * @param isLowerHalo2 Indicates whether it is the lower halo cell for the second axis (true) or upper halo cell (false).
     */
    void handleBoundariesTwoAxes(int boundaryCellIndex, int axisIndex1, int axisIndex2, bool isLowerHalo1, bool isLowerHalo2);

    /**
     * @brief Handles periodic boundary conditions for three axes simultaneously.
     *
     * Copies particles from a boundary cell to the corresponding halo cell along the specified axes. Useful when all 3
     * axes have periodic boundary conditions.
     *
     * @param boundaryCellIndex Index of the boundary cell to be handled.
     * @param isLowerHalo1 Indicates whether it is the lower halo cell for the first axis (true) or upper halo cell (false).
     * @param isLowerHalo2 Indicates whether it is the lower halo cell for the second axis (true) or upper halo cell (false).
     * @param isLowerHalo3 Indicates whether it is the lower halo cell for the third axis (true) or upper halo cell (false).
     */

    void handleBoundariesThreeAxes(int boundaryCellIndex, bool isLowerHalo1, bool isLowerHalo2, bool isLowerHalo3);


    /**
     * @brief Updates halo cells in the linked cell particle container.
     *
     * This method iterates through the boundary cells and, based on the specified boundary behaviors,
     * copies particles to the corresponding halo cells in a periodic manner, ensuring the periodicity
     * of the simulation for maintaining accurate boundary conditions.
     */
    void updateHaloCells();




    /**
     * @brief Apply a lambda function to all particles including the particles in the halo layer (cells)
     *
     * This method iterates through all cells in the particle container and applies the given lambda function
     * to all particles in each cell, including the particles in the halo layer. This method is useful for
     * testing purposes, as we used it with the plotParticles method inside the VTKWriter to see on Paraview if
     * our algorithm for updating halo cells for periodic boundaries work correctly.
     *
     * @param function
     */
    void applyToAllHalo(const std::function<void(Particle &)> &function);



    /**
     * @brief Remove all particles from halo cells.
     * Clear halo cells.
     */
    void deleteParticlesInHaloCells();



    /**
     * @brief Check and reflect a particle
     * @note Delegates the necessary reflections along all axis to reflectIfNecessaryOnAxis method.
     *
     * @param particle The particle to be checked and reflected.
     */
    void vectorReverseReflection(Particle &particle);



    /**
    * @brief Reflects the particle if it has crossed the specified axis boundaries.
    *
    * This method checks if a particle has crossed the boundaries along a specific axis and reflects
    * the particle's position and velocity accordingly. It considers the specified axis boundaries
    * (axisMin and axisMax) and the reflective behavior defined for each boundary condition. In essence
    * the function checks independently along each axis, if the new position of a particle is outside the
    * simulation boundaries. If it is, the function calculates how far the particle has traveled beyond the
    * boundary in the last iteration. It then allows the particle to move up to the boundary and reflects the
    * remaining distance in the opposite direction. Additionally, the function reverses the dimension component
    * in the particle's velocity that corresponds to the boundary it reached.
    *
    * @param particle The particle to be reflected if necessary.
    * @param axisMin The minimum boundary value along the specified axis.
    * @param axisMax The maximum boundary value along the specified axis.
    * @param axisIndex The index of the axis along which the reflection is applied (0 for x, 1 for y, 2 for z).
    *
    */
    void reflectIfNecessaryOnAxis(Particle &particle, double axisMin, double axisMax, int axisIndex);

    void applyToAllPairsOnceMembrane(const std::function<void(Particle &, Particle &)> &function);

};


