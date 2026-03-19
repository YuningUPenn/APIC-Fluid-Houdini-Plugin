#pragma once
#include "apic_types.h"
#include "ApicGrid.h"
#include "ApicParticleSet.h"
#include <functional>

namespace apic {

// -------------------------------------------------------
// CollisionHandler
//
// Handles two types of boundaries:
//   1. Domain boundary (solid walls around the grid)
//   2. Static obstacle geometry represented as an SDF
//
// The SDF is sampled at grid nodes to mark solid cells.
// Velocity is zeroed (or reflected) on solid-fluid faces
// after the grid update and after G2P.
// -------------------------------------------------------
class CollisionHandler {
public:
    // SDF function type: returns signed distance at world position
    // Negative = inside solid
    using SDFFunc = std::function<Scalar(const Vec3&)>;

    CollisionHandler() = default;

    // Register a static collision SDF (can call multiple times
    // to union multiple obstacles -- min of SDFs)
    void addObstacle(SDFFunc sdf);

    // Evaluate combined SDF at pos
    Scalar evalSDF(const Vec3& pos) const;

    // Mark solid nodes on the grid based on registered SDFs
    // and the domain boundary.
    void markSolidNodes(ApicGrid& grid) const;

    // Enforce velocity boundary conditions on the grid:
    //   - Zero normal velocity on solid nodes
    //   - Apply friction to tangential velocity
    void enforceGridBoundary(ApicGrid& grid, const SimParams& params) const;

    // Resolve particle-solid collisions:
    //   - Push particles out of solid regions
    //   - Zero/reflect their velocity component into the solid
    void resolveParticleCollisions(ApicParticleSet& particles,
                                   const SimParams& params) const;

    // Convenience: build domain-boundary SDF for a box [0,1]^3 * scale
    static SDFFunc makeBoxSDF(const Vec3& minPt, const Vec3& maxPt);

private:
    std::vector<SDFFunc> obstacles_;
};

} // namespace apic
