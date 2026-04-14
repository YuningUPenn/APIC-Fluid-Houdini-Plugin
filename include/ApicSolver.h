#pragma once
#include "apic_types.h"
#include "ApicParticleSet.h"
#include "ApicGrid.h"
#include "PressureSolver.h"
#include "CollisionHandler.h"

namespace apic {

// -------------------------------------------------------
// ApicSolver
//
// Central controller.  Owns or holds references to all
// simulation data and executes one full APIC timestep:
//
//   1. P2G  -- transfer mass and momentum from particles
//              to the grid (APIC uses affine B matrix)
//   2. Grid update -- apply external forces (gravity, etc.)
//   3. Pressure solve -- enforce incompressibility
//   4. G2P  -- interpolate updated grid velocities back to
//              particles and update B matrix
//   5. Advect -- integrate particle positions
//
// The method selector (APIC / FLIP / PIC / Hybrid) changes
// steps 1 and 4 only; everything else is identical.
// -------------------------------------------------------
class ApicSolver {
public:
    // ---- Construction ----
    explicit ApicSolver(const SimParams& params);

    // ---- Setup ----
    void setParams(const SimParams& params);
    void setCollisionHandler(CollisionHandler* handler) { collision_ = handler; }

    // Grid and particle access (owned by this solver)
    ApicGrid&         grid()      { return grid_; }
    ApicParticleSet&  particles() { return particles_; }
    const ApicGrid&   grid()      const { return grid_; }
    const ApicParticleSet& particles() const { return particles_; }
    const SimParams& params() const { return params_; }

    // ---- Main API ----
    // Advance simulation by one substep (dt = params.dt / params.substeps)
    void step();

    // Run a full frame (substeps iterations of step())
    void stepFrame();

    // Reset grid + particles
    void reset();

    // ---- Diagnostics ----
    Scalar lastPressureResidual() const { return pressureResidual_; }
    Scalar totalKineticEnergy() const;
    Vec3   totalAngularMomentum() const;
    

private:
    // ---- Simulation sub-steps ----

    // Step 1: Particle-to-Grid (P2G) transfer
    // Rasterises particle mass and momentum (+ affine contribution) to grid.
    void p2gTransfer();

    // Step 1a: APIC P2G -- includes affine B contribution
    void p2gAPIC();

    // Step 1b: PIC P2G -- pure mass/velocity transfer (no B)
    void p2gPIC();

    // Step 2: Grid update -- apply forces, handle boundary
    void gridUpdate(Scalar substepDt);

    // Step 3: Pressure projection
    void pressureSolve(Scalar substepDt);

    // Step 4: Grid-to-Particle (G2P) transfer
    // Interpolates grid velocities back and updates B.
    void g2pTransfer();

    // Step 4a: APIC G2P -- updates v_p and B_p
    void g2pAPIC();

    // Step 4b: PIC G2P -- updates v_p, zeros B_p
    void g2pPIC();

    // Step 4c: FLIP G2P -- v_p += delta_v, preserves old B
    void g2pFLIP();

    // Step 5: Particle advection
    void advectParticles(Scalar substepDt);

    // ---- Data ----
    SimParams       params_;
    ApicGrid        grid_;
    ApicParticleSet particles_;
    PressureSolver  pressureSolver_;
    CollisionHandler* collision_ = nullptr;  // not owned

    Scalar pressureResidual_ = 0.f;
};

} // namespace apic
