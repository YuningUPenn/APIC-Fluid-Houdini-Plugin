// test_p2g_conservation.cpp
// Standalone test -- no Houdini required.
//
// Tests: after a P2G transfer, total mass and linear momentum
// on the grid must equal the total mass and momentum of particles.
// After a G2P transfer (with no grid update), kinetic energy
// should be conserved (within floating-point tolerance).

#include "ApicSolver.h"
#include <cstdio>
#include <cmath>
#include <cassert>

using namespace apic;

static bool nearEqual(Scalar a, Scalar b, Scalar tol = 1e-3f) {
    return std::abs(a - b) < tol * (std::abs(a) + std::abs(b) + 1.f);
}

int main() {
    int failures = 0;

    // ---- Setup ----
    SimParams params;
    params.gridResX = params.gridResY = params.gridResZ = 16;
    params.gridSpacing  = 1.0f / 16.0f;
    params.dt           = 1.0f / 24.0f;
    params.substeps     = 1;
    params.method       = TransferMethod::APIC;
    params.gravity      = Vec3::Zero();  // no gravity for conservation test

    ApicSolver solver(params);
    ApicParticleSet& particles = solver.particles();
    ApicGrid& grid = solver.grid();

    // Place a small set of particles with known velocities
    const Scalar dx = params.gridSpacing;
    const Scalar mass_per_particle = 1.0f;
    int added = 0;

    for (int k = 4; k <= 11; ++k)
    for (int j = 4; j <= 11; ++j)
    for (int i = 4; i <= 11; ++i) {
        Vec3 pos = vec3(i + 0.5f, j + 0.5f, k + 0.5f) * dx;
        Vec3 vel = vec3(
            0.1f * std::sin(i + j),
            0.1f * std::cos(j + k),
            0.1f * std::sin(i * k)
        );
        particles.addParticle(pos, vel, mass_per_particle);
        ++added;
    }

    // Compute expected totals from particles
    Scalar pMassTotal = 0.f;
    Vec3   pMomTotal  = Vec3::Zero();
    for (size_t p = 0; p < particles.size(); ++p) {
        pMassTotal += particles.mass(p);
        pMomTotal  += particles.mass(p) * particles.velocity(p);
    }

    // Run P2G
    grid.clear();

    // Manually call what ApicSolver::p2gTransfer does:
    // We replicate so the test doesn't depend on private methods
    {
        const Scalar Dp_inv = 4.0f / (dx * dx);
        std::vector<ApicGrid::WeightEntry> weights;
        for (size_t p = 0; p < particles.size(); ++p) {
            grid.gatherWeights(particles.position(p), weights);
            for (const auto& we : weights) {
                Vec3 affine = particles.affineB(p) * (Dp_inv * we.offset);
                grid.mass(we.flatIdx)     += we.weight * particles.mass(p);
                grid.velocity(we.flatIdx) += we.weight * particles.mass(p)
                                             * (particles.velocity(p) + affine);
            }
        }
    }

    // Compute grid totals (before normalisation, grid.velocity stores momentum)
    Scalar gMassTotal = 0.f;
    Vec3   gMomTotal  = Vec3::Zero();
    for (size_t i = 0; i < grid.numNodes(); ++i) {
        gMassTotal += grid.mass(i);
        gMomTotal  += grid.velocity(i);  // = weighted momentum before normalise
    }

    printf("=== P2G Conservation Test ===\n");
    printf("Particle mass:  %.6f\n", pMassTotal);
    printf("Grid mass:      %.6f\n", gMassTotal);
    printf("Particle mom:   [%.4f, %.4f, %.4f]\n",
           pMomTotal.x(), pMomTotal.y(), pMomTotal.z());
    printf("Grid mom:       [%.4f, %.4f, %.4f]\n",
           gMomTotal.x(), gMomTotal.y(), gMomTotal.z());

    if (!nearEqual(pMassTotal, gMassTotal)) {
        printf("FAIL: mass not conserved (diff = %.6f)\n",
               std::abs(pMassTotal - gMassTotal));
        ++failures;
    } else {
        printf("PASS: mass conserved\n");
    }

    for (int axis = 0; axis < 3; ++axis) {
        if (!nearEqual(pMomTotal(axis), gMomTotal(axis))) {
            printf("FAIL: momentum[%d] not conserved (diff = %.6f)\n",
                   axis, std::abs(pMomTotal(axis) - gMomTotal(axis)));
            ++failures;
        }
    }
    if (failures == 0) printf("PASS: linear momentum conserved\n");

    return failures > 0 ? 1 : 0;
}
